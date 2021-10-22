#pragma once

#include "push_relabel_commons.h"

#include <tbb/parallel_for.h>
#include "../datastructure/buffered_vector.h"


namespace whfc {

class ParallelPushRelabel : public PushRelabelCommons {
public:
	using Type = ParallelPushRelabel;
	static constexpr bool log = false;
	static constexpr bool capacitate_incoming_edges_of_in_nodes = true;

	ParallelPushRelabel(FlowHypergraph& hg) : PushRelabelCommons(hg), next_active(0) { }

	Flow computeFlow(Node s, Node t) {
		source = s; target = t;
		clearDatastructures();
		timer.start("push relabel");
		saturateSourceEdges();
		size_t num_tries = 0;
		do {
			while (!next_active.empty()) {
				num_active = next_active.size();
				next_active.swap_container(active);
				if (work_since_last_global_relabel > 2 * global_relabel_work_threshold) {
					globalRelabel();
				}
				dischargeActiveNodes();
				applyUpdates();
			}

			bool high_label_excess_left = false;
			tbb::parallel_for(0, max_level, [&](int i) {
				Node u(i);
				if (u != source && u != target && excess[u] > 0 && level[u] >= max_level) {
					high_label_excess_left = true;
					// can cancel, but we expect that this branch is rarely hit
				}
			});

			// no more nodes with level < n and excess > 0 left.
			// however labels might be broken from parallelism
			// --> run global relabeling to check if done.

			if (high_label_excess_left) {
				num_active = 0;
				globalRelabel();
				// plug queue back in (regular loop picks it out again)
				next_active.swap_container(active);
				next_active.set_size(num_active);
			}
			num_tries++;
		} while (!next_active.empty());

		// target node is never pushed to active set --> apply update separately.
		// TODO once we get multiple target nodes, this may require changes.
		//  the values are only needed to determine the flow diff
		excess[target] += excess_diff[target];
		excess_diff[target] = 0;
		timer.stop("push relabel");
		return excess[target];
	}

	void dischargeActiveNodes() {
		if (++round == 0) {
			last_activated.assign(max_level, 0);
			++round;
		}
		next_active.clear();
		tbb::enumerable_thread_specific<size_t> work(0);
		auto task = [&](size_t i) {
			const Node u = active[i];
			assert(excess[u] > 0);
			if (level[u] >= max_level) { return; }
			if (isHypernode(u)) { work.local() += dischargeHypernode(u); }
			else if (isOutNode(u)) { work.local() += dischargeOutNode(u); }
			else { work.local() += dischargeInNode(u); }
		};
		tbb::parallel_for(0UL, num_active, task);
		next_active.finalize();
		work_since_last_global_relabel += work.combine(std::plus<>());
	}

	void applyUpdates() {
		tbb::parallel_for(0UL, num_active, [&](size_t i) {
			const Node u = active[i];
			if (level[u] >= max_level) { assert(excess_diff[u] == 0); return; }
			level[u] = next_level[u];
			excess[u] += excess_diff[u];
			excess_diff[u] = 0;
		});
		tbb::parallel_for(0UL, next_active.size(), [&](size_t i) {
			const Node u = next_active[i];
			excess[u] += excess_diff[u];
			excess_diff[u] = 0;
		});
	}

	size_t dischargeHypernode(Node u) {
		auto next_active_handle = next_active.local_buffer();
		auto push = [&](Node v) { if (v != target && activate(v)) next_active_handle.push_back(v); };
		size_t work = 0;
		Flow my_excess = excess[u];
		int my_level = level[u];

		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;
			bool skipped = false;

			// push to in-nodes of incident nets
			auto i = hg.beginIndexHyperedges(u);
			for ( ; my_excess > 0 && i < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_in = edgeToInNode(e);
				Flow d = my_excess;
				if constexpr (capacitate_incoming_edges_of_in_nodes) {
					d = std::min(d, hg.capacity(e) - flow[inNodeIncidenceIndex(i)]);
				}
				if (my_level == level[e_in] + 1) {
					if (excess[e_in] > 0 && !winEdge(u, e_in)) {
						skipped = true;
					} else if (d > 0) {
						flow[inNodeIncidenceIndex(i)] += d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_in], d, __ATOMIC_RELAXED);
						push(e_in);
					}
				} else if (my_level <= level[e_in] && d > 0) {
					new_level = std::min(new_level, level[e_in]);
				}
			}
			work += i - hg.beginIndexHyperedges(u);
			if (my_excess == 0)
				break;

			// push back to out-nodes of incident nets
			for (i = hg.beginIndexHyperedges(u); my_excess > 0 && i < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_out = edgeToOutNode(e);
				if (my_level == level[e_out] + 1) {
					if (excess[e_out] > 0 && !winEdge(u, e_out)) {
						skipped = true;
					} else {
						const Flow d = std::min(my_excess, flow[outNodeIncidenceIndex(i)]);
						if (d > 0) {
							assert(flow[outNodeIncidenceIndex(i)] <= hg.capacity(e));
							flow[outNodeIncidenceIndex(i)] -= d;
							my_excess -= d;
							__atomic_fetch_add(&excess_diff[e_out], d, __ATOMIC_RELAXED);
							push(e_out);
						}
					}
				} else if (my_level <= level[e_out] && flow[outNodeIncidenceIndex(i)] > 0) {
					new_level = std::min(new_level, level[e_out]);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0 || skipped) {
				break;
			}
			my_level = new_level + 1;	// relabel
		}

		next_level[u] = my_level;	// make relabel visible
		if (my_excess > 0 && my_level < max_level) {	// go again in the next round if excess left
			push(u);
		}
		__atomic_fetch_sub(&excess_diff[u], (excess[u] - my_excess), __ATOMIC_RELAXED); // excess[u] serves as indicator for other nodes that u is active --> update later
		return work;
	}

	size_t dischargeInNode(Node e_in) {
		auto next_active_handle = next_active.local_buffer();
		auto push = [&](Node v) { if (v != target && activate(v)) next_active_handle.push_back(v); };
		size_t work = 0;
		Flow my_excess = excess[e_in];
		int my_level = level[e_in];
		next_level[e_in] = my_level;
		Hyperedge e = inNodeToEdge(e_in); assert(e < hg.numHyperedges());
		Node e_out = edgeToOutNode(e);
		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;
			bool skipped = false;

			// push through bridge edge
			if (my_level == level[e_out] + 1) {
				if (excess[e_out] > 0 && !winEdge(e_in, e_out)) {
					skipped = true;
				} else {
					const Flow d = std::min(hg.capacity(e) - flow[bridgeEdgeIndex(e)], my_excess);
					if (d > 0) {
						flow[bridgeEdgeIndex(e)] += d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_out], d, __ATOMIC_RELAXED);
						push(e_out);
					}
				}
				work++;
			} else if (my_level <= level[e_out] && flow[bridgeEdgeIndex(e)] < hg.capacity(e)) {
				new_level = std::min(new_level, level[e_out]);
			}

			// push back to pins
			for (const auto& p : hg.pinsOf(e)) {
				if (my_excess == 0) {
					break;
				}
				Node v = p.pin;
				size_t j = inNodeIncidenceIndex(p.he_inc_iter);
				Flow d = flow[j];
				if (my_level == level[v] + 1) {
					if (excess[v] > 0 && !winEdge(e_in, v)) {
						skipped = true;
					} else if (d > 0) {
						d = std::min(d, my_excess);
						flow[j] -= d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[v], d, __ATOMIC_RELAXED);
						push(v);
					}
				} else if (d > 0 && my_level <= level[v]) {
					new_level = std::min(new_level, level[v]);
				}
				work++;
			}

			if (my_excess == 0 || skipped) {
				break;
			}
			my_level = new_level + 1; 	// relabel
		}

		next_level[e_in] = my_level; 	// make relabel visible
		if (my_excess > 0 && my_level < max_level) {	// go again in the next round if excess left
			push(e_in);
		}
		__atomic_fetch_sub(&excess_diff[e_in], (excess[e_in] - my_excess), __ATOMIC_RELAXED); // excess[u] serves as indicator for other nodes that u is active --> update later
		return work;
	}

	size_t dischargeOutNode(Node e_out) {
		auto next_active_handle = next_active.local_buffer();
		auto push = [&](Node v) { if (v != target && activate(v)) next_active_handle.push_back(v); };
		size_t work = 0;
		Flow my_excess = excess[e_out];
		int my_level = level[e_out];
		Hyperedge e = outNodeToEdge(e_out); assert(e < hg.numHyperedges());
		Node e_in = edgeToInNode(e);
		assert(my_excess <= hg.capacity(e));

		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;
			bool skipped = false;

			// push out to pins
			for (const auto& p : hg.pinsOf(e)) {
				if (my_excess == 0) {
					break;
				}
				Node v = p.pin;
				Flow d = my_excess;
				if (my_level == level[v] + 1) {
					if (excess[v] > 0 && !winEdge(e_out, v)) {
						skipped = true;
					} else {
						assert(d > 0 && d <= hg.capacity(e) - flow[outNodeIncidenceIndex(p.he_inc_iter)]);
						flow[outNodeIncidenceIndex(p.he_inc_iter)] += d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[v], d, __ATOMIC_RELAXED);
						push(v);
					}
				} else if (my_level <= level[v]) {
					new_level = std::min(new_level, level[v]);
				}
				work++;
			}

			if (my_excess == 0) {
				break;
			}

			// push back through bridge edge
			if (my_level == level[e_in] + 1) {
				if (excess[e_in] > 0 && !winEdge(e_out, e_in)) {
					skipped = true;
				} else {
					Flow d = std::min(flow[bridgeEdgeIndex(e)], my_excess);
					if (d > 0) {
						flow[bridgeEdgeIndex(e)] -= d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_in], d, __ATOMIC_RELAXED);
						push(e_in);
					}
					work++;
				}
			} else if (my_level <= level[e_in] && flow[bridgeEdgeIndex(e)] > 0) {
				new_level = std::min(new_level, level[e_in]);
			}

			if (my_excess == 0 || skipped) {
				break;
			}
			my_level = new_level + 1; 	// relabel
		}

		next_level[e_out] = my_level; 	// make relabel visible
		if (my_excess > 0 && my_level < max_level) {	// go again in the next round if excess left
			push(e_out);
		}
		__atomic_fetch_sub(&excess_diff[e_out], (excess[e_out] - my_excess), __ATOMIC_RELAXED); // excess[u] serves as indicator for other nodes that u is active --> update later
		return work;
	}

	void globalRelabel() {
		work_since_last_global_relabel = 0;
		tbb::parallel_for(0, max_level, [&](size_t i) { level[i] = max_level; }, tbb::static_partitioner());

		next_active.clear();
		next_active.push_back_atomic(target);	// parallel special case for target/high degree nodes?
		level[target] = 0;
		int dist = 1;
		size_t first = 0, last = 1;
		while (first != last) {
			tbb::parallel_for(first, last, [&](size_t i) {
				auto next_layer = next_active.local_buffer();
				auto push = [&](const Node v) {
					if (level[v] == max_level && __atomic_exchange_n(&level[v], dist, __ATOMIC_ACQ_REL) == max_level) {
						next_layer.push_back(v);
					}
				};
				const Node u = next_active[i];

				// bfs iteration
				if (isHypernode(u)) {
					for (InHeIndex incnet_ind : hg.incidentHyperedgeIndices(u)) {
						const Hyperedge e = hg.getInHe(incnet_ind).e;
						if (flow[inNodeIncidenceIndex(incnet_ind)] > 0) {								/* scan */
							push(edgeToInNode(e));
						}
						push(edgeToOutNode(e));
					}
				} else if (isOutNode(u)) {
					const Hyperedge e = outNodeToEdge(u);
					if (flow[bridgeEdgeIndex(e)] < hg.capacity(e)) { // to e_in if flow(e_in, e_out) < capacity(e)
						push(edgeToInNode(e));
					}
					for (const auto& pin : hg.pinsOf(e)) { // to v if flow(e_out, v) > 0
						if (pin.pin != source && flow[outNodeIncidenceIndex(pin.he_inc_iter)] > 0) {	/* random access */
							push(pin.pin);
						}
					}
				} else {
					assert(isInNode(u));
					const Hyperedge e = inNodeToEdge(u);
					if (flow[bridgeEdgeIndex(e)] > 0) { // to e_out if flow(e_in, e_out) > 0
						push(edgeToOutNode(e));
					}
					for (const auto& pin : hg.pinsOf(e)) { // to v always
						if constexpr (capacitate_incoming_edges_of_in_nodes) {
							if (pin.pin != source && flow[inNodeIncidenceIndex(pin.he_inc_iter)] < hg.capacity(e)) {
								push(pin.pin);
							}
						} else {
							if (pin.pin != source) {
								push(pin.pin);
							}
						}
					}
				}

				// add previously mis-labeled nodes to active queue, if not already contained
				// expected to happen rarely for regular global relabeling
				// bit more frequently for termination checks
				// however even there only few nodes were affected (very little flow missing)
				if (excess[u] > 0 && last_activated[u] != round) {
					size_t pos = __atomic_fetch_add(&num_active, 1, __ATOMIC_RELAXED);
					active[pos] = u;
				}
			});
			next_active.finalize();
			first = last;
			last = next_active.size();
			dist++;
		}
	}

	void saturateSourceEdges() {
		level[source] = max_level;
		for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(source)) {
			const Hyperedge e = hg.getInHe(inc_iter).e;
			const Flow d = hg.capacity(e);
			excess[source] -= d;
			excess[edgeToInNode(e)] += d;
			flow[inNodeIncidenceIndex(inc_iter)] += d;
			next_active.push_back_atomic(edgeToInNode(e));
		}
	}

	void clearDatastructures() {
		PushRelabelCommons::clearDatastructures();

		excess_diff.assign(max_level, 0);
		next_level.assign(max_level, 0);

		next_active.adapt_capacity(max_level);
		active.resize(max_level);
		last_activated.assign(max_level, 0);
		round = 0;
	}


private:
	vec<Flow> excess_diff;
	vec<int> next_level;
	size_t num_active = 0;
	BufferedVector<Node> next_active;
	vec<Node> active;

	Node source, target;

	vec<uint32_t> last_activated;
	uint32_t round = 0;
	bool activate(Node u) {
		return last_activated[u] != round && __atomic_exchange_n(&last_activated[u], round, __ATOMIC_ACQ_REL) != round;
	}
};

}
