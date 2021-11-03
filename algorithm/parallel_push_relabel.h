#pragma once

#include "push_relabel_commons.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include "../datastructure/buffered_vector.h"

#include "cutter_state.h"

namespace whfc {

class ParallelPushRelabel : public PushRelabelCommons {
public:
	using Type = ParallelPushRelabel;
	static constexpr bool log = false;
	static constexpr bool capacitate_incoming_edges_of_in_nodes = true;

	ParallelPushRelabel(FlowHypergraph& hg) : PushRelabelCommons(hg), next_active(0) { }

	bool findMinCuts() {
		saturateSourceEdges();
		size_t num_tries = 0;
		do {
			while (!next_active.empty()) {
				if (flow_value > upper_flow_bound || shall_terminate) {
					return false;
				}
				num_active = next_active.size();
				next_active.swap_container(active);
				if (work_since_last_global_relabel > 2 * global_relabel_work_threshold) {
					globalRelabel();
				}
				dischargeActiveNodes();
				applyUpdates();
			}

			// no more nodes with level < n and excess > 0 left.
			// however labels might be broken from parallelism
			// --> run global relabeling to check if done.
			num_active = 0;
			globalRelabel<true>();	// setting the template parameter to true means the function sets reachability info, since we expect to be finished
			// plug queue back in (regular loop picks it out again)
			next_active.swap_container(active);
			next_active.set_size(num_active);
			num_tries++;
		} while (!next_active.empty());

		deriveSourceSideCut();
		return true;
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
			if (level[u] >= max_level || isTarget(u)) { return; }		// target nodes can now be pushed to consume the updates
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
			if (isTarget(u) && isHypernode(u)) {
				__atomic_fetch_add(&flow_value, excess_diff[u], __ATOMIC_RELAXED);
			}
			excess_diff[u] = 0;
		});
		tbb::parallel_for(0UL, next_active.size(), [&](size_t i) {
			const Node u = next_active[i];
			excess[u] += excess_diff[u];
			if (isTarget(u) && isHypernode(u)) {
				__atomic_fetch_add(&flow_value, excess_diff[u], __ATOMIC_RELAXED);
			}
			excess_diff[u] = 0;
		});
	}

	size_t dischargeHypernode(Node u) {
		auto next_active_handle = next_active.local_buffer();
		auto push = [&](Node v) { if (activate(v)) next_active_handle.push_back(v); };
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
		auto push = [&](Node v) { if (activate(v)) next_active_handle.push_back(v); };
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
		auto push = [&](Node v) { if (activate(v)) next_active_handle.push_back(v); };
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

	template<bool set_reachability = false>
	void globalRelabel() {
		work_since_last_global_relabel = 0;
		tbb::parallel_for(0, max_level, [&](size_t i) { level[i] = max_level; }, tbb::static_partitioner());
		next_active.clear();
		for (const Node& s : target_piercing_nodes) {
			level[s] = 0;
			next_active.push_back_atomic(s);
		}

		if (set_reachability) {
			resetReachability(false);
		}

		auto scan = [&](Node u, int dist) {
			auto next_layer = next_active.local_buffer();
			scanBackward(u, [&](const Node v) {
				if (!isSource(v) && !isTarget(v) && level[v] == max_level
						&& __atomic_exchange_n(&level[v], dist, __ATOMIC_ACQ_REL) == max_level) {
					next_layer.push_back(v);
					if constexpr (set_reachability) {
						reach[v] = target_reachable_stamp;
					}
				}
			});

			if (excess[u] > 0 && last_activated[u] != round) { // add previously mis-labeled nodes to active queue, if not already contained
				size_t pos = __atomic_fetch_add(&num_active, 1, __ATOMIC_RELAXED);
				active[pos] = u;
			}
		};

		parallelBFS(0, scan);

		if (set_reachability) {
			last_target_side_queue_entry = next_active.size();
		}
	}

	void deriveSourceSideCut() {
		// after global relabel with termination check its container is swapped out --> this function doesn't swap
		resetReachability(true);
		for (const Node& s : source_piercing_nodes) {
			next_active.push_back_atomic(s);
		}

		auto scan = [&](Node u, int ) {
			scanForward(u, [&](const Node v) {
				auto next_layer = next_active.local_buffer();
				if (!isSourceReachable(v) && __atomic_exchange_n(&reach[v], source_reachable_stamp, __ATOMIC_ACQ_REL)) {	/* TODO get the right function! */
					next_layer.push_back(v);
				}
			});
		};

		parallelBFS(0, scan);
		last_source_side_queue_entry = next_active.size();
	}

	void deriveTargetSideCut() {
		next_active.swap_container(active);		// don't overwrite contents of source side
		next_active.clear();

		resetReachability(false);
		for (const Node& t : target_piercing_nodes) {
			next_active.push_back_atomic(t);
		}

		auto scan = [&](Node u, int ) {
			scanForward(u, [&](const Node v) {
				auto next_layer = next_active.local_buffer();
				if (!isSourceReachable(v) && __atomic_exchange_n(&reach[v], source_reachable_stamp, __ATOMIC_ACQ_REL)) {	/* TODO get the right function! */
					next_layer.push_back(v);
				}
			});
		};

		parallelBFS(0, scan);

		last_target_side_queue_entry = next_active.size();
		next_active.swap_container(active); 	// go back
	}

	template<typename ScanFunc>
	void parallelBFS(size_t first, ScanFunc&& scan) {
		size_t last = next_active.size();
		int dist = 1;
		while (first != last) {
			tbb::parallel_for(first, last, [&](size_t i) { scan(next_active[i], dist); });
			next_active.finalize();
			first = last;
			last = next_active.size();
			dist++;
		}
	}

	std::pair<NodeWeight, NodeWeight> computeReachableWeights() {
		// next_active container is for source side. active for target side
		NodeWeight extra_source_weight = tbb::parallel_reduce(
				tbb::blocked_range<size_t>(0, last_source_side_queue_entry), 0, [&](const auto& r, NodeWeight sum) -> NodeWeight {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Node u = next_active[i];
				if (isHypernode(u) && !isSource(u)) {
					sum += hg.nodeWeight(u);
				}
			}
			return sum;
		}, std::plus<>());

		NodeWeight extra_target_weight = tbb::parallel_reduce(
				tbb::blocked_range<size_t>(0, last_target_side_queue_entry), 0, [&](const auto& r, NodeWeight sum) -> NodeWeight {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Node u = active[i];
				if (isHypernode(u) && !isTarget(u)) {
					sum += hg.nodeWeight(u);
				}
			}
			return sum;
		}, std::plus<>());

		return std::make_pair(extra_source_weight, extra_target_weight);
	}

	void assimilate(CutterState<Type>& cs, bool source_side) {
		// TODO if target side is assimilated, the distance labels break. we can just keep running and in case it breaks the termination check global relabel strikes
		if (source_side) {
			for (size_t i = 0; i < last_source_side_queue_entry; ++i) {
				Node u = next_active[i];
				if (!isSource(u) && isInNode(u)) {
					Hyperedge e = inNodeToEdge(u);
					Node out_node = edgeToOutNode(e);
					if (!isSourceReachable(out_node)) {		// in node visited but not out node --> cut hyperedge
						cs.addToSourceSideCut(e);
					}
				}
				makeSource(u);
			}
		} else {
			for (size_t i = 0; i < last_target_side_queue_entry; ++i) {
				Node u = active[i];
				if (!isTarget(u) && isOutNode(u)) {
					Hyperedge e = outNodeToEdge(u);
					Node in_node = edgeToInNode(e);
					if (!isTargetReachable(in_node)) {		// out node visited but not in node --> cut hyperedge
						cs.addToTargetSideCut(e);
					}
				}
				makeTarget(u);
			}
		}
	}

	void saturateSourceEdges() {
		// TODO parallelize?
		next_active.clear();
		for (const Node& source : source_piercing_nodes) {
			level[source] = max_level;
			for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(source)) {
				const Hyperedge e = hg.getInHe(inc_iter).e;
				const Flow d = hg.capacity(e);		// TODO adapt for residual capacity
				excess[source] -= d;
				excess[edgeToInNode(e)] += d;
				flow[inNodeIncidenceIndex(inc_iter)] += d;
				next_active.push_back_atomic(edgeToInNode(e));
			}
		}
	}

	void reset() {
		PushRelabelCommons::reset();

		excess_diff.assign(max_level, 0);
		next_level.assign(max_level, 0);

		next_active.clear();
		next_active.adapt_capacity(max_level);
		active.resize(max_level);
		last_activated.assign(max_level, 0);
		round = 0;

		last_source_side_queue_entry = 0;
		last_target_side_queue_entry = 0;
	}


private:
	vec<Flow> excess_diff;
	vec<int> next_level;
	size_t num_active = 0;
	BufferedVector<Node> next_active;
	vec<Node> active;

	vec<uint32_t> last_activated;
	uint32_t round = 0;
	bool activate(Node u) {
		return last_activated[u] != round && __atomic_exchange_n(&last_activated[u], round, __ATOMIC_ACQ_REL) != round;
	}

	size_t last_source_side_queue_entry = 0, last_target_side_queue_entry = 0;
};

}
