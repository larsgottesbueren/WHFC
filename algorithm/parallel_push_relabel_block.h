#pragma once

#include <vector>

#include <tbb/scalable_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/buffered_vector.h"


namespace whfc {

template<typename T>
using vec = std::vector<T, tbb::scalable_allocator<T> >;

class ParallelPushRelabel {

  enum class LevelState : uint8_t {
		NOT_MODIFIED,
		EXPECT_STABLE,
		IS_RELABELED
	};

public:
	using Type = ParallelPushRelabel;
	static constexpr bool log = false;
	static constexpr bool capacitate_incoming_edges_of_in_nodes = true;

	explicit ParallelPushRelabel(FlowHypergraph& hg) : hg(hg), next_active(0) { }

	TimeReporter timer;

	Flow computeFlow(Node s, Node t) {
		source = s; target = t;
		clearDatastructures();
		timer.start("push relabel");
		saturateSourceEdges();
		while (!next_active.empty()) {
			size_t num_active = next_active.size();
			next_active.swap_container(active);
			if (work_since_last_global_relabel > 2 * global_relabel_work_threshold) {
				globalRelabel();
				work_since_last_global_relabel = 0;
			}
			dischargeActiveNodes(num_active);
			applyUpdates(num_active);
		}
		excess[target] += excess_diff[target];
		excess_diff[target] = 0;
		timer.stop("push relabel");
		return excess[target];
	}

	void dischargeActiveNodes(size_t num_active) {
		if (++round == 0) {
			last_activated.assign(max_level, 0);
			++round;
		}
		next_active.clear();
		tbb::enumerable_thread_specific<size_t> work;
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

	void applyUpdates(size_t num_active) {
		tbb::parallel_for(0UL, num_active, [&](size_t i) {
			const Node u = active[i];
			if (level[u] >= max_level) { assert(excess_diff[u] == 0); return; }
			level[u] = next_level[u];
			excess[u] += excess_diff[u];
			excess_diff[u] = 0;
			node_state[u] = LevelState::NOT_MODIFIED;
		});
		tbb::parallel_for(0UL, next_active.size(), [&](size_t i) {
			const Node u = next_active[i];
			assert(node_state[u] == LevelState::NOT_MODIFIED);
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
				if (my_level == level[e_in] + 1 && d > 0) {
					if ( excess[e_in] == 0 || updateNodeState(e_in, LevelState::EXPECT_STABLE) ) {
						flow[inNodeIncidenceIndex(i)] += d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_in], d, __ATOMIC_RELAXED);
						push(e_in);
					} else {
						skipped = true;
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
				const Flow d = std::min(my_excess, flow[outNodeIncidenceIndex(i)]);
				if (my_level == level[e_out] + 1 && d > 0) {
					if ( excess[e_out] == 0 || updateNodeState(e_out, LevelState::EXPECT_STABLE) ) {
						assert(flow[outNodeIncidenceIndex(i)] <= hg.capacity(e));
						flow[outNodeIncidenceIndex(i)] -= d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_out], d, __ATOMIC_RELAXED);
						push(e_out);
					} else {
						skipped = true;
					}
				} else if (my_level <= level[e_out] && flow[outNodeIncidenceIndex(i)] > 0) {
					new_level = std::min(new_level, level[e_out]);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0 || skipped) {
				break;
			}

			if ( updateNodeState(u, LevelState::IS_RELABELED) ) {
				my_level = new_level + 1;	// relabel
			} else {
				break;
			}
		}

		next_level[u] = my_level;	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
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
			const Flow d = std::min(hg.capacity(e) - flow[bridgeEdgeIndex(e)], my_excess);
			if (my_level == level[e_out] + 1 && d > 0) {
				if ( excess[e_out] == 0 || updateNodeState(e_out, LevelState::EXPECT_STABLE) ) {
					flow[bridgeEdgeIndex(e)] += d;
					my_excess -= d;
					__atomic_fetch_add(&excess_diff[e_out], d, __ATOMIC_RELAXED);
					push(e_out);
				} else {
					skipped = true;
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
				if (my_level == level[v] + 1 && d > 0) {
					if ( excess[v] == 0 || updateNodeState(v, LevelState::EXPECT_STABLE) ) {
						d = std::min(d, my_excess);
						flow[j] -= d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[v], d, __ATOMIC_RELAXED);
						push(v);
					} else {
						skipped = true;
					}
				} else if (d > 0 && my_level <= level[v]) {
					new_level = std::min(new_level, level[v]);
				}
				work++;
			}

			if (my_excess == 0 || skipped) {
				break;
			}

			if ( updateNodeState(e_in, LevelState::IS_RELABELED) ) {
				my_level = new_level + 1;	// relabel
			} else {
				break;
			}
		}

		next_level[e_in] = my_level; 	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
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
					if ( excess[v] == 0 || updateNodeState(v, LevelState::EXPECT_STABLE) ) {
						assert(d > 0 && d <= hg.capacity(e) - flow[outNodeIncidenceIndex(p.he_inc_iter)]);
						flow[outNodeIncidenceIndex(p.he_inc_iter)] += d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[v], d, __ATOMIC_RELAXED);
						push(v);
					} else {
						skipped = true;
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
			Flow d = std::min(flow[bridgeEdgeIndex(e)], my_excess);
			if (my_level == level[e_in] + 1 && d > 0) {
				if ( excess[e_in] == 0 || updateNodeState(e_in, LevelState::EXPECT_STABLE) ) {
					flow[bridgeEdgeIndex(e)] -= d;
					my_excess -= d;
					__atomic_fetch_add(&excess_diff[e_in], d, __ATOMIC_RELAXED);
					push(e_in);
					work++;
				} else {
					skipped = true;
				}
			} else if (my_level <= level[e_in] && flow[bridgeEdgeIndex(e)] > 0) {
				new_level = std::min(new_level, level[e_in]);
			}

			if (my_excess == 0 || skipped) {
				break;
			}

			if ( updateNodeState(e_out, LevelState::IS_RELABELED) ) {
				my_level = new_level + 1;	// relabel
			} else {
				break;
			}
		}

		next_level[e_out] = my_level; 	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
			push(e_out);
		}
		__atomic_fetch_sub(&excess_diff[e_out], (excess[e_out] - my_excess), __ATOMIC_RELAXED); // excess[u] serves as indicator for other nodes that u is active --> update later
		return work;
	}

	void globalRelabel() {
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
		out_node_offset = hg.numPins();
		bridge_node_offset = 2 * hg.numPins();

		max_level = hg.numNodes() + 2 * hg.numHyperedges();

		flow.assign(2 * hg.numPins() + hg.numHyperedges(), 0);
		excess.assign(max_level, 0);
		excess_diff.assign(max_level, 0);
		level.assign(max_level, 0);
		next_level.assign(max_level, 0);

		next_active.adapt_capacity(max_level);
		active.resize(max_level);
		node_state.assign(max_level, LevelState::NOT_MODIFIED);
		last_activated.assign(max_level, 0);
		round = 0;

		work_since_last_global_relabel = std::numeric_limits<size_t>::max();
		global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;
	}

private:
	int max_level = 0;
	FlowHypergraph& hg;
	vec<Flow> flow, excess, excess_diff;
	vec<int> level, next_level;
	BufferedVector<Node> next_active;
	vec<Node> active;
	vec<LevelState> node_state;

	Node source, target;

	static constexpr size_t global_relabel_alpha = 6;
	static constexpr size_t global_relabel_frequency = 5;
	size_t work_since_last_global_relabel = 0, global_relabel_work_threshold = 0;

	size_t out_node_offset = 0, bridge_node_offset = 0;

	// position where flow going from vertex into hyperedge is stored
	size_t inNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind; }
	// position where flow going from hyperedge into vertex is stored
	size_t outNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind + out_node_offset; }
	// position where flow on hyperedge is stored
	size_t bridgeEdgeIndex(Hyperedge he) const { return he + bridge_node_offset; }

	// hypernodes | in-nodes | out-nodes
	// [0..n - 1][n..n+m-1][n+m..n+2m]
	bool isHypernode(Node u) const { return u < hg.numNodes(); }
	bool isInNode(Node u) const { return u >= hg.numNodes() && u < hg.numNodes() + hg.numHyperedges(); }
	bool isOutNode(Node u) const { assert(u < hg.numNodes() + 2 * hg.numHyperedges()); return u >= hg.numNodes() + hg.numHyperedges(); }
	Hyperedge inNodeToEdge(Node u) const { assert(isInNode(u)); return Hyperedge(u - hg.numNodes()); }
	Hyperedge outNodeToEdge(Node u) const { assert(isOutNode(u)); return Hyperedge(u - hg.numNodes() - hg.numHyperedges()); }
	Node edgeToInNode(Hyperedge e) const { assert(e < hg.numHyperedges()); return Node(e + hg.numNodes()); }
	Node edgeToOutNode(Hyperedge e) const { assert(e < hg.numHyperedges()); return Node(e + hg.numNodes() + hg.numHyperedges()); }

	bool winEdge(Node u, Node v) {
		return level[u] == level[v] + 1 || level[u] < level[v] - 1 || (level[u] == level[v] && u < v);
	}

	bool updateNodeState(const Node u, const LevelState desired) {
		LevelState expected = LevelState::NOT_MODIFIED;
		return node_state[u] == desired ||
					 __atomic_compare_exchange_n(&node_state[u], &expected,
						static_cast<uint8_t>(desired), false,
						__ATOMIC_ACQ_REL, __ATOMIC_RELAXED);
	}

	vec<uint32_t> last_activated;
	uint32_t round = 0;
	bool activate(Node u) {
		return last_activated[u] != round && __atomic_exchange_n(&last_activated[u], round, __ATOMIC_ACQ_REL) != round;
	}

};

}
