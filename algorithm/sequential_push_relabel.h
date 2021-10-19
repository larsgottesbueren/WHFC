#pragma once

#include <vector>
#include <queue>

#include <tbb/scalable_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/buffered_vector.h"


namespace whfc {

template<typename T>
using vec = std::vector<T, tbb::scalable_allocator<T> >;

class SequentialPushRelabel {
public:
	using Type = SequentialPushRelabel;
	static constexpr bool log = false;
	static constexpr bool capacitate_incoming_edges_of_in_nodes = true;

	explicit SequentialPushRelabel(FlowHypergraph& hg) : hg(hg) { }

	TimeReporter timer;

	Flow computeFlow(Node s, Node t) {
		source = s; target = t;
		clearDatastructures();
		timer.start("push relabel");
		saturateSourceEdges();
		size_t work = 0;
		while (!active.empty()) {
			if (work > global_relabel_work_threshold) {
				work = 0;
				globalRelabel();
			}
			const Node u = active.front();
			active.pop();
			if (excess[u] == 0 || level[u] >= max_level) { continue; }
			if (isHypernode(u)) { work += dischargeHypernode(u); }
			else if (isOutNode(u)) { work += dischargeOutNode(u); }
			else { work += dischargeInNode(u); }
		}
		timer.stop("push relabel");
		return excess[target];
	}

	size_t dischargeHypernode(Node u) {
		size_t work = 0;
		Flow my_excess = excess[u];
		int my_level = level[u];

		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;

			auto i = hg.beginIndexHyperedges(u);
			// push to in-nodes of incident nets
			for ( ; my_excess > 0 && i < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_in = edgeToInNode(e);
				Flow d = my_excess;
				if constexpr (capacitate_incoming_edges_of_in_nodes) {
					// (u, e_in) has infinite capacity but it never makes sense to push more flow into e_in than can be sent on (e_in, e_out)
					d = std::min(d, hg.capacity(e) - flow[inNodeIncidenceIndex(i)]);
				}
				if (my_level == level[e_in] + 1) {
					if (d > 0) {
						if (excess[e_in] == 0) { active.push(e_in); }
						flow[inNodeIncidenceIndex(i)] += d;
						my_excess -= d;
						excess[e_in] += d;
					}
				} else if (my_level <= level[e_in] && d > 0) {
					new_level = std::min(new_level, level[e_in]);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0) {
				break;
			}

			// push back to out-nodes of incident nets
			for (i = hg.beginIndexHyperedges(u); my_excess > 0 && i < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_out = edgeToOutNode(e);
				if (my_level == level[e_out] + 1) {
					assert(flow[outNodeIncidenceIndex(i)] <= hg.capacity(e));
					const Flow d = std::min(my_excess, flow[outNodeIncidenceIndex(i)]);
					if (d > 0) {
						if (excess[e_out] == 0) { active.push(e_out); }
						flow[outNodeIncidenceIndex(i)] -= d;
						my_excess -= d;
						excess[e_out] += d;
					}
				} else if (my_level <= level[e_out] && flow[outNodeIncidenceIndex(i)] > 0) {
					new_level = std::min(new_level, level[e_out]);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0) {
				break;
			}
			my_level = new_level + 1;	// relabel
		}

		level[u] = my_level;	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
			active.push(u);
		}
		excess[u] = my_excess;

		return work;
	}

	size_t dischargeInNode(Node e_in) {
		size_t work = 0;
		Flow my_excess = excess[e_in];
		int my_level = level[e_in];
		Hyperedge e = inNodeToEdge(e_in); assert(e < hg.numHyperedges());
		Node e_out = edgeToOutNode(e);

		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;

			// push through bridge edge
			if (my_level == level[e_out] + 1) {
				Flow d = std::min(hg.capacity(e) - flow[bridgeEdgeIndex(e)], my_excess);
				if (d > 0) {
					if (excess[e_out] == 0) { active.push(e_out); }
					flow[bridgeEdgeIndex(e)] += d;
					my_excess -= d;
					excess[e_out] += d;
				}
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
				if constexpr (capacitate_incoming_edges_of_in_nodes) {
					assert(d <= hg.capacity(e));
				}
				if (my_level == level[v] + 1) {
					if (d > 0) {
						d = std::min(d, my_excess);
						flow[j] -= d;
						my_excess -= d;
						if (v != target && excess[v] == 0) { active.push(v); }
						excess[v] += d;
					}
				} else if (my_level <= level[v] && d > 0) {
					new_level = std::min(new_level, level[v]);
				}
			}
			work += hg.pinCount(e) + 6;

			if (my_excess == 0) {
				break;
			}
			my_level = new_level + 1; 	// relabel
		}

		level[e_in] = my_level; 	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
			active.push(e_in);
		}
		excess[e_in] = my_excess;
		return work;
	}

	size_t dischargeOutNode(Node e_out) {
		size_t work = 0;
		Flow my_excess = excess[e_out];
		int my_level = level[e_out];
		Hyperedge e = outNodeToEdge(e_out); assert(e < hg.numHyperedges());
		Node e_in = edgeToInNode(e);
		assert(my_excess <= hg.capacity(e));

		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;

			// push out to pins
			for (const auto& p : hg.pinsOf(e)) {
				if (my_excess == 0) {
					break;
				}
				Node v = p.pin;
				Flow d = my_excess;
				if (my_level == level[v] + 1) {
					assert(d <= hg.capacity(e) - flow[outNodeIncidenceIndex(p.he_inc_iter)]);
					flow[outNodeIncidenceIndex(p.he_inc_iter)] += d;
					my_excess -= d;
					if (v != target && excess[v] == 0) { active.push(v); }
					excess[v] += d;
				} else if (my_level <= level[v] && d > 0) {
					new_level = std::min(new_level, level[v]);
				}
			}
			work += hg.pinCount(e) + 6;

			if (my_excess == 0) {
				break;
			}

			// push back through bridge edge
			if (my_level == level[e_in] + 1) {
				Flow d = std::min(flow[bridgeEdgeIndex(e)], my_excess);
				if (d > 0) {
					flow[bridgeEdgeIndex(e)] -= d;
					my_excess -= d;
					if (excess[e_in] == 0) { active.push(e_in); }
					excess[e_in] += d;
				}
			} else if (my_level <= level[e_in] && flow[bridgeEdgeIndex(e)] > 0) {
				new_level = std::min(new_level, level[e_in]);
			}

			if (my_excess == 0) {
				break;
			}
			my_level = new_level + 1; 	// relabel
		}

		level[e_out] = my_level; 	// make relabel visible
		if (my_level < max_level && my_excess > 0) {	// go again in the next round if excess left
			active.push(e_out);
		}
		excess[e_out] = my_excess;
		return work;
	}

	void globalRelabel() {
		level.assign(max_level, max_level);
		relabel_queue.clear();
		relabel_queue.push_back(target);
		level[target] = 0;
		int dist = 1;
		size_t first = 0, last = 1;
		while (first != last) {
			for (size_t i = first; i < last; ++i) {
				auto push = [&](const Node v) {
					if (level[v] == max_level) {
						level[v] = dist;
						relabel_queue.push_back(v);
					}
				};
				const Node u = relabel_queue[i];
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
			}
			first = last;
			last = relabel_queue.size();
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
			active.push(edgeToInNode(e));
		}
	}

	void clearDatastructures() {
		out_node_offset = hg.numPins();
		bridge_node_offset = 2 * hg.numPins();

		max_level = hg.numNodes() + 2 * hg.numHyperedges();

		relabel_queue.reserve(max_level);

		flow.assign(2 * hg.numPins() + hg.numHyperedges(), 0);
		excess.assign(max_level, 0);
		level.assign(max_level, 0);

		work_since_last_global_relabel = std::numeric_limits<size_t>::max();
		global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;
	}

private:
	int max_level = 0;
	FlowHypergraph& hg;
	vec<Flow> flow, excess;
	vec<int> level;
	std::queue<Node> active;
	vec<Node> relabel_queue;

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
};

}
