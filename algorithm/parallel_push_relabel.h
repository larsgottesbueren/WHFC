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
public:
	using Type = ParallelPushRelabel;
	static constexpr bool log = true;

	explicit ParallelPushRelabel(FlowHypergraph& hg) : hg(hg), next_active(0) { }


	Flow computeFlow(Node s, Node t) {
		source = s; target = t;
		clearDatastructures();
		work_since_last_global_relabel = std::numeric_limits<size_t>::max();
		saturateSourceEdges();
		while (!next_active.empty()) {
			size_t num_active = next_active.size();
			next_active.swap_container(active);
			if (work_since_last_global_relabel > global_relabel_work_threshold) {
				globalRelabel();
				work_since_last_global_relabel = 0;
			}
			dischargeActiveNodes(num_active);
			applyUpdates(num_active);
		}
		flowDecomposition();
		return excess[target];
	}

	void dischargeActiveNodes(size_t num_active) {
		if (++round == 0) {
			last_activated.assign(max_level, 0);
		}
		next_active.clear();
		tbb::parallel_for(0UL, num_active, [&](size_t i) {
			const Node u = active[i];
			if (level[u] >= max_level) {
				return;
			}
			if (isHypernode(u)) {
				dischargeHypernode(u);
			} else if (isOutNode(u)) {
				LOGGER << "discharge out node";
			} else {
				LOGGER << "discharge in node";
			}
		});
		next_active.finalize();
	}

	void applyUpdates(size_t num_active) {
		tbb::parallel_for(0UL, num_active, [&](size_t i) {
			const Node u = active[i];
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

	void dischargeHypernode(Node u) {
		auto next_active_handle = next_active.local_buffer();
		auto push = [&](Node v) { if (v != target && activate(v)) next_active_handle.push_back(v); };
		size_t work = 0;

		Flow my_excess = excess[u];
		int my_level = level[u];
		while (my_excess > 0 && my_level < max_level) {
			int new_level = max_level;
			bool skipped = false;

			auto i = hg.beginIndexHyperedges(u);
			for ( ; my_excess > 0 && i < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_in = edgeToInNode(e);
				if (my_level == level[e_in] + 1) {
					if (excess[e_in] > 0 && !winEdge(u, e_in)) {
						skipped = true;
						continue;
					}
					const Flow d = my_excess; 	// inf cap. TODO but does it make sense to push more than the hyperedge can take??? my sequential code has this as well but commented out
					flow[inNodeIncidenceIndex(i)] += d;
					__atomic_fetch_add(&excess_diff[e_in], d, __ATOMIC_RELAXED);
					my_excess -= d;
					push(e_in);
				} else if (my_level <= level[e_in]) {
					new_level = std::min(new_level, level[e_in] + 1);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0) {
				break;
			}

			for (i = hg.beginIndexHyperedges(u); my_excess > 0 && u < hg.endIndexHyperedges(u); ++i) {
				Hyperedge e = hg.getInHe(i).e; Node e_out = edgeToOutNode(e);
				if (my_level == level[e_out] + 1) {
					if (excess[e_out] > 0 && !winEdge(u, e_out)) {
						skipped = true;
						continue;
					}
					const Flow d = std::min(my_excess, flow[outNodeIncidenceIndex(i)]);
					if (d > 0) {
						flow[outNodeIncidenceIndex(i)] -= d;
						my_excess -= d;
						__atomic_fetch_add(&excess_diff[e_out], d, __ATOMIC_RELAXED);
						push(e_out);
					}
				} else if (my_level <= level[e_out] && flow[outNodeIncidenceIndex(i)] > 0) {
					new_level = std::min(new_level, level[e_out] + 1);
				}
			}
			work += i - hg.beginIndexHyperedges(u);

			if (my_excess == 0 || skipped) {
				break;
			}

			my_level = new_level;	// relabel
		}

		if (my_level != level[u]) {		// make relabel visible
			next_level[u] = my_level;
		}
		if (my_level < max_level && my_excess < excess[u]) {	// go again in the next round if excess left
			push(u);
		}
		__atomic_fetch_sub(&excess_diff[u], (excess[u] - my_excess), __ATOMIC_RELAXED); // excess[u] serves as indicator for other nodes that u is active --> update later
	}

	void dischargeInNode(Node e_in) {
		Flow my_excess = excess[e_in];
		while (my_excess > 0) {

		}

	}

	void dischargeOutNode(Node e_out) {
		while (excess[e_out] > 0) {

		}

	}

	void globalRelabel() {
		level.assign(max_level, max_level);

		next_active.clear();
		next_active.push_back_atomic(target);	// parallel special case for target/high degree nodes?
		level[target] = 0;
		int dist = 1;
		size_t first = 0, last = 1;
		while (first != last) {
			tbb::parallel_for(first, last, [&](size_t i) {
				auto next_layer = next_active.local_buffer();
				auto push = [&](const Node v) {
					if (level[v] != max_level && __atomic_exchange_n(&level[v], dist, __ATOMIC_ACQ_REL) != max_level) {
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
					if (flow[bridgeNodeIndex(e)] < hg.capacity(e)) { // to e_in if flow(e_in, e_out) < capacity(e)
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
					if (flow[bridgeNodeIndex(e)] > 0) { // to e_out if flow(e_in, e_out) > 0
						push(edgeToOutNode(e));
					}
					for (const auto& pin : hg.pinsOf(e)) { // to v always
						if (pin.pin != source) {
							push(pin.pin);
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

	void flowDecomposition() {

	}

	void clearDatastructures() {
		out_node_offset = hg.numPins();
		bridge_node_offset = 2 * hg.numPins();

		max_level = hg.numNodes() + 2 * hg.numHyperedges();

		flow.assign(2 * hg.numPins() + hg.numHyperedges(), 0);
		excess.assign(max_level, 0);
		level.assign(max_level, 0);

		next_active.adapt_capacity(max_level);
		active.resize(max_level);
		last_activated.assign(max_level, 0);
		round = 0;

		global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;
	}

private:
	int max_level = 0;
	FlowHypergraph& hg;
	vec<Flow> flow, excess, excess_diff;
	vec<int> level, next_level;
	BufferedVector<Node> next_active;
	vec<Node> active;

	Node source, target;

	static constexpr size_t global_relabel_alpha = 6;
	static constexpr size_t global_relabel_frequency = 5;
	size_t work_since_last_global_relabel = 0, global_relabel_work_threshold = 0;

	size_t out_node_offset = 0, bridge_node_offset = 0;

	size_t inNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind; }
	size_t outNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind + out_node_offset; }
	size_t bridgeNodeIndex(Hyperedge he) const { return he + bridge_node_offset; }

	bool isHypernode(Node u) const { return u < hg.numNodes(); }
	bool isInNode(Node u) const { return u >= hg.numNodes() && u < hg.numNodes() + hg.numHyperedges(); }
	bool isOutNode(Node u) const { return u >= hg.numNodes() + hg.numHyperedges(); }
	Hyperedge inNodeToEdge(Node u) const { return Hyperedge(u - hg.numNodes()); }
	Hyperedge outNodeToEdge(Node u) const { return Hyperedge(u - hg.numNodes() - hg.numHyperedges()); }
	Node edgeToInNode(Hyperedge e) const { return Node(e + hg.numNodes()); }
	Node edgeToOutNode(Hyperedge e) const { return Node(e + hg.numNodes() + hg.numHyperedges()); }

	bool winEdge(Node u, Node v) {
		return level[u] == level[v] + 1 || level[u] < level[v] - 1 || (level[u] == level[v] && u < v);
	}

	vec<uint32_t> last_activated;
	uint32_t round = 0;
	bool activate(Node u) {
		return last_activated[u] != round && __atomic_exchange_n(&last_activated[u], round, __ATOMIC_ACQ_REL) != round;
	}

};

}
