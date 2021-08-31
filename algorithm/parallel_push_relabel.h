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

	explicit ParallelPushRelabel(FlowHypergraph& hg) : hg(hg), next_active(0) { }


	Flow computeFlow(Node source, Node target) {
		clearDatastructures();
		work_since_last_global_relabel = std::numeric_limits<size_t>::max();
		saturateSourceEdges(source);
		while (true) {
			if (work_since_last_global_relabel > global_relabel_work_threshold) {
				globalRelabel(source, target);
				work_since_last_global_relabel = 0;
			}
			if (next_active.empty()) {
				break;
			}
			dischargeActiveNodes();
		}
		flowDecomposition();
		return excess[target];
	}

	void dischargeActiveNodes() {

	}

	void globalRelabel(Node source, Node target) {
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
					if (level[v] != max_level && __atomic_exchange_n(&level[v], dist, __ATOMIC_ACQ_REL) == max_level) {
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

	void saturateSourceEdges(Node source) {
		level[source] = max_level;
		for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(source)) {
			const Hyperedge e = hg.getInHe(inc_iter).e;
			const Flow d = hg.capacity(e);
			excess[source] -= d;
			excess[edgeToInNode(e)] += d;
			flow[inNodeIncidenceIndex(inc_iter)] += d;
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
		global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;
	}

private:
	int max_level = 0;
	FlowHypergraph& hg;
	vec<Flow> flow, excess;
	vec<int> level;
	BufferedVector<Node> next_active;
	vec<Node> active;

	static constexpr size_t global_relabel_alpha = 6;
	static constexpr size_t global_relabel_frequency = 5;
	size_t work_since_last_global_relabel = 0, global_relabel_work_threshold = 0;

	size_t out_node_offset = 0, bridge_node_offset = 0;

	size_t inNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind; }
	size_t outNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind + out_node_offset; }
	size_t bridgeNodeIndex(Hyperedge he) const { return he + 2 * bridge_node_offset; }

	bool isHypernode(Node u) const { return u < hg.numNodes(); }
	bool isInNode(Node u) const { return u >= hg.numNodes() && u < hg.numNodes() + hg.numHyperedges(); }
	bool isOutNode(Node u) const { return u >= hg.numNodes() + hg.numHyperedges(); }
	Hyperedge inNodeToEdge(Node u) const { return Hyperedge(u - hg.numNodes()); }
	Hyperedge outNodeToEdge(Node u) const { return Hyperedge(u - hg.numNodes() - hg.numHyperedges()); }
	Node edgeToInNode(Hyperedge e) const { return Node(e + hg.numNodes()); }
	Node edgeToOutNode(Hyperedge e) const { return Node(e + hg.numNodes() + hg.numHyperedges()); }
};

}
