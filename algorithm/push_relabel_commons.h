#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/distance_reachable_sets.h"
#include "../datastructure/queue.h"

#include <tbb/scalable_allocator.h>

namespace whfc {
	template<typename T>
	using vec = std::vector<T, tbb::scalable_allocator<T> >;

	class PushRelabelCommons {
	public:
		PushRelabelCommons(FlowHypergraph& hg) : hg(hg) {

		}

		FlowHypergraph& hg;
		TimeReporter timer;
		Flow upper_flow_bound = std::numeric_limits<Flow>::max();
		bool shall_terminate = false;
		bool has_cut = false;

		/** mapping between ID types */
		// hypernodes | in-nodes | out-nodes
		// [0..n - 1][n..n+m-1][n+m..n+2m]
		bool isHypernode(Node u) const { return u < hg.numNodes(); }
		bool isInNode(Node u) const { return u >= hg.numNodes() && u < hg.numNodes() + hg.numHyperedges(); }
		bool isOutNode(Node u) const { assert(u < hg.numNodes() + 2 * hg.numHyperedges()); return u >= hg.numNodes() + hg.numHyperedges(); }
		Hyperedge inNodeToEdge(Node u) const { assert(isInNode(u)); return Hyperedge(u - hg.numNodes()); }
		Hyperedge outNodeToEdge(Node u) const { assert(isOutNode(u)); return Hyperedge(u - hg.numNodes() - hg.numHyperedges()); }
		Node edgeToInNode(Hyperedge e) const { assert(e < hg.numHyperedges()); return Node(e + hg.numNodes()); }
		Node edgeToOutNode(Hyperedge e) const { assert(e < hg.numHyperedges()); return Node(e + hg.numNodes() + hg.numHyperedges()); }

		/** flow assignment */
		Flow flow_value = 0;
		vec<Flow> flow;
		vec<Flow> excess;
		size_t out_node_offset = 0, bridge_node_offset = 0;

		// position where flow going from vertex into hyperedge is stored
		size_t inNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind; }
		// position where flow going from hyperedge into vertex is stored
		size_t outNodeIncidenceIndex(InHeIndex inc_he_ind) const { return inc_he_ind + out_node_offset; }
		// position where flow on hyperedge is stored
		size_t bridgeEdgeIndex(Hyperedge he) const { return he + bridge_node_offset; }


		/** levels */
		int max_level = 0;
		vec<int> level;
		// to avoid concurrently pushing the same edge in different directions
		bool winEdge(Node v, Node u) { return level[u] == level[v] + 1 || level[u] < level[v] - 1 || (level[u] == level[v] && u < v); }

		/** reachability */
		vec<uint32_t> reach;
		uint32_t source_reachable_stamp = 0, target_reachable_stamp = 0, running_timestamp = 0;
		bool isSource(Node u) const { return reach[u] == 1; }
		void makeSource(Node u) { reach[u] = 1; }
		bool isSourceReachable(Node u) const { return isSource(u) || reach[u] == source_reachable_stamp; }
		void reachFromSource(Node u) { reach[u] = source_reachable_stamp; }
		bool isTarget(Node u) const { return reach[u] == 2; }
		void makeTarget(Node u) { reach[u] = 2; }
		bool isTargetReachable(Node u) const { return isTarget(u) || reach[u] == target_reachable_stamp; }
		void reachFromTarget(Node u) { reach[u] = target_reachable_stamp; }
		void unreach(Node u) { reach[u] = 0; }
		void resetReachability(bool forward) {
			if (++running_timestamp == 0) {
				reach.assign(max_level, 0);
				running_timestamp = 3;	// initial value is 2, but the increment is already executed
			}
			if (forward) {
				source_reachable_stamp = running_timestamp;
			} else {
				target_reachable_stamp = running_timestamp;
			}
		}

		/** global relabeling */
		static constexpr size_t global_relabel_alpha = 6;
		static constexpr size_t global_relabel_frequency = 5;
		size_t work_since_last_global_relabel = 0, global_relabel_work_threshold = 0;

		/** source / sink */
		vec<Node> source_piercing_nodes, target_piercing_nodes;

		void reset() {
			out_node_offset = hg.numPins();
			bridge_node_offset = 2 * hg.numPins();

			max_level = hg.numNodes() + 2 * hg.numHyperedges();

			flow_value = 0;
			flow.assign(2 * hg.numPins() + hg.numHyperedges(), 0);
			excess.assign(max_level, 0);
			level.assign(max_level, 0);

			reach.assign(max_level, 0);
			running_timestamp = 2;

			work_since_last_global_relabel = std::numeric_limits<size_t>::max();
			global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;

			shall_terminate = false;
		}

		/** BFS stuff */
		template<typename PushFunc>
		void scanBackward(Node u, PushFunc&& push) {
			if (isHypernode(u)) {
				for (InHeIndex incnet_ind : hg.incidentHyperedgeIndices(u)) {
					const Hyperedge e = hg.getInHe(incnet_ind).e;
					if (flow[inNodeIncidenceIndex(incnet_ind)] > 0) {
						push(edgeToInNode(e));
					}
					push(edgeToOutNode(e));
				}
			} else if (isOutNode(u)) {
				const Hyperedge e = outNodeToEdge(u);
				if (flow[bridgeEdgeIndex(e)] < hg.capacity(e)) {
					push(edgeToInNode(e));
				}
				for (const auto& pin : hg.pinsOf(e)) {
					if (flow[outNodeIncidenceIndex(pin.he_inc_iter)] > 0) {
						push(pin.pin);
					}
				}
			} else {
				assert(isInNode(u));
				const Hyperedge e = inNodeToEdge(u);
				if (flow[bridgeEdgeIndex(e)] > 0) {
					push(edgeToOutNode(e));
				}
				for (const auto& pin : hg.pinsOf(e)) {
					if (flow[inNodeIncidenceIndex(pin.he_inc_iter)] < hg.capacity(e)) {
						push(pin.pin);
					}
				}
			}
		}

		template<typename PushFunc>
		void scanForward(Node u, PushFunc&& push) {
			if (isHypernode(u)) {
				for (InHeIndex incnet_ind : hg.incidentHyperedgeIndices(u)) {
					const Hyperedge e = hg.getInHe(incnet_ind).e;
					if (flow[inNodeIncidenceIndex(incnet_ind)] < hg.capacity(e) || excess[edgeToInNode(e)] > 0) {
						push(edgeToInNode(e));
					}
					if (flow[outNodeIncidenceIndex(incnet_ind)] > 0 || excess[edgeToOutNode(e)] > 0) {
						push(edgeToOutNode(e));
					}
				}
			} else if (isOutNode(u)) {
				const Hyperedge e = outNodeToEdge(u);
				if (flow[bridgeEdgeIndex(e)] > 0 || excess[edgeToInNode(e)] > 0) {
					push(edgeToInNode(e));
				}
				for (const auto& pin : hg.pinsOf(e)) {
					push(pin.pin);
				}
			} else {
				assert(isInNode(u));
				const Hyperedge e = inNodeToEdge(u);
				if (flow[bridgeEdgeIndex(e)] < hg.capacity(e) || excess[edgeToOutNode(e)] > 0) {
					push(edgeToOutNode(e));
				}
				for (const auto& pin : hg.pinsOf(e)) {
					if (flow[inNodeIncidenceIndex(pin.he_inc_iter)] > 0 || excess[pin.pin] > 0) {
						push(pin.pin);
					}
				}
			}
		}

	};
}
