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

		using ScanList = LayeredQueue<Node>;
		using ReachableNodes = DistanceReachableNodes;
		using ReachableHyperedges = DistanceReachableHyperedges;

		FlowHypergraph& hg;

		TimeReporter timer;

		int max_level = 0;
		vec<Flow> flow;		// move to cutter state? or make cutter state own flow algorithm?
		vec<Flow> excess;
		vec<int> level;

		static constexpr size_t global_relabel_alpha = 6;
		static constexpr size_t global_relabel_frequency = 5;
		size_t work_since_last_global_relabel = 0, global_relabel_work_threshold = 0;

		size_t out_node_offset = 0, bridge_node_offset = 0;

		Flow upperFlowBound = std::numeric_limits<Flow>::max();

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

		bool winEdge(Node v, Node u) {
			return level[u] == level[v] + 1 || level[u] < level[v] - 1 || (level[u] == level[v] && u < v);
		}

		void clearDatastructures() {
			out_node_offset = hg.numPins();
			bridge_node_offset = 2 * hg.numPins();

			max_level = hg.numNodes() + 2 * hg.numHyperedges();

			flow.assign(2 * hg.numPins() + hg.numHyperedges(), 0);
			excess.assign(max_level, 0);
			level.assign(max_level, 0);

			work_since_last_global_relabel = std::numeric_limits<size_t>::max();
			global_relabel_work_threshold = (global_relabel_alpha * max_level + 2 * hg.numPins() + hg.numHyperedges()) / global_relabel_frequency;
		}
	};
}
