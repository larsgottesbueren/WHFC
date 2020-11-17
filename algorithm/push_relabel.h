#include "cutter_state.h"

#include "../datastructure/queue.h"
#include "../datastructure/stack.h"

#include <boost/circular_buffer.hpp>

#include "../datastructure/bitset_reachable_sets.h"

namespace whfc {

	class PushRelabel{
	public:
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;

		using Type = PushRelabel;
		using ScanList = LayeredQueue<Node>;
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		using ReachableNodes = BitsetReachableNodes;			// TODO temporary so CutterState won't complain.
		using ReachableHyperedges = BitsetReachableHyperedges;

		FlowHypergraph& hg;

		LayeredQueue<Node> bfs_queue;		// for global relabeling
		FixedCapacityStack<Pin> stack;	// for flow decomposition
		int direction = 0, previous_cutter_state_direction = 0;

		std::vector<PinIndex> current_pin;
		std::vector<InHeIndex> current_hyperedge;

		std::vector<Flow> excess;
		boost::circular_buffer<Node> active_vertices_and_edges;
		BitVector in_active_queue;
		std::vector<uint32_t> level;
		uint32_t max_level;

		Flow upperFlowBound = std::numeric_limits<Flow>::max();

		PushRelabel(FlowHypergraph& hg) :
				hg(hg), bfs_queue(hg.numNodes()), stack(hg.numNodes()),
				current_pin(hg.numHyperedges(), PinIndex::Invalid()),
				current_hyperedge(hg.numNodes(), InHeIndex::Invalid()),
				excess(hg.numNodes() + hg.numHyperedges(), 0),
				in_active_queue(hg.numNodes() + hg.numHyperedges(), 0),
				level(hg.numNodes() + hg.numHyperedges(), 0)
		{

		}

		void flipViewDirection() {
			direction = 1 - direction;
		}


		void reset() {

		}

		ScanList& getScanList() {
			return bfs_queue;
		}

		bool exhaustFlow(CutterState<Type>& cs) {
			prepare(cs);
			max_level = hg.numNodes() + hg.numHyperedges();
			const Node target = cs.targetPiercingNodes.front().node;
			Flow old_excess = excess[target];

			while ( -excess[target] <= upperFlowBound && !active_vertices_and_edges.empty() ) {
				const Node x = active_vertices_and_edges.front();
				active_vertices_and_edges.pop_front();
				discharge(x);
				if (excess[x] > 0 && level[x] < max_level) {
					active_vertices_and_edges.push_back(x);
				} else {
					in_active_queue.reset(x);
				}
			}

			for (Node x : active_vertices_and_edges) {
				in_active_queue.reset(x);
			}
			active_vertices_and_edges.clear();

			finish(cs);
			timer.report(std::cout);

			Flow flow_delta = old_excess - excess[target];
			assert(flow_delta >= 0);
			cs.flowValue += flow_delta;
			return flow_delta > 0;
		}

		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			prepare(cs);
			Flow f = 0;
			finish(cs);
			return f;
		}

		void growReachable(CutterState<Type>& cs) {
			prepare(cs);
			finish(cs);
		}

		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			return 0;		// can't do anything here
		}

	private:
		static constexpr bool log = false;

		void discharge(Node u) {
			if (excess[u] > 0 && level[u] < max_level) {
				if (u < hg.numNodes()) {
					dischargeNode(u);
				} else {
					dischargeEdge(Hyperedge(u - hg.numNodes()), u);
				}
			}
		}

		void dischargeNode(Node u) {

		}

		void dischargeEdge(Hyperedge e, Node e_node) {

		}

		void prepare(CutterState<Type>& cs) {
			previous_cutter_state_direction = cs.currentViewDirection();
			if (previous_cutter_state_direction != 0) {
				cs.flipViewDirection();
			}
			assert(cs.currentViewDirection() == 0);
		}

		void finish(CutterState<Type>& cs) {

		}

		TimeReporter timer;
	};
}