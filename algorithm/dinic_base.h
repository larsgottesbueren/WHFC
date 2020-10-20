#pragma once

#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"

namespace whfc {
	class DinicBase {
	public:
		using ScanList = LayeredQueue<Node>;
		
		using ReachableNodes = DistanceReachableNodes;
		using ReachableHyperedges = DistanceReachableHyperedges;
		
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		
		FlowHypergraph& hg;
		LayeredQueue<Node> queue;
		struct StackFrame {
			Node u;
			InHeIndex parent_he_it;
		};
		FixedCapacityStack<StackFrame> stack;
		int direction = 0;
		std::vector<PinIndex> current_flow_sending_pin, current_flow_receiving_pin, current_pin;
		std::vector<InHeIndex> current_hyperedge;
		
		Flow upperFlowBound = std::numeric_limits<Flow>::max();
		
		DinicBase(FlowHypergraph& hg) : hg(hg), queue(hg.numNodes()), stack(hg.numNodes()), current_flow_sending_pin(hg.numHyperedges(), PinIndex::Invalid()), current_flow_receiving_pin(hg.numHyperedges(), PinIndex::Invalid()), current_pin(hg.numHyperedges(), PinIndex::Invalid()), current_hyperedge(hg.numNodes(), InHeIndex::Invalid()) {
		
		}
		
		void flipViewDirection() {
			std::swap(current_flow_sending_pin, current_flow_receiving_pin);
			direction = 1 - direction;
		}
		
	};
}