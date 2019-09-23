#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"

namespace whfc {
	class ReachableNodesBase {
	public:
		explicit ReachableNodesBase(const FlowHypergraph& _hg) : hg(_hg) { }
		inline void reach(const Node u) {
			sourceReachableWeight += hg.nodeWeight(u); sourceReachableSize++;
		}

		inline void unreachSource(const Node u) {
			sourceReachableWeight -= hg.nodeWeight(u); sourceReachableSize--;
		}

		inline void unreachTarget(const Node u) {
			targetReachableWeight -= hg.nodeWeight(u); targetReachableSize--;
		}

		inline void settle(const Node u) {
			sourceWeight += hg.nodeWeight(u); sourceSize++;
		}

		void resetSourceReachableToSource() {
			sourceReachableWeight = sourceWeight;
			sourceReachableSize = sourceSize;
		}

		void flipViewDirection() {
			std::swap(sourceReachableWeight, targetReachableWeight);
			std::swap(sourceWeight, targetWeight);
			std::swap(sourceReachableSize, targetReachableSize);
			std::swap(sourceSize, targetSize);
		}

		NodeWeight sourceReachableWeight = NodeWeight(0), sourceWeight = NodeWeight(0), targetReachableWeight = NodeWeight(0), targetWeight = NodeWeight(0);
		NodeIndex sourceReachableSize = NodeIndex(0), sourceSize = NodeIndex(0), targetReachableSize = NodeIndex(0), targetSize = NodeIndex(0);

	protected:
		const FlowHypergraph& hg;

	};


	//ReachableHyperedgesInterface
	//bool allPinsReachable(const Hyperedge e)
	//bool flowSendingPinsReachable(const Hyperedge e)
	//bool reachAllPins(const Hyperedge e)
	//bool reachFlowSendingPins(const Hyperedge e)
}