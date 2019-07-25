#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"

namespace whfc {
	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		NodeWeight maxBlockWeight;
		std::vector<HopDistance> distanceFromCut;

		bool distancesFromCutAvailable() const { return distanceFromCut.size() == cs.flow_hg.numNodes(); }

		bool isBalanced() const {
			const NodeWeight sw = cs.n.sourceReachableWeight, tw = cs.n.targetReachableWeight, total = cs.flow_hg.totalNodeWeight();
			return (sw <= maxBlockWeight && total - sw <= maxBlockWeight) || (tw <= maxBlockWeight && total - tw <= maxBlockWeight);
		}

		int moreBalancedSide() const {
			const NodeWeight sw = cs.n.sourceReachableWeight, tw = cs.n.targetReachableWeight, total = cs.flow_hg.totalNodeWeight();
			const NodeWeight s_disc = std::max(maxBlockWeight - (total - sw), maxBlockWeight - sw);
			const NodeWeight t_disc = std::max(maxBlockWeight - (total - tw), maxBlockWeight - tw);
			return s_disc <= t_disc ? cs.currentViewDirection() : cs.oppositeViewDirection();
		}

		void initialize(Node s, Node t) {
			cs.sourcePiercingNodes = {s};
			cs.targetPiercingNodes = {t};
			grow();
		}

		void pierce() {
			//TODO avoid augmenting paths + distance from cut if available + random
		}

		void advance() {
			pierce();
			grow();
		}

		void grow() {
			if (cs.augmentingPathAvailable) {
				cs.cutSize += flow_algo.exhaustFlow(cs);
				cs.flipViewDirection();
				flow_algo.growReachable(cs);
			}

			if (cs.n.targetReachableWeight > cs.n.sourceReachableWeight)
				cs.flipViewDirection();
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getQueue());
		}

	};
}