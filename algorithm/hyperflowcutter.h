#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"


namespace whfc {

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		NodeWeight maxBlockWeight;
		Piercer piercer;



		bool isBalanced() const {
			const NodeWeight sw = cs.n.sourceReachableWeight, tw = cs.n.targetReachableWeight, total = cs.hg.totalNodeWeight();
			return (sw <= maxBlockWeight && total - sw <= maxBlockWeight) || (tw <= maxBlockWeight && total - tw <= maxBlockWeight);
		}

		int moreBalancedSide() const {
			const NodeWeight sw = cs.n.sourceReachableWeight, tw = cs.n.targetReachableWeight, total = cs.hg.totalNodeWeight();
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
			cs.filterCut();
			cs.filterBorder();

			//For now, we only consider single piercing nodes
			const Node piercingNode = piercer.findPiercingNode(cs.n, cs.borderNodes);
			cs.augmentingPathAvailable = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.push_back(piercingNode);
			cs.settleNode(piercingNode);
		}

		void grow() {
			if (cs.augmentingPathAvailable) {
				cs.cutSize += flow_algo.exhaustFlow(cs);
				cs.flipViewDirection();
				flow_algo.growReachable(cs);
			}
			else {
				flow_algo.growReachable(cs);
			}

			if (cs.n.targetReachableWeight > cs.n.sourceReachableWeight)
				cs.flipViewDirection();
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getQueue());
		}

		void advance() {
			//TODO verify stopping condiitions
			pierce();
			grow();
		}

		void runUntilBalanced() {
			while (!isBalanced())
				advance();
		}

	};

	//TODO

	//TODO multicutter ? potentially parallel
}