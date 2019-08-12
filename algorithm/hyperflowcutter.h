#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"


namespace whfc {

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		Flow upperFlowBound;
		Piercer piercer;

		void initialize(Node s, Node t) {
			cs.sourcePiercingNodes = {s};
			cs.targetPiercingNodes = {t};
			exhaustFlowAndGrow();
		}

		bool pierce() {
			cs.filterCut();
			cs.filterBorder();
			if (cs.unclaimedNodeWeight() == 0)
				return false;

			Node piercingNode = piercer.findPiercingNode(cs.n, cs.borderNodes, cs.maxBlockWeight);

			if (piercingNode == invalidNode) {
				for (const Node u : hg.nodeIDs())	//Optimization options, if this occurs too often: Track unclaimed nodes in a list that gets filtered on demand.
					if (cs.canBeSettled(u) && hg.nodeWeight(u) + cs.n.sourceWeight <= cs.maxBlockWeight) {
						piercingNode = u;
						break;
					}
				if (piercingNode == invalidNode)
					return false;
			}

			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.push_back(piercingNode);
			cs.settleNode(piercingNode);
			cs.hasCut = false;
			return true;
		}

		void exhaustFlowAndGrow() {
			if (cs.augmentingPathAvailableFromPiercing) {
				cs.flowValue += flow_algo.exhaustFlow(cs);
				cs.flipViewDirection();
				flow_algo.growReachable(cs);
			}
			else {
				flow_algo.growReachable(cs);
			}

			if (cs.n.targetReachableWeight < cs.n.sourceReachableWeight)
				cs.flipViewDirection();
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getQueue());
			cs.hasCut = true;
		}

		//for cut-based interleaving
		//This function could be implemented via advanceOneFlowIteration() But the interface to the flow algorithm is so much nicer that I want to keep both
		void advanceUntilCut() {
			AssertMsg(cs.hasCut, "Advancing until cut, but hasCut flag not set");
			pierce();
			exhaustFlowAndGrow();
		}

		//for flow-based interleaving
		void advanceOneFlowIteration() {
			//currently we return after every cut found (potentially multiple with the same cutsize due to aap).
			//when aap applies, this does not incur linear work, i.e. when parallelizing/interleaving, all aap rounds should be bundled with its previous flow-increasing round
			//we want to report all of these aap cuts. therefore the surrounding caller has to perform the bundling
			//could do an example implementation with a report callback here.

			const bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration) {
				pierce();
			}

			if (cs.augmentingPathAvailableFromPiercing) {
				if (pierceInThisIteration)
					cs.flowValue += flow_algo.recycleDatastructuresFromGrowReachablePhase(cs);	//the flow due to recycled datastructures does not matter when deciding whether we have a cut available
				Flow flow_diff = flow_algo.growFlowOrSourceReachable(cs);
				cs.flowValue += flow_diff;
				cs.hasCut = flow_diff == 0;
				if (cs.hasCut) {
					cs.flipViewDirection();
					flow_algo.growReachable(cs);
				}
			}
			else {
				flow_algo.growReachable(cs);		//don't grow target reachable
				cs.hasCut = true;
			}

			if (cs.hasCut) {
				if (cs.n.targetReachableWeight < cs.n.sourceReachableWeight)
					cs.flipViewDirection();
				GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getQueue());
			}
		}

		void runUntilBalanced() {
			while (!cs.isBalanced() && cs.flowValue < upperFlowBound)
				advanceUntilCut();
		}

		void runUntilBalancedAndReportMostBalanced() {
			runUntilBalanced();
			if (cs.flowValue >= upperFlowBound)
				return;
			//either run until cut would have to be broken. or track most balanced and revert to that
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here.
	};

}