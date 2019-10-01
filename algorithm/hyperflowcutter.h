#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"


namespace whfc {
	
	//TODO flow assertions

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		Flow upperFlowBound;
		Piercer piercer;

		static constexpr bool debug = true;
		
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight) : hg(hg), cs(hg, maxBlockWeight), flow_algo(hg), upperFlowBound(maxFlow), piercer(hg) { }

		void initialize(Node s, Node t) {
			cs.initialize(s,t);
			exhaustFlowAndGrow();
		}
		
		void clear() {
			cs.clear();
			upperFlowBound = maxFlow;
			piercer.clear();
		}

		bool pierce() {
			cs.cleanUpCut();
			cs.cleanUpBorder();
			Assert(static_cast<uint32_t>(cs.flowValue) == cs.cut.weight(hg));
			
			if (cs.notSettledNodeWeight() == 0)
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
			cs.sourcePiercingNodes.emplace_back(piercingNode, cs.n.isTargetReachable(piercingNode));
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
			
			cs.n.verifyDisjoint();
			cs.n.verifySettledIsSubsetOfReachable();
			cs.h.verifyDisjoint();
			cs.h.verifySettledIsSubsetOfReachable();
			
			
			if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight)
				cs.flipViewDirection();
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			cs.hasCut = true;
			
			
			cs.n.verifyDisjoint();
			cs.n.verifySettledIsSubsetOfReachable();
			cs.h.verifyDisjoint();
			cs.h.verifySettledIsSubsetOfReachable();
			LOG << cs.toString();
		}

		//for cut-based interleaving
		//This function could be implemented via advanceOneFlowIteration() But the interface to the flow algorithm is so much nicer that I want to keep both
		void advanceUntilCut() {
			AssertMsg(cs.hasCut, "Advancing until cut, but hasCut flag not set");
			if (!pierce())
				throw std::runtime_error("Piercing went wrong.");
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
			LOG << V(cs.isBalanced()) << V(cs.maxBlockWeight) << V(cs.flowValue);
			LOG << V(cs.n.sourceReachableWeight) << V(cs.n.targetReachableWeight) << V(cs.isolatedNodes.weight) << V(cs.unclaimedNodeWeight()) << V(hg.totalNodeWeight());
			
		}

		void runUntilBalancedAndReportMostBalanced() {
			runUntilBalanced();
			if (cs.flowValue >= upperFlowBound)
				return;
			//either run until cut would have to be broken. or track most balanced (just as a partition) and revert to that
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here, if desired.
	};

}