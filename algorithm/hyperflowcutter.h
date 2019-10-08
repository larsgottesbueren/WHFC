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

		static constexpr bool log = true;
		
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight) : hg(hg), cs(hg, maxBlockWeight), flow_algo(hg), upperFlowBound(maxFlow), piercer(hg) { }

		void reset() {
			cs.reset();
			flow_algo.reset();
			upperFlowBound = maxFlow;
			piercer.clear();
		}

		Node selectPiercingNode() {
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;
			
			Assert(cs.n.sourceWeight == cs.n.sourceReachableWeight);
			Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
			cs.cleanUpCut();
			cs.cleanUpBorder();
			cs.verifyCutPostConditions();
			
			Node piercingNode = piercer.findPiercingNode(cs.n, cs.borderNodes, cs.maxBlockWeight);
			if (piercingNode == invalidNode) {
				for (const Node u : hg.nodeIDs())	//Optimization options, if this occurs too often: Track unclaimed nodes in a list that gets filtered on demand.
					if (cs.canBeSettled(u) && hg.nodeWeight(u) + cs.n.sourceWeight <= cs.maxBlockWeight) {
						piercingNode = u;
						break;
					}
			}
			return piercingNode;
		}

		void setPiercingNode(const Node piercingNode) {
			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.emplace_back(piercingNode, cs.n.isTargetReachable(piercingNode));
			cs.settleNode(piercingNode);
			cs.hasCut = false;
		}
		
		bool pierce() {
			Node piercingNode = selectPiercingNode();
			if (piercingNode == invalidNode)
				return false;
			if (cs.flowValue == upperFlowBound && cs.n.isTargetReachable(piercingNode))
				return false;
			setPiercingNode(piercingNode);
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
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			
			cs.hasCut = true;
			if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight) {
				cs.flipViewDirection();
			}
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			LOGGER << cs.toString();
		}

		//for flow-based interleaving
		bool advanceOneFlowIteration() {
			const bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration)
				if (!pierce())
					return false;

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
			
			cs.verifyFlowConstraints();
			
			if (cs.hasCut) {
				cs.verifySetInvariants();
				if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight)
					cs.flipViewDirection();
				GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
				cs.verifyFlowConstraints();
				cs.verifySetInvariants();
				LOGGER << cs.toString();
			}
			
			return true;
		}

		bool runUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			//TODO add most balanced minimum cut
			
			cs.initialize(s,t);
			bool piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = false;
			bool is_balanced = false;
			
																												//no cut ==> run and don't check for balance.
			while (cs.flowValue <= upperFlowBound && (!cs.hasCut || !is_balanced)) {
				piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = !advanceOneFlowIteration();
				if (piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode)
					break;
				is_balanced = cs.isBalanced();
			}
			
			Assert(!cs.hasCut || cs.isBalanced() || cs.flowValue > upperFlowBound);
			
			LOGGER << V(cs.maxBlockWeight) << V(cs.flowValue);
			if (cs.hasCut)
				LOGGER << "has cut" << V(cs.isBalanced());
			else
				LOGGER << "no cut available ==> balance irrelevant";
			LOGGER << V(cs.n.sourceReachableWeight) << V(cs.n.targetReachableWeight) << V(cs.isolatedNodes.weight) << V(cs.unclaimedNodeWeight()) << V(hg.totalNodeWeight());
		
			return !piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode && cs.flowValue <= upperFlowBound && cs.hasCut && is_balanced;
		}
		
		void findCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			exhaustFlowAndGrow();
			while (cs.flowValue <= upperFlowBound && !cs.isBalanced()) {
				if (pierce())
					exhaustFlowAndGrow();
				else
					break;
			}
			LOGGER << V(cs.isBalanced()) << V(cs.maxBlockWeight) << V(cs.flowValue);
			LOGGER << V(cs.n.sourceReachableWeight) << V(cs.n.targetReachableWeight) << V(cs.isolatedNodes.weight) << V(cs.unclaimedNodeWeight()) << V(hg.totalNodeWeight());
			
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here, if desired.
	};

}
