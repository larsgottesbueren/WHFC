#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"



namespace whfc {

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		static constexpr const char *algo_name = "HyperFlowCutter";		//this looks weird
		TimeReporter timer;
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		Flow upperFlowBound;
		Piercer piercer;

		static constexpr bool log = true;
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight) : timer(algo_name), hg(hg), cs(hg, maxBlockWeight, timer), flow_algo(hg), upperFlowBound(maxFlow), piercer(hg) { }

		void reset() {
			cs.reset();
			flow_algo.reset();
			upperFlowBound = maxFlow;
			piercer.clear();
			timer.clear();
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
			timer.start("Piercing");
			Node piercingNode = selectPiercingNode();
			timer.stop("Piercing");
			if (piercingNode == invalidNode)
				return false;
			if (cs.flowValue == upperFlowBound && cs.n.isTargetReachable(piercingNode))
				return false;
			setPiercingNode(piercingNode);
			return true;
		}
		
		void exhaustFlowAndGrow() {
			timer.start("Flow");
			if (cs.augmentingPathAvailableFromPiercing) {
				timer.start("Augment", "Flow");
				cs.flowValue += flow_algo.exhaustFlow(cs);
				timer.stop("Augment");
				cs.flipViewDirection();
				timer.start("Grow Backward Reachable", "Flow");
				flow_algo.growReachable(cs);
				timer.stop("Grow Backward Reachable");
			}
			else {
				timer.start("Grow Reachable due to AAP", "Flow");
				flow_algo.growReachable(cs);
				timer.stop("Grow Reachable due to AAP");
			}
			timer.stop("Flow");
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			
			cs.hasCut = true;
			if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight) {
				cs.flipViewDirection();
			}
			timer.start("Grow Assimilated");
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			timer.stop("Grow Assimilated");
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
			
			timer.start("Flow");
			if (cs.augmentingPathAvailableFromPiercing) {
				timer.start("Augment", "Flow");
				if (pierceInThisIteration)
					cs.flowValue += flow_algo.recycleDatastructuresFromGrowReachablePhase(cs);	//the flow due to recycled datastructures does not matter when deciding whether we have a cut available
				Flow flow_diff = flow_algo.growFlowOrSourceReachable(cs);
				timer.stop("Augment");
				cs.flowValue += flow_diff;
				cs.hasCut = flow_diff == 0;
				if (cs.hasCut) {
					cs.flipViewDirection();
					timer.start("Grow Backward Reachable", "Flow");
					flow_algo.growReachable(cs);
					timer.stop("Grow Backward Reachable");
				}
			}
			else {
				timer.start("Grow Reachable due to AAP", "Flow");
				flow_algo.growReachable(cs);		//don't grow target reachable
				timer.stop("Grow Reachable due to AAP");
				cs.hasCut = true;
			}
			timer.stop("Flow");
			
			cs.verifyFlowConstraints();
			
			if (cs.hasCut) {
				cs.verifySetInvariants();
				if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight)
					cs.flipViewDirection();
				timer.start("Grow Assimilated");
				GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
				timer.stop("Grow Assimilated");
				cs.verifyFlowConstraints();
				cs.verifySetInvariants();
				LOGGER << cs.toString();
			}
			
			return true;
		}

		bool runUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			//TODO add most balanced minimum cut, once is_balanced becomes true
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
			while (cs.flowValue <= upperFlowBound && !cs.isBalanced() && pierce())
				exhaustFlowAndGrow();
			LOGGER << V(cs.isBalanced()) << V(cs.maxBlockWeight) << V(cs.flowValue);
			LOGGER << V(cs.n.sourceReachableWeight) << V(cs.n.targetReachableWeight) << V(cs.isolatedNodes.weight) << V(cs.unclaimedNodeWeight()) << V(hg.totalNodeWeight());
			
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here, if desired.
	};

}
