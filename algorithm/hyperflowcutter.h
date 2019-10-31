#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"

namespace whfc {

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		TimeReporter timer;
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		Flow upperFlowBound;
		Piercer<FlowAlgorithm> piercer;
		bool find_most_balanced = true;

		static constexpr bool log = true;
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight, int seed) :
				timer("HyperFlowCutter"),
				hg(hg),
				cs(hg, maxBlockWeight, timer),
				flow_algo(hg),
				upperFlowBound(maxFlow),
				piercer(hg, cs, timer)
		{
			Random::setSeed(seed);
		}

		void reset() {
			cs.reset();
			flow_algo.reset();
			upperFlowBound = maxFlow;
			//timer.clear();
		}
		
		void setPiercingNode(const Node piercingNode) {
			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.emplace_back(piercingNode, cs.n.isTargetReachable(piercingNode));
			cs.settleNode(piercingNode);
			cs.hasCut = false;
		}
		
		bool pierce(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			Node piercingNode = piercer.findPiercingNode();
			if (piercingNode == invalidNode)
				return false;
			if (reject_piercing_if_it_creates_an_augmenting_path && cs.n.isTargetReachable(piercingNode))
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
			cs.verifyCutPostConditions();
			LOGGER << cs.toString();
		}

		//for flow-based interleaving
		bool advanceOneFlowIteration(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			const bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration) {
				timer.start("Piercing");
				const bool early_reject = !pierce(reject_piercing_if_it_creates_an_augmenting_path);
				timer.stop("Piercing");
				if (early_reject)
					return false;
			}
			
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
				cs.verifyCutPostConditions();
				LOGGER << cs.toString();
			}
			
			return true;
		}

		bool runUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			bool piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = false;
			bool has_balanced_cut = false;
			
			while (cs.flowValue <= upperFlowBound && !has_balanced_cut) {
				piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = !advanceOneFlowIteration(cs.flowValue == upperFlowBound);
				if (piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode)
					break;
				has_balanced_cut = cs.hasCut && cs.isBalanced(); //no cut ==> run and don't check for balance.
			}
			
			LOGGER << V(has_balanced_cut) << V(cs.flowValue) << V(upperFlowBound);
			
			if (has_balanced_cut && cs.flowValue <= upperFlowBound) {
				// S + U + ISO <= T ==> will always add U and ISO completely to S, i.e. take target-side cut (we know S <= T)
				const bool better_balance_impossible = hg.totalNodeWeight() - cs.n.targetReachableWeight <= cs.n.targetReachableWeight || cs.unclaimedNodeWeight() == 0;
				LOGGER << V(better_balance_impossible);
				if (find_most_balanced && !better_balance_impossible)
					mostBalancedCut();
				else
					cs.writePartition();
				
				LOGGER << cs.toString(true);
				cs.verifyCutInducedByPartitionMatchesFlowValue();
			}
			return !piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode && cs.flowValue <= upperFlowBound && has_balanced_cut;
		}
		
		void mostBalancedCut() {
			timer.start("MBMC");
			
			LOGGER << "MBC Mode";
			
			//settle target reachable nodes, so we don't have to track them in the moves
			cs.flipViewDirection();
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			cs.verifyCutPostConditions();
			cs.flipViewDirection();
			
			Assert(cs.n.sourceReachableWeight == cs.n.sourceWeight);
			Assert(cs.n.targetReachableWeight == cs.n.targetWeight);
			
			NonDynamicCutterState first_balanced_state = cs.enterMostBalancedCutMode();
			SimulatedIsolatedNodesAssignment initial_sol = cs.mostBalancedIsolatedNodesAssignment();
			std::vector<Move> best_moves;
			SimulatedIsolatedNodesAssignment best_sol = initial_sol;
			
			const size_t mbc_iterations = 7;
			for (size_t i = 0; i < mbc_iterations && best_sol.blockWeightDiff > 0; ++i) {
				LOGGER << "MBC it" << i;
				Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
				SimulatedIsolatedNodesAssignment sol = best_sol;
				while (sol.blockWeightDiff > 0 && pierce(true)) {
					flow_algo.growReachable(cs);		//could consolidate for factor 2 speedup. but avoid any code duplication.
					GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
					cs.hasCut = true;
					cs.verifyCutPostConditions();
					
					LOGGER << cs.toString();
					
					if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight)
						cs.flipViewDirection();
					
					SimulatedIsolatedNodesAssignment sim = cs.mostBalancedIsolatedNodesAssignment();
					if (sim.blockWeightDiff < sol.blockWeightDiff)
						sol = sim;
				}
				
				if (sol.blockWeightDiff < best_sol.blockWeightDiff) {
					best_sol = sol;
					cs.revertMoves(sol.numberOfTrackedMoves);
					LOGGER << "improved";
					best_moves = cs.trackedMoves;
				}
				
				cs.resetToFirstBalancedState(first_balanced_state);
			}
			
			cs.applyMoves(best_moves);
			cs.writePartition(best_sol);
			
			timer.stop("MBMC");
		}

		
		void findCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			exhaustFlowAndGrow();
			while (cs.flowValue <= upperFlowBound && !cs.isBalanced() && pierce())
				exhaustFlowAndGrow();
		}
		
	};

}
