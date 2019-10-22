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
		Piercer<FlowAlgorithm> piercer;
		bool find_most_balanced = true;

		static constexpr bool log = true;
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight, int seed) :
				timer(algo_name),
				hg(hg),
				cs(hg, maxBlockWeight, timer),
				flow_algo(hg),
				upperFlowBound(maxFlow),
				piercer(hg, cs)
		{
			Random::setSeed(seed);
		}

		void reset() {
			cs.reset();
			flow_algo.reset();
			upperFlowBound = maxFlow;
			piercer.clear();
			timer.clear();
		}
		
		void flipViewDirection() {
			cs.flipViewDirection();
			piercer.flipViewDirection();
		}

		void setPiercingNode(const Node piercingNode) {
			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.emplace_back(piercingNode, cs.n.isTargetReachable(piercingNode));
			cs.settleNode(piercingNode);
			cs.hasCut = false;
		}
		
		bool pierce(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			timer.start("Piercing");
			Node piercingNode = piercer.findPiercingNode();
			timer.stop("Piercing");
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
				flipViewDirection();
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
				flipViewDirection();
			}
			timer.start("Grow Assimilated");
			LOGGER << "grow assimilated";
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			timer.stop("Grow Assimilated");
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			LOGGER << cs.toString();
		}

		//for flow-based interleaving
		bool advanceOneFlowIteration(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			const bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration)
				if (!pierce(reject_piercing_if_it_creates_an_augmenting_path))
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
					flipViewDirection();
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
					flipViewDirection();
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
			cs.initialize(s,t);
			bool piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = false;
			bool has_balanced_cut = false;
			
			while (cs.flowValue <= upperFlowBound && !has_balanced_cut) {
				piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = !advanceOneFlowIteration(cs.flowValue == upperFlowBound);
				if (piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode)
					break;
				has_balanced_cut = cs.hasCut && cs.isBalanced(); //no cut ==> run and don't check for balance.
			}
			
			if (has_balanced_cut && cs.flowValue <= upperFlowBound) {
				// S + U + ISO <= T ==> will always add U and ISO completely to S
				const bool better_balance_impossible = hg.totalNodeWeight() - cs.n.targetReachableWeight <= cs.n.targetReachableWeight;
				LOGGER << V(better_balance_impossible);
				if (find_most_balanced && !better_balance_impossible)
					mostBalancedMinimumCut();
				else
					cs.mostBalancedIsolatedNodesAssignment();
			}
			
			Assert(!cs.hasCut || cs.isBalanced() || cs.flowValue > upperFlowBound);
			return !piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode && cs.flowValue <= upperFlowBound && has_balanced_cut;
		}
		
		void mostBalancedMinimumCut() {
			timer.start("MBMC");
			
			//settle target reachable nodes, so we don't have to track them in the moves
			flipViewDirection();
			timer.start("Grow Assimilated");
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			timer.stop("Grow Assimilated");
			flipViewDirection();
			
			NodeWeight best_bwd = cs.mostBalancedIsolatedNodesAssignment(false);
			cs.mostBalancedMinimumCutMode = true;	//enables move tracker
			
			// Option 1: implement MBMC as KaHyPar does. always grow only one side. they do target side because of reverse topo sweep. we could do either.
			// Option 2:
			
			// Options:
			// 1) save the state every time we would flip the search direction.
			// starting with the first time flipping, this would probably just alternate. so that's expensive
			//
			// 2) save before the first flip and again at the end.
			// 3) make fully dynamic, i.e. store only deltas and then roll back. this is quite nasty
			
			// This implements option 2) from above. not so pretty.
			
			// Option 2) is still too slow (slower than computing flow) --> so think about how to make it fully dynamic or find cases where we can output less often
			
			// instead of the partition saved in ReachableNodes : provide a set of node moves.
			
			auto saved_state = cs.saveState();
			std::vector<typename CutterState<FlowAlgorithm>::SolutionState> solutions;
			size_t restore_to = 0;
			NodeWeight bwd_at_first_partition = saveState_OutputPartition_Restore(solutions);
			
			if (bwd_at_first_partition > 0) {
				bool changed = false;
				int dir = cs.currentViewDirection();
				while (cs.currentViewDirection() == dir && advanceOneFlowIteration(true /* reject piercing if it creates an augmenting path */))
					changed = true;
				
				if (changed && cs.isBalanced()) {
					NodeWeight bwd_at_first_flip = saveState_OutputPartition_Restore(solutions);
					if (bwd_at_first_flip < bwd_at_first_partition)
						restore_to = 1;
					
					changed = false;
					while (advanceOneFlowIteration(true))
						changed = true;
					
					if (changed && cs.isBalanced()) {
						NodeWeight end_bwd = saveState_OutputPartition_Restore(solutions);    //one state save too many. probably not problematic
						if (end_bwd < std::min(bwd_at_first_partition, bwd_at_first_flip))
							restore_to = 2;
					}
				}
			}
			cs.restoreState(solutions[restore_to]);
			timer.stop("MBMC");
		}
		
		NodeWeight saveState_OutputPartition_Restore(std::vector<typename CutterState<FlowAlgorithm>::SolutionState>& solutions) {
			auto saved = cs.saveState();
			NodeWeight bwd = cs.mostBalancedIsolatedNodesAssignment();
			solutions.emplace_back(cs.saveState());
			cs.restoreState(saved);
			cs.partitionWrittenToNodeSet = false;
			return bwd;
		}
		
		void findCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			exhaustFlowAndGrow();
			while (cs.flowValue <= upperFlowBound && !cs.isBalanced() && pierce())
				exhaustFlowAndGrow();
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here, if desired.
	};

}
