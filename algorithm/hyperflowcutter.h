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
		Piercer<FlowAlgorithm> piercer;
		bool find_most_balanced = true;

		static constexpr bool log = false;

		HyperFlowCutter(FlowHypergraph& hg, int seed) :
				timer("HyperFlowCutter"),
				hg(hg),
				cs(hg, timer),
				piercer(hg, cs, timer)
		{
			cs.rng.setSeed(seed);
			reset();
		}

		void reset() {
			cs.reset();
			piercer.reset();
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


		bool findNextCut(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			if (cs.hasCut) {	// false on the first call, true on subsequent calls.
				const bool early_reject = !pierce(reject_piercing_if_it_creates_an_augmenting_path);
				if (early_reject) {
					return false;
				}
			}

			if (cs.augmentingPathAvailableFromPiercing) {
				cs.hasCut = cs.flow_algo.findMinCuts();
			}
			else {
				if (cs.side_to_pierce == 0) {
					cs.flow_algo.deriveSourceSideCut();
				} else {
					cs.flow_algo.deriveTargetSideCut();
				}
				cs.hasCut = true;	// no flow increased
			}

			if (cs.hasCut) {
				cs.assimilate();
			}

			return cs.hasCut && cs.flowValue <= cs.flow_algo.upper_flow_bound;
		}


		/*
		 * Equivalent to runUntilBalancedOrFlowBoundExceeded(s,t) except that it does not use the flow-based interleaving that is necessary when running multiple HFC instances
		 */
		bool enumerateCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			bool has_balanced_cut_below_flow_bound = false;
			while (!has_balanced_cut_below_flow_bound && findNextCut(cs.flowValue == cs.flow_algo.upper_flow_bound)) {
				has_balanced_cut_below_flow_bound |= cs.isBalanced();
			}

			if (has_balanced_cut_below_flow_bound) {
				// TODO this doesnt work without directions
				const double imb_S_U_ISO = static_cast<double>(hg.totalNodeWeight() - cs.n.targetReachableWeight) / static_cast<double>(cs.maxBlockWeight(cs.currentViewDirection()));
				const double imb_T = static_cast<double>(cs.n.targetReachableWeight) / static_cast<double>(cs.maxBlockWeight(cs.oppositeViewDirection()));
				const bool better_balance_impossible = cs.unclaimedNodeWeight() == 0 || imb_S_U_ISO <= imb_T;
				if (find_most_balanced && !better_balance_impossible) {
					mostBalancedCut();
				}
				else {
					cs.writePartition();
				}

				LOGGER << cs.toString(true);
				cs.verifyCutInducedByPartitionMatchesFlowValue();
			}

			return has_balanced_cut_below_flow_bound;
		}

		void mostBalancedCut() {
			timer.start("MBMC");
			LOGGER << "MBC Mode";
			// assimilate the missing side, so we don't have to track it in the moves
			if (cs.side_to_pierce == 0) {
				cs.assimilateTargetSide();
			} else {
				cs.assimilateSourceSide();
			}
			assert(cs.source_reachable_weight == cs.source_weight);
			assert(cs.target_reachable_weight == cs.target_weight);

			NonDynamicCutterState first_balanced_state = cs.enterMostBalancedCutMode();
			SimulatedNodeAssignment initial_sol = cs.mostBalancedAssignment();
			std::vector<Move> best_moves;
			SimulatedNodeAssignment best_sol = initial_sol;

			const size_t mbc_iterations = 7;
			for (size_t i = 0; i < mbc_iterations && !best_sol.isPerfectlyBalanced(); ++i) {
				LOGGER << "MBC it" << i;
				SimulatedNodeAssignment sol = best_sol;
				while (!sol.isPerfectlyBalanced() && pierce(true)) {
					if (cs.side_to_pierce == 0) {
						cs.flow_algo.deriveSourceSideCut();
						cs.computeSourceReachableWeight();
						cs.assimilateSourceSide();
					} else {
						cs.flow_algo.deriveTargetSideCut();
						cs.computeTargetReachableWeight();
						cs.assimilateTargetSide();
					}
					cs.side_to_pierce = cs.sideToGrow();
					cs.hasCut = true;
					cs.verifyCutPostConditions();
					LOGGER << cs.toString();

					SimulatedNodeAssignment sim = cs.mostBalancedAssignment();
					if (sim.imbalance() < sol.imbalance()) {
						sol = sim;
					}
				}

				if (sol.imbalance() < best_sol.imbalance()) {
					best_sol = sol;
					cs.revertMoves(sol.numberOfTrackedMoves);
					best_moves = cs.trackedMoves;
				}
				cs.resetToFirstBalancedState(first_balanced_state);
			}

			cs.applyMoves(best_moves);
			cs.writePartition(best_sol);

			timer.stop("MBMC");
		}

		void signalTermination() {
			cs.flow_algo.shall_terminate = true;
		}

		void setFlowBound(Flow bound) {
			cs.flow_algo.upper_flow_bound = bound;
		}
	};

}
