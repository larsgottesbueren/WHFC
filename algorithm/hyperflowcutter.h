#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
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

		static constexpr bool log = true;

		HyperFlowCutter(FlowHypergraph& hg, int seed) :
				timer("HyperFlowCutter"),
				hg(hg),
				cs(hg, timer),
				piercer(hg, cs)
		{
			cs.rng.setSeed(seed);
			reset();
		}

		void reset() {
			cs.reset();
			piercer.reset();
			//timer.clear();
		}

		bool pierce() {
			Node piercingNode = piercer.findPiercingNode();
			LOGGER << V(piercingNode);
			if (piercingNode == invalidNode)
				return false;
			if (cs.rejectPiercingIfAugmenting() && cs.reachableFromSideNotToPierce(piercingNode))
				return false;
			cs.setPiercingNode(piercingNode);
			return true;
		}


		bool findNextCut() {
			if (cs.hasCut && !pierce()) {
				return false;
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

			return cs.hasCut && cs.flow_algo.flow_value <= cs.flow_algo.upper_flow_bound;
		}


		/*
		 * Equivalent to runUntilBalancedOrFlowBoundExceeded(s,t) except that it does not use the flow-based interleaving that is necessary when running multiple HFC instances
		 */
		bool enumerateCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			bool has_balanced_cut_below_flow_bound = false;
			while (!has_balanced_cut_below_flow_bound && findNextCut()) {
				has_balanced_cut_below_flow_bound |= cs.isBalanced();
			}

			if (has_balanced_cut_below_flow_bound) {
				if (find_most_balanced && !cs.betterBalanceImpossible()) {
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
				while (!sol.isPerfectlyBalanced() && pierce()) {
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
