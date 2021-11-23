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

		static constexpr bool log = false;

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
			return piercer.findPiercingNode() && (!cs.rejectPiercingIfAugmenting() || !cs.augmentingPathAvailableFromPiercing);
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
					cs.flow_algo.deriveSourceSideCut(false);  // no flow changed --> no new excesses created
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
		std::tuple<bool, bool, size_t>
		enumerateCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			bool time_limit_exceeded = false;
			size_t measure_step = 0;
			auto start_time = std::chrono::high_resolution_clock::now();
			int time_limit = 3600; // seconds
			size_t num_cuts = 0;
			Flow last_cut = 0;

			cs.initialize(s,t);
			bool has_balanced_cut_below_flow_bound = false;
			while (!has_balanced_cut_below_flow_bound && !time_limit_exceeded && findNextCut()) {
				has_balanced_cut_below_flow_bound |= cs.isBalanced();
				if (cs.flow_algo.flow_value != last_cut) {
					last_cut = cs.flow_algo.flow_value;
					num_cuts++;
				}
				if (++measure_step == 500) {
					measure_step = 0;
					auto now = std::chrono::high_resolution_clock::now();
					time_limit_exceeded = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count() > time_limit;
				}
			}

			if (has_balanced_cut_below_flow_bound) {
				if (find_most_balanced && !cs.addingAllUnreachableNodesDoesNotChangeHeavierBlock()) {
					mostBalancedCut();
				}
				else {
					cs.writePartition();
				}
				LOGGER << cs.toString();
			}

			return std::make_tuple(has_balanced_cut_below_flow_bound, time_limit_exceeded, num_cuts);
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
						cs.flow_algo.deriveSourceSideCut(false);
						cs.computeSourceReachableWeight();
						cs.assimilateSourceSide();
					} else {
						cs.flow_algo.deriveTargetSideCut();
						cs.computeTargetReachableWeight();
						cs.assimilateTargetSide();
					}
					cs.side_to_pierce = cs.sideToGrow();
					cs.hasCut = true;
					LOGGER << cs.toString() << V(cs.side_to_pierce);
					cs.verifyCutPostConditions();

					SimulatedNodeAssignment sim = cs.mostBalancedAssignment();
					if (sim.balance() > sol.balance()) {
						sol = sim;
					}
				}

				if (sol.balance() > best_sol.balance()) {
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
