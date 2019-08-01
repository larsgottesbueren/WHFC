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

		Piercer piercer;

		void initialize(Node s, Node t) {
			cs.sourcePiercingNodes = {s};
			cs.targetPiercingNodes = {t};
			exhaustFlowAndGrow();
		}

		void pierce() {
			cs.filterCut();
			cs.filterBorder();

			//For now, we only consider single piercing nodes
			const Node piercingNode = piercer.findPiercingNode(cs.n, cs.borderNodes);
			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.push_back(piercingNode);
			cs.settleNode(piercingNode);
			cs.hasCut = false;
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
		void advanceUntilCut() {
			//TODO verify stopping conditions in advance() functions.
			//This function could be implemented via advanceOneFlowIteration() But the interface to the flow algorithm is so much nicer that I want to keep both
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

			bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration) {
				pierce();
			}

			if (cs.augmentingPathAvailableFromPiercing) {
				if (pierceInThisIteration)
					cs.flowValue += flow_algo.takeFreebie(cs);	//the flow due to freebies does not matter when deciding whether we have a cut available
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
			while (!cs.isBalanced())
				advanceUntilCut();
		}

		void runUntilBalancedAndReportMostBalanced() {
			runUntilBalanced();
			//TODO implement this Pareto point as a sequence of events describing at which step a node joined what side. further which side is the smaller in which step.


			//TODO really figure out the exact details, including IsolatedNodes. This is really overoptimizing stuff in places that won't really matter.
			//But I couldn't bear doing them in a stupid way.
		}

	};

}