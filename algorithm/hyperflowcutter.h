#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"

namespace whfc {
	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm

		void pierce() {

		}

		void advance() {
			pierce();
			flow_algo.exhaustFlow(cs);
			flow_algo.growReachable(cs);
			//output cut
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getQueue());
		}

	};
}