
#include "cutter_state.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"

namespace whfc {
	class Dinic {
	public:
		using Type = Dinic;
		using ScanList = FixedCapacityStack<Node>;

		using ReachableNodes = DistanceReachableNodes;
		using ReachableHyperedges = DistanceReachableHyperedges;

	};
}