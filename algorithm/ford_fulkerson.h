#pragma once

#include "cutter_state.h"
#include "../datastructure/stack.h"
#include "../datastructure/timestamp_reachable_sets.h"

namespace whfc {

	//template<class ScanList>
	class FordFulkerson /* : public FlowAlgorithm */ {
	public:
		using Type = FordFulkerson;
		using ScanList = FixedCapacityStack<Node>;

		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = TimestampReachableNodes;
		//using ReachableHyperedges = TimestampReachableHyperedges;


		using FlowHypergraph::Pin;
		using FlowHypergraph::InHe;

		FlowHypergraph& hg;

		ScanList nodes_to_scan;

		std::vector<InHe> parent;
		Flow capacity = 1 << 24;


		Flow exhaustFlow(CutterState<Type>& cs) {
			return -1;
		}

		template<bool augment_flow>
		Node growWithoutScaling(CutterState<Type>& cs) {
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (Node s : cs.sourcePiercingNodes) nodes_to_scan.push(s);

			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
				for (const InHe& inc_u : hg.hyperedgesOf(u)) {
					const Hyperedge e = inc_u.e;

					if (!h.areAllPinsSourceReachable(e)) {
						const bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(inc_u) > 0;
						if (scanAllPins)
							h.reachAllPins(e);
						else if (h.areFlowSendingPinsSourceReachable(e))
							continue;
						else
							h.reachFlowSendingPins(e);

						for (const Pin& pv : scanAllPins ? hg.pinsOf(e) : hg.pinsSendingFlowInto(e)) {
							const Node v = pv.pin;
							if constexpr (augment_flow)
								if (n.isTarget(v))
									return v;
							AssertMsg(augment_flow || !n.isTargetReachable(v), "Scanned node " + std::to_string(v) + " is reachable from target-side when growing reachable");
							if (!n.isSourceReachable(v)) {		//don't do VD label propagation
								n.reach(v, hg.nodeWeight(u));
								parent[v] = hg.getInHe(pv.he_inc_iter);
								nodes_to_scan.push(v);
							}
						}
					}
				}
			}
			return invalidNode;
		}

		Flow exhaustFlowWithoutScaling(CutterState<Type>& cs) {
			Flow flow = 0;
		}

		Flow growFlowsScaling() {
			return -1;
		}

	};

	//using PseudoDepthFirstFordFulkerson = FordFulkerson<FixedCapacityStack>;
	//using EdmondsKarp = FordFulkerson<LayeredQueue>;
}