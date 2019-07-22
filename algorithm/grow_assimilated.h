#pragma once

#include "../datastructure/queue.h"
#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "../datastructure/stack.h"


namespace whfc {

template<typename FlowAlgorithm>
class GrowAssimilated {
public:
	//using ScanList = FlowAlgorithm::ScanList;
	using ScanList = FixedCapacityStack<Node>;
	using Pin = FlowHypergraph::Pin;

	//using ReachableNodes = typename FlowAlgorithm::ReachableNodes;
	//using ReachableHyperedges = typename FlowAlgorithm::ReachableHyperedges;

	using ReachableNodes = BitsetReachableNodes;
	using ReachableHyperedges = BitsetReachableHyperedges;


	static void grow(CutterState<FlowAlgorithm>& cs, ScanList& nodes_to_scan) {
		ReachableNodes& n = cs.n;
		ReachableHyperedges& h = cs.h;
		nodes_to_scan.clear();
		FlowHypergraph& hg = cs.flow_hg;

		for (Node s : cs.sourcePiercingNodes) {
			nodes_to_scan.push(s);
			AssertMsg(n.isSource(s), "source-side piercing node " + std::to_string(s) + " not a source");
			AssertMsg(!n.isTarget(s), "source-side piercing node " + std::to_string(s) + "  is target");
			AssertMsg(!n.isTargetReachable(s), "source-side piercing node " + std::to_string(s) + " is reachable from target-side");
		}

		while (!nodes_to_scan.empty()) {
			const Node u = nodes_to_scan.pop();
			for (const auto& he_inc : hg.hyperedgesOf(u)) {
				const Hyperedge e = he_inc.e;

				if (!h.areAllPinsSources(e)) {
					const bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(he_inc) > 0;
					if (scanAllPins)
						h.settleAllPins(e);
					else if (h.areFlowSendingPinsSources(e))
						continue;
					else
						h.settleFlowSendingPins(e);

					for (const Pin& pv : scanAllPins ? hg.pinsOf(e) : hg.pinsSendingFlowInto(e)) {
						const Node v = pv.pin;
						AssertMsg(!n.isTargetReachable(v), "Scanned node " + std::to_string(v) + " is reachable from target-side");
						if (!n.isSource(v)) {
							cs.settleNode(v);
							nodes_to_scan.push(v);
						}
					}
				}
			}
		}
	}


};

}