#pragma once

#include "../datastructure/queue.h"
#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "../datastructure/stack.h"


namespace whfc {

template<typename FlowAlgorithm>
class GrowAssimilated {
public:
	using Pin = FlowHypergraph::Pin;
	
	using ScanList = typename FlowAlgorithm::ScanList;
	using ReachableNodes = typename FlowAlgorithm::ReachableNodes;
	using ReachableHyperedges = typename FlowAlgorithm::ReachableHyperedges;

	
	static void grow(CutterState<FlowAlgorithm>& cs, ScanList& nodes_to_scan, const bool reach_and_settle = false) {
		unused(reach_and_settle);
		
		ReachableNodes& n = cs.n;
		ReachableHyperedges& h = cs.h;
		nodes_to_scan.clear();
		FlowHypergraph& hg = cs.hg;

		for (auto& ps : cs.sourcePiercingNodes) {
			const Node s = ps.node;
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
					if (scanAllPins) {
						Assert(h.areAllPinsSourceReachable(e) || reach_and_settle);
						cs.settleAllPins(e);
						
						if (FlowAlgorithm::grow_reachable_marks_flow_sending_pins_when_marking_all_pins) {
							Assert(h.areFlowSendingPinsSourceReachable(e) || reach_and_settle);
							if (!h.areFlowSendingPinsSources(e))
								cs.settleFlowSendingPins(e);
						}
					}
					else {
						if (cs.shouldBeAddedToCut(e))
							cs.addToCut(e);
						if (h.areFlowSendingPinsSources(e))
							continue;
						
						Assert(h.areFlowSendingPinsSourceReachable(e) || reach_and_settle || (!FlowAlgorithm::same_traversal_as_grow_assimilated && h.areAllPinsSourceReachable(e)));
#ifndef NDEBUG
						if (!FlowAlgorithm::same_traversal_as_grow_assimilated && h.areAllPinsSourceReachable(e) && !h.areFlowSendingPinsSourceReachable(e)) {
							h.reachFlowSendingPins(e);
						}
#endif
						cs.settleFlowSendingPins(e);
					}

					for (const Pin& pv : scanAllPins ? hg.pinsOf(e) : hg.pinsSendingFlowInto(e)) {
						const Node v = pv.pin;
						AssertMsg(!n.isTargetReachable(v), "Settled node " + std::to_string(v) + " is reachable from target-side");
						AssertMsg(n.isSourceReachable(v) || cs.isIsolated(v), "Settled node " + std::to_string(v) + " must be reachable from source side. It may be isolated, but in this case we don't settle it.");
						if (!n.isSource(v) && !cs.isIsolated(v)) {
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