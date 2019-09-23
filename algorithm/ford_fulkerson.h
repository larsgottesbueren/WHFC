#pragma once

#include "cutter_state.h"
#include "../datastructure/stack.h"
#include "../datastructure/bitset_reachable_sets.h"
//#include "../datastructure/timestamp_reachable_sets.h"

namespace whfc {

	//template<class ScanList>
	template<bool capacityScaling, bool alwaysSetParent = true>
	class FordFulkerson /* : public FlowAlgorithm */ {
	public:
		using Type = FordFulkerson;
		using ScanList = FixedCapacityStack<Node>;

		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = TimestampReachableNodes;
		//using ReachableHyperedges = TimestampReachableHyperedges;


		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		FlowHypergraph& hg;
		ScanList nodes_to_scan;
		std::vector<InHeIndex> parent;

		static constexpr Flow InitialScalingCapacity = 1 << 24; //NOTE choose sensibly
		static constexpr Flow ScalingCutOff = 4; //NOTE choose sensibly
		Flow scalingCapacity = InitialScalingCapacity;


		Flow recycleDatastructuresFromGrowReachablePhase(CutterState <Type> &cs) {
			Flow flow = 0;
			if constexpr (alwaysSetParent) {
				if (cs.augmentingPathAvailableFromPiercing) {
					for (Node s : cs.sourcePiercingNodes) {
						if (cs.n.isTargetReachable(s)) {
							cs.flipViewDirection();
							flow += augmentFromTarget(cs.n, parent[s]);
							cs.flipViewDirection();
							break;	//no VD label propagation --> only one path
						}
					}
				}
			}
			return flow;
		}

		Flow exhaustFlow(CutterState<Type>& cs) {
			Flow flow = 0;
			flow += recycleDatastructuresFromGrowReachablePhase(cs);
			Flow diff = -1;
			if constexpr (capacityScaling) {
				while (scalingCapacity > ScalingCutOff) {
					while (diff != 0) {
						diff = growWithScaling(cs);
						flow += diff;
					}
					scalingCapacity /= 2;
					diff = -1;
				}
			}
			while (diff != 0) {
				diff = growWithoutScaling<true>(cs);
				flow += diff;
			}
			scalingCapacity = InitialScalingCapacity;
			return flow;
		}

		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			if constexpr (capacityScaling) {
				while (scalingCapacity > ScalingCutOff) {
					Flow flow = growWithScaling(cs);
					if (flow != 0)
						return flow;
					else
						scalingCapacity /= 2;
				}
			}
			return growWithoutScaling<true>(cs);
		}

		Flow augmentFromTarget(ReachableNodes& n, const InHeIndex inc_target_index) {
			Flow bottleneckCapacity = maxFlow;
			const Node target = hg.getPin(hg.getInHe(inc_target_index)).pin;

			Node v = target;
			InHeIndex inc_v_index = inc_target_index;
			while (!n.isSource(v)) {
				InHeIndex inc_u_index = parent[hg.getPin(hg.getInHe(inc_v_index)).pin];
				bottleneckCapacity = std::min(bottleneckCapacity, hg.residualCapacity(hg.getInHe(inc_u_index), hg.getInHe(inc_v_index)));
				inc_v_index = inc_u_index;
				v = hg.getPin(hg.getInHe(inc_v_index)).pin;
			}
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
			v = target;
			inc_v_index = inc_target_index;
			while (!n.isSource(v)) {
				InHeIndex inc_u_index = parent[hg.getPin(hg.getInHe(inc_v_index)).pin];
				hg.routeFlow(hg.getInHe(inc_u_index), hg.getInHe(inc_v_index), bottleneckCapacity);
				inc_v_index = inc_u_index;
				v = hg.getPin(hg.getInHe(inc_v_index)).pin;
			}
			return bottleneckCapacity;
		}

		template<bool augment_flow>
		Flow growWithoutScaling(CutterState<Type>& cs) {
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
									return augmentFromTarget(n, pv.he_inc_iter);
							AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
							if (!n.isSourceReachable(v)) {		//don't do VD label propagation
								n.reach(v);
								nodes_to_scan.push(v);
								if constexpr (augment_flow || alwaysSetParent)
									parent[v] = pv.he_inc_iter;
							}
						}
					}
				}
			}
			return 0;
		}


		/*
		 * Note: capacity scaling is implemented separately from search without capacity scaling, as capacity scaling pruning requires more memory accesses than plain search
		 */
		Flow growWithScaling(CutterState<Type>& cs) {
			AssertMsg(scalingCapacity > 1, "Don't call this method with ScalingCapacity <= 1. Use growWithoutScaling instead.");
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (Node s : cs.sourcePiercingNodes) nodes_to_scan.push(s);

			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
				for (const InHe& inc_u : hg.hyperedgesOf(u)) {
					const Hyperedge e = inc_u.e;
					//can push at most flow(e) back into flow-sending pin and at most residual(e) = capacity(e) - flow(e) further flow.
					//other pins can receive at most residual(e) <= capacity(e). so checking capacity(e) < scalingCapacity is a good pruning rule
					if (hg.capacity(e) < scalingCapacity)
						continue;

					Flow residualCapacity = hg.absoluteFlowReceived(inc_u) + hg.residualCapacity(e);
					if (!h.areFlowSendingPinsSourceReachable(e)) {
						h.reachFlowSendingPins(e);
						for (const Pin& pv : hg.pinsSendingFlowInto(e)) {
							if (residualCapacity + hg.absoluteFlowSent(pv) >= scalingCapacity) {//residual = flow received by u + residual(e) + flow sent by v
								const Node v = pv.pin;
								if (n.isTarget(v))
									return augmentFromTarget(n, pv.he_inc_iter);
								if (!n.isSourceReachable(v)) {
									n.reach(v);
									nodes_to_scan.push(v);
									parent[v] = pv.he_inc_iter;
								}
							}
						}
					}

					if (residualCapacity >= scalingCapacity && !h.areAllPinsSourceReachable(e) && (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0)) {
						h.reachAllPins(e);
						for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
							const Node v = pv.pin;
							if (n.isTarget(v))
								return augmentFromTarget(n, pv.he_inc_iter);
							if (!n.isSourceReachable(v)) {
								n.reach(v);
								nodes_to_scan.push(v);
								parent[v] = pv.he_inc_iter;
							}
						}
					}
				}
			}
			return 0;
		}

		void growReachable(CutterState<Type>& cs) {
			growWithoutScaling<false>(cs);
		}
	};

	//using PseudoDepthFirstFordFulkerson = FordFulkerson<FixedCapacityStack>;
	//using EdmondsKarp = FordFulkerson<LayeredQueue>;

	class TrueDepthFirstFordFulkerson {
	public:
		using Type = TrueDepthFirstFordFulkerson;
		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = TimestampReachableNodes;
		//using ReachableHyperedges = TimestampReachableHyperedges;

		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		FlowHypergraph& hg;

		Flow growWithoutScaling(CutterState<Type>& cs) {
			/*
			for (Node s : cs.sourcePiercingNodes) {
				const InHeIndex e_it = hg.beginIndexHyperedges(s);
				const PinIndex pin_it = hg.getInHe(e_it).pin_iter;
				stack.push( { e_it, pin_it } );
			}
			 */
			return 0;
		}

		struct StackElement {
			InHeIndex e_it;
			PinIndex pin_it;
		};

		FixedCapacityStack<StackElement> stack;


	};
}