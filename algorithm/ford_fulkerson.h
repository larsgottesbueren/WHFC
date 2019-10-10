#pragma once

#include "cutter_state.h"
#include "../datastructure/stack.h"
#include "../datastructure/queue.h"
#include "../datastructure/bitset_reachable_sets.h"
#include "../datastructure/timestamp_reachable_sets.h"
#include "../datastructure/reachable_checker.h"

namespace whfc {
	namespace FlowCommons {
		template<typename Type>
		bool incidentToPiercingNodes(const Hyperedge e, CutterState<Type>& cs) {
			return std::any_of(cs.sourcePiercingNodes.begin(), cs.sourcePiercingNodes.end(), [&](const auto& sp) {
				return std::any_of(cs.hg.pinsOf(e).begin(), cs.hg.pinsOf(e).end(), [&](const auto& pe) {
					return pe.pin == sp.node;
				});
			});
		}
	}
	
	template<typename ScanListType, bool capacityScaling, bool alwaysSetParent = true>
	class FordFulkerson /* : public FlowAlgorithm */ {
	public:
		static constexpr bool log = true;
		
		using Type = FordFulkerson<ScanListType, capacityScaling, alwaysSetParent>;
		using ScanList = ScanListType;

		//using ReachableNodes = BitsetReachableNodes;
		//using ReachableHyperedges = BitsetReachableHyperedges;
		
		using ReachableNodes = ReachableNodesChecker;
		using ReachableHyperedges = ReachableHyperedgesChecker;
		
		using Timestamp = uint8_t;
		//using ReachableNodes = TimestampReachableNodes<Timestamp>;
		//using ReachableHyperedges = TimestampReachableHyperedges<Timestamp>;


		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		FlowHypergraph& hg;
		ScanList nodes_to_scan;
		
		struct Parent {
			InHeIndex 	parentIncidenceIterator = InHeIndex::Invalid(),
						currentIncidenceIterator = InHeIndex::Invalid();
		};
		std::vector<Parent> parent;

		explicit FordFulkerson(FlowHypergraph& hg) : hg(hg), nodes_to_scan(hg.numNodes()), parent(hg.numNodes()) {
			initalizeScalingCapacity(hg.maxHyperedgeCapacity);
		}

		static constexpr Flow DefaultInitialScalingCapacity = 1 << 24;
		static constexpr Flow ScalingCutOff = 4; //NOTE choose sensibly
		Flow initialScalingCapacity = DefaultInitialScalingCapacity;
		Flow scalingCapacity = initialScalingCapacity;


		void reset() {
			initalizeScalingCapacity(hg.maxHyperedgeCapacity);
		}
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			Flow flow = 0;
			if constexpr (alwaysSetParent) {
				if (cs.augmentingPathAvailableFromPiercing) {
					for (auto& s : cs.sourcePiercingNodes) {
						if (s.isReachableFromOppositeSide) {
							cs.flipViewDirection();
							flow += augmentFromTarget(cs.n, s.node);
							cs.flipViewDirection();
							break;	//only one path available
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
					reduceScalingCapacity();
					diff = -1;
				}
			}
			while (diff != 0) {
				diff = growWithoutScaling<true>(cs);
				flow += diff;
			}
			scalingCapacity = initialScalingCapacity;
			return flow;
		}

		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			if constexpr (capacityScaling) {
				while (scalingCapacity > ScalingCutOff) {
					Flow flow = growWithScaling(cs);
					if (flow != 0)
						return flow;
					else
						reduceScalingCapacity();
				}
			}
			Flow flow = growWithoutScaling<true>(cs);
			if (flow == 0)
				scalingCapacity = initialScalingCapacity;
			return flow;
		}

		Flow augmentFromTarget(ReachableNodes& n, const Node target) {
			Assert(n.isTarget(target));
			Flow bottleneckCapacity = maxFlow;
			Node v = target;
			while (!n.isSource(v)) {
				Parent p = parent[v];
				const Flow residual = hg.residualCapacity(hg.getInHe(p.parentIncidenceIterator), hg.getInHe(p.currentIncidenceIterator));
				bottleneckCapacity = std::min(bottleneckCapacity, residual);
				v = hg.getPin(hg.getInHe(p.parentIncidenceIterator)).pin;
			}
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
			v = target;
			while (!n.isSource(v)) {
				Parent p = parent[v];
				auto& inc_u = hg.getInHe(p.parentIncidenceIterator);
				auto& inc_v = hg.getInHe(p.currentIncidenceIterator);
				hg.routeFlow(inc_u, inc_v, bottleneckCapacity);
				v = hg.getPin(hg.getInHe(p.parentIncidenceIterator)).pin;
			}
			return bottleneckCapacity;
		}

		template<bool augment_flow>
		Flow growWithoutScaling(CutterState<Type>& cs) {
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (auto& s : cs.sourcePiercingNodes)
				nodes_to_scan.push(s.node);
			
			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
				for (InHeIndex inc_u_iter : hg.incidentHyperedgeIndices(u)) {
					const InHe& inc_u = hg.getInHe(inc_u_iter);
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
							AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
							Assert(!cs.isIsolated(v) || (augment_flow && FlowCommons::incidentToPiercingNodes(e, cs)));
							if (!n.isSourceReachable(v)) {		//don't do VD label propagation
								if constexpr (augment_flow || alwaysSetParent)
									parent[v] = { inc_u_iter, pv.he_inc_iter };
								if constexpr (augment_flow)
									if (n.isTarget(v))
										return augmentFromTarget(n, v);
								n.reach(v);		//this has to be after the return if v is target, to avoid overwriting the targetSettled timestamp with sourceReachable
								nodes_to_scan.push(v);
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
			LOGGER << "Grow with scaling " << V(scalingCapacity);
			AssertMsg(scalingCapacity > 1, "Don't call this method with ScalingCapacity <= 1. Use growWithoutScaling instead.");
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (auto& s : cs.sourcePiercingNodes)
				nodes_to_scan.push(s.node);

			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
				for (InHeIndex inc_u_iter : hg.incidentHyperedgeIndices(u)) {
					const InHe& inc_u = hg.getInHe(inc_u_iter);
					const Hyperedge e = inc_u.e;
					//can push at most flow(e) back into flow-sending pin and at most residual(e) = capacity(e) - flow(e) further flow.
					//other pins can receive at most residual(e) <= capacity(e). so checking capacity(e) < scalingCapacity is a good pruning rule
					//Note that this is not residual capacity
					if (hg.capacity(e) < scalingCapacity)
						continue;

					Flow residualCapacity = hg.absoluteFlowReceived(inc_u) + hg.residualCapacity(e);
					if (!h.areFlowSendingPinsSourceReachable(e)) {
						h.reachFlowSendingPins(e);		//Note: this is only fine because we're not using growWithScaling to determine the reachable sets!
						for (const Pin& pv : hg.pinsSendingFlowInto(e)) {
							if (residualCapacity + hg.absoluteFlowSent(pv) >= scalingCapacity) {//residual = flow received by u + residual(e) + flow sent by v
								const Node v = pv.pin;
								Assert(!cs.isIsolated(v) || FlowCommons::incidentToPiercingNodes(e, cs));
								if (!n.isSourceReachable(v)) {
									parent[v] = parent[v] = { inc_u_iter, pv.he_inc_iter };
									if (n.isTarget(v))
										return augmentFromTarget(n, v);
									n.reach(v);		//this has to be after the return if v is target, to avoid overwriting the targetSettled timestamp with sourceReachable
									nodes_to_scan.push(v);
								}
							}
						}
					}

					if (residualCapacity >= scalingCapacity && !h.areAllPinsSourceReachable(e) && (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0)) {
						h.reachAllPins(e);
						for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
							const Node v = pv.pin;
							Assert(!cs.isIsolated(v) || FlowCommons::incidentToPiercingNodes(e, cs));
							if (!n.isSourceReachable(v)) {
								parent[v] = parent[v] = { inc_u_iter, pv.he_inc_iter };
								if (n.isTarget(v))
									return augmentFromTarget(n, v);
								n.reach(v);		//this has to be after the return if v is target, to avoid overwriting the targetSettled timestamp with sourceReachable
								nodes_to_scan.push(v);
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
		
		ScanList& getScanList() {
			return nodes_to_scan;
		}
		
	private:
		void reduceScalingCapacity() {
			scalingCapacity /= 2;
		}
		
		void initalizeScalingCapacity(Flow cap) {
			cap = std::min(DefaultInitialScalingCapacity, cap);
			initialScalingCapacity = 1;
			while (2 * initialScalingCapacity <= cap) {
				initialScalingCapacity *= 2;
			}
			scalingCapacity = initialScalingCapacity;
		}
	};

	using ScalingFordFulkerson = FordFulkerson<FixedCapacityStack<Node>, true>;
	using BasicFordFulkerson = FordFulkerson<FixedCapacityStack<Node>, false>;

	using ScalingEdmondsKarp = FordFulkerson<LayeredQueue<Node>, true>;
	using BasicEdmondsKarp = FordFulkerson<LayeredQueue<Node>, false>;

	
	class DepthFirstFordFulkerson {
	public:
		static constexpr bool log = true;
		using Type = DepthFirstFordFulkerson;
		
		using ReachableNodes = ReachableNodesChecker;
		using ReachableHyperedges = ReachableHyperedgesChecker;

		using Pin = FlowHypergraph::Pin;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using InHe = FlowHypergraph::InHe;
		
		DepthFirstFordFulkerson(FlowHypergraph& hg) : hg(hg), stack(hg.numNodes()), scan_list(hg.numNodes()) {
		
		}
		
		void reset() {
		
		}
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			return 0;	//we don't know which side worked on the stack the last time, and we don't have the incoming hyperedge iterator at the target
		}
		
		Flow exhaustFlow(CutterState<Type>& cs) {
			Flow flow = 0;
			Flow diff = -1;
			while (diff != 0) {
				diff = growWithoutScaling<true>(cs);
				flow += diff;
			}
			return flow;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			return growWithoutScaling<true>(cs);
		}
		
		void growReachable(CutterState<Type>& cs) {
			growWithoutScaling<false>(cs);
		}
		
		Flow augmentFromTarget(InHeIndex inc_target_it) {
			Flow bottleneckCapacity = maxFlow;
			InHeIndex inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const StackElement& t = stack.at(stack_pointer);
				const Flow residual = hg.residualCapacity(hg.getInHe(t.out_he_it), hg.getInHe(inc_v_it));
				bottleneckCapacity = std::min(bottleneckCapacity, residual);
				inc_v_it = t.parent_he_it;
			}
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
			inc_v_it = inc_target_it;
			while (!stack.empty()) {
				const StackElement& t = stack.top();
				hg.routeFlow(hg.getInHe(t.out_he_it), hg.getInHe(inc_v_it), bottleneckCapacity);
				inc_v_it = t.parent_he_it;
				stack.pop();
			}
			return bottleneckCapacity;
		}
		
		template<bool augment_flow>
		Flow growWithoutScaling(CutterState<Type>& cs) {
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			stack.clear();
			
			for (auto& s : cs.sourcePiercingNodes) {
				Assert(stack.empty());
				stack.push({ s.node, hg.beginIndexHyperedges(s.node), InHeIndex::Invalid(), PinIndexRange::Invalid() });
				while (!stack.empty()) {
					InHeIndex& he_it = stack.top().out_he_it;
					PinIndexRange& pins_to_scan = stack.top().pins;
					const Node u = stack.top().u;
					Node v = invalidNode;
					InHeIndex inc_v_it = InHeIndex::Invalid();
					
					while (he_it < hg.endIndexHyperedges(u) && v == invalidNode) {
						if (pins_to_scan.isInvalid()) {		//start new hyperedge
							auto& inc_u = hg.getInHe(he_it);
							const Hyperedge e = inc_u.e;
							pins_to_scan = PinIndexRange();		//empty
							if (!h.areAllPinsSourceReachable(e)) {
								if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
									h.reachAllPins(e);
									pins_to_scan = hg.pinIndices(e);
								}
								else if (!h.areFlowSendingPinsSourceReachable(e)) {
									h.reachFlowSendingPins(e);
									pins_to_scan = hg.pinsSendingFlowIndices(e);
								}
							}
						}
						
						for ( ; v == invalidNode && !pins_to_scan.empty(); pins_to_scan.advance_begin()) {
							if (!n.isSourceReachable(hg.getPin(pins_to_scan.begin()).pin)) {
								v = hg.getPin(pins_to_scan.begin()).pin;
								inc_v_it = hg.getPin(pins_to_scan.begin()).he_inc_iter;
							}
						}
						
						if (v == invalidNode) {
							he_it++;
							pins_to_scan.invalidate();
						}
					}
					
					if (v == invalidNode)
						stack.pop();
					else {
						AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
						Assert(!cs.isIsolated(v) || (augment_flow && FlowCommons::incidentToPiercingNodes(hg.getInHe(he_it).e, cs)));
						if constexpr (augment_flow)
							if (n.isTarget(v))
								return augmentFromTarget(inc_v_it);
						n.reach(v);
						stack.push({ v, hg.beginIndexHyperedges(v), inc_v_it, PinIndexRange::Invalid() });
					}
				}
			}
			return 0;
		}

		struct StackElement {
			Node u;
			InHeIndex out_he_it;
			InHeIndex parent_he_it;
			PinIndexRange pins;
		};
		
		FlowHypergraph& hg;
		FixedCapacityStack<StackElement> stack;
		
		using ScanList = LayeredQueue<Node>;
		ScanList scan_list;
		ScanList& getScanList() {
			return scan_list;
		}
		
	};
	
	 
}