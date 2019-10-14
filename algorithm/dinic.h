
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "ford_fulkerson.h"

namespace whfc {
	
	//TODO scaling
	class Dinic {
	public:
		using Type = Dinic;

		using ReachableNodes = DistanceReachableNodes;	//TODO also write checkers for the distance sets
		using ReachableHyperedges = DistanceReachableHyperedges;
		
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		static constexpr bool log = true;
		
		FlowHypergraph& hg;
		LayeredQueue<Node> queue;
		struct StackFrame {
			Node u;
			InHeIndex parent_he_it;
		};
		FixedCapacityStack<StackFrame> stack;
		std::vector<PinIndexRange> pin_range;
		std::vector<InHeIndex> current_hyperedge;
		
		Dinic(FlowHypergraph& hg) : hg(hg), queue(hg.numNodes()), stack(hg.numNodes()),
									pin_range(hg.numHyperedges(), PinIndexRange::createEmpty()), current_hyperedge(hg.numNodes(), InHeIndex::Invalid())
		{
			
		}
		
		
		
		Flow exhaustFlow(CutterState<Type>& cs) {
			Flow f = 0;
			f += recycleDatastructuresFromGrowReachablePhase(cs);
			while (buildLayeredNetwork<true>(cs)) {
				f += augmentFlowInLayeredNetwork(cs);
			}
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			Flow f = 0;
			if (buildLayeredNetwork<true>(cs))
				f += augmentFlowInLayeredNetwork(cs);
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			cs.flipViewDirection();
			resetSourcePiercingNodeDistances(cs, false);
			//Flow f = augmentFlowInLayeredNetwork(cs);
			Flow f = 0;
			resetSourcePiercingNodeDistances(cs);
			cs.flipViewDirection();
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			bool found_target = buildLayeredNetwork<false>(cs);
			Assert(!found_target);
			resetSourcePiercingNodeDistances(cs);
		}
		
	private:
		
		void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
			for (auto& sp: cs.sourcePiercingNodes)
				cs.n.setPiercingNodeDistance(sp.node, reset);
		}
		
		template<bool augment_flow>
		bool buildLayeredNetwork(CutterState<Type>& cs) {
			LOGGER << "bfs";
			cs.clearForSearch();
			auto& n = cs.n;
			auto& h = cs.h;
			queue.clear();
			bool found_target = false;
			
			for (auto& sp : cs.sourcePiercingNodes) {
				n.setPiercingNodeDistance(sp.node, false);
				queue.push(sp.node);
				current_hyperedge[sp.node] = hg.beginIndexHyperedges(sp.node);
			}
			n.hop(); h.hop(); queue.finishNextLayer();
			
			while (!queue.empty()) {
				while (!queue.currentLayerEmpty()) {
					const Node u = queue.pop();
					for (InHe& inc_u : hg.hyperedgesOf(u)) {
						const Hyperedge e = inc_u.e;
						if (!h.areAllPinsSourceReachable(e)) {
							if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
								h.reachAllPins(e);
								Assert(n.distance[u] + 1 == h.outDistance[e]);
								pin_range[e] = hg.pinIndices(e);
							}
							else if (!h.areFlowSendingPinsSourceReachable(e)) {
								h.reachFlowSendingPins(e);
								Assert(n.distance[u] + 1 == h.inDistance[e]);
								pin_range[e] = hg.pinsSendingFlowIndices(e);
							}
							else {
								continue;
							}
							
							for (const Pin& pv : hg.pinsInRange(pin_range[e])) {
								const Node v = pv.pin;
								AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
								Assert(!cs.isIsolated(v) || (augment_flow && FlowCommons::incidentToPiercingNodes(e, cs)));
								found_target |= n.isTarget(v);
								if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
									n.reach(v);
									Assert(n.distance[u] + 1 == n.distance[v]);
									queue.push(v);
									current_hyperedge[v] = hg.beginIndexHyperedges(v);
								}
							}
						}
					}
				}
				
				n.hop(); h.hop(); queue.finishNextLayer();
			}
			
			n.lockInSourceDistance(); h.lockInSourceDistance();
			h.compareDistances(n);
			LOGGER << V(hg.numNodes()) << V(n.sourceReachableWeight) << V(n.sourceWeight);
			return found_target;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			LOGGER << "dfs";
			auto& n = cs.n;
			auto& h = cs.h;
			Flow f = 0;
			for (auto& sp : cs.sourcePiercingNodes) {
				Assert(stack.empty());
				stack.push({ sp.node, InHeIndex::Invalid() });
				
				while (!stack.empty()) {
					const Node u = stack.top().u;
					Node v = invalidNode;
					DistanceT req_dist = n.distance[u] + 1;
					Assert(!n.isDistanceStale(u));
					Assert(stack.size() + n.sourceBaseDistance() == req_dist);
					
					InHeIndex& he_it = current_hyperedge[u];
					for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
						InHe& inc_u = hg.getInHe(he_it);
						const Hyperedge e = inc_u.e;
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
						const bool scanAll = req_dist == h.outDistance[e];
						const bool scanFlowSending = req_dist == h.inDistance[e];
						if (scanAll || scanFlowSending) {	//TODO check for the appropriate ranges with the actual conditions, not the distances
							for ( ; v == invalidNode && !pin_range[e].empty(); pin_range[e].advance_begin()) {
								const Pin& pv = hg.getPin(pin_range[e].begin());
								if (residual + hg.absoluteFlowSent(pv) > 0 && (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist))
									v = pv.pin;
							}
							if (v != invalidNode)
								break;		//don't advance hyperedge iterator
						}
					}
					
					if (v == invalidNode)
						stack.pop();
					else {
						PinIndexRange pinsOfE = pin_range[ hg.getInHe(he_it).e ];
						pinsOfE.retreat_begin();
						const InHeIndex inc_v_it = hg.getPin(pinsOfE.begin()).he_inc_iter;
						if (n.isTarget(v))
							f += augmentFromTarget(cs, inc_v_it);
						else
							stack.push( { v, inc_v_it } );
					}
					
				}
			}
			Assert(f > 0);
			return f;
		}
		
		Flow augmentFromTarget(CutterState<Type>& cs, InHeIndex inc_target_it) {
			LOGGER << "found target";
			Flow bottleneckCapacity = maxFlow;
			int64_t lowest_bottleneck = std::numeric_limits<int64_t>::max();
			InHeIndex inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const StackFrame& t = stack.at(stack_pointer);
				const Flow residual = hg.residualCapacity(hg.getInHe(current_hyperedge[t.u]), hg.getInHe(inc_v_it));
				if (residual <= bottleneckCapacity) {
					bottleneckCapacity = residual;
					lowest_bottleneck = stack_pointer;
				}
				inc_v_it = t.parent_he_it;
			}
			LOGGER << V(bottleneckCapacity) << V(lowest_bottleneck);
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
			inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const StackFrame& t = stack.at(stack_pointer);
				hg.routeFlow(hg.getInHe(current_hyperedge[t.u]), hg.getInHe(inc_v_it), bottleneckCapacity);
				inc_v_it = t.parent_he_it;
			}
			stack.popDownTo(lowest_bottleneck);
			return bottleneckCapacity;
		}
		
		
	};
}