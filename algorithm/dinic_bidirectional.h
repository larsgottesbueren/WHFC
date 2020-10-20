#include "cutter_state.h"

#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/bidirectional_distance_reachable_sets.h"

namespace whfc {
	
	class BidirectionalDinicBase {
	public:
		using ScanList = LayeredQueue<Node>;
		
		using ReachableNodes = BidirectionalDistanceReachableNodes;
		using ReachableHyperedges = BidirectionalDistanceReachableHyperedges;
		
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		
		FlowHypergraph& hg;
		LayeredQueue<Node> fqueue, bqueue;
		struct StackFrame {
			Node u;
			InHeIndex parent_he_it;
		};
		FixedCapacityStack<StackFrame> stack;
		int direction = 0, previous_cutter_state_direction = 0;
		std::vector<PinIndex> current_flow_sending_pin, current_flow_receiving_pin, current_pin;
		std::vector<InHeIndex> current_hyperedge;
		
		Flow upperFlowBound = std::numeric_limits<Flow>::max();
		
		BidirectionalDinicBase(FlowHypergraph& hg) :
			hg(hg), fqueue(hg.numNodes()), bqueue(hg.numNodes()), stack(hg.numNodes()),
			current_flow_sending_pin(hg.numHyperedges(), PinIndex::Invalid()),
			current_flow_receiving_pin(hg.numHyperedges(), PinIndex::Invalid()),
			current_pin(hg.numHyperedges(), PinIndex::Invalid()),
			current_hyperedge(hg.numNodes(), InHeIndex::Invalid())
		{
		
		}
		
		void flipViewDirection() {
			std::swap(current_flow_sending_pin, current_flow_receiving_pin);
			direction = 1 - direction;
		}
		
	};
	
	
	class BidirectionalDinic : public BidirectionalDinicBase {
	public:
		using Type = BidirectionalDinic;
		
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
		static constexpr bool log = false;
		
		BidirectionalDinic(FlowHypergraph& hg) : BidirectionalDinicBase(hg)
		{
			reset();
		}
		
		void reset() {
		
		}
		
		ScanList& getScanList() {
			return fqueue;
		}
		
		bool exhaustFlow(CutterState<Type>& cs) {
			prepare(cs);
			cs.flowValue += recycleDatastructuresFromGrowReachablePhase(cs);
			bool hasCut = false;
			while (cs.flowValue <= upperFlowBound) {
				hasCut = !buildLayeredNetwork(cs, true);
				if (hasCut || cs.flowValue >= upperFlowBound) {
					break;
				}
				else {
					cs.flowValue += augmentFlowInLayeredNetwork(cs);
				}
			}
			finish(cs);
			return hasCut;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			prepare(cs);
			Flow f = 0;
			if (buildLayeredNetwork(cs, true))
				f += augmentFlowInLayeredNetwork(cs);
			finish(cs);
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			// TODO do uni-directional only here!
			prepare(cs);
			bool found_target = buildLayeredNetwork(cs, false);
			assert(!found_target); unused(found_target);
			finish(cs);
		}
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			/*
			if (!cs.augmentingPathAvailableFromPiercing || std::none_of(cs.sourcePiercingNodes.begin(), cs.sourcePiercingNodes.end(),
																		[](const auto& sp) { return sp.isReachableFromOppositeSide; })) {
				return 0;
			}
			cs.flipViewDirection();
			resetSourcePiercingNodeDistances(cs, false);
			Flow f = augmentFlowInLayeredNetwork(cs);
			resetSourcePiercingNodeDistances(cs);
			cs.flipViewDirection();
			return f;
			 */
			return 0;		// ignore for now
		}
		
	private:
		
		void prepare(CutterState<Type>& cs) {
			previous_cutter_state_direction = cs.currentViewDirection();
			if (previous_cutter_state_direction != 0) {
				cs.flipViewDirection();
			}
			assert(cs.currentViewDirection() == 0);
		}
		
		void finish(CutterState<Type>& cs) {
			for (auto& sp : cs.sourcePiercingNodes) {
				cs.n.distance[sp.node] = cs.n.sourceSettledDistance;
			}
			for (auto& tp : cs.targetPiercingNodes) {
				cs.n.distance[tp.node] = cs.n.targetSettledDistance;
			}
			if (previous_cutter_state_direction != 0) {
				cs.flipViewDirection();
			}
			assert(cs.currentViewDirection() == previous_cutter_state_direction);
		}
		
		bool buildLayeredNetwork(CutterState<Type>& cs, const bool augment_flow) {
			unused(augment_flow);	//for debug builds only
			cs.clearForSearch();
			assert(cs.currentViewDirection() == 0);
			
			BidirectionalDistanceReachableNodes& n = cs.n;
			std::vector<DistanceT>& dist = n.distance;
			BidirectionalDistanceReachableHyperedges& h = cs.h;
			std::vector<DistanceT>& inDist = h.inDistance, outDist = h.outDistance;
			fqueue.clear(); bqueue.clear();
			bool searches_met = false;
			
			DistanceT flayer = n.s.base, blayer = n.t.upper_bound - 1;
			const DistanceT f_lb = flayer, b_ub = blayer, meeting_dist = n.meetingDistance;
			size_t fdeg = 0, bdeg = 0;
			Node source = cs.sourcePiercingNodes.begin()->node, target = cs.targetPiercingNodes.begin()->node;
			
			auto is_source_reachable = [&](Node v) {
				return (dist[v] >= f_lb && dist[v] <= flayer) || dist[v] == meeting_dist || n.isSource(v);
			};
			
			auto is_target_reachable = [&](Node v) {
				return (dist[v] <= b_ub && dist[v] >= blayer) || dist[v] == meeting_dist || n.isTarget(v);
			};
			
			auto forward_visit = [&](Node v) {
				if (!is_source_reachable(v)) {
					searches_met |= is_target_reachable(v);
					if (!searches_met) {
						dist[v] = flayer;
						fdeg += hg.degree(v);
						fqueue.push(v);
					} else {
						if (!n.isTarget(v)) {
							dist[v] = meeting_dist;
						} else {
							// TODO actually n.isTarget(target) will return false, since we overwrote its distance
							// the assertion should check whether we meet any targets actually --> so we can do that just always
							assert(v == target);
						}
					}
				}
			};
			
			auto backward_visit = [&](Node v) {
				if (!is_target_reachable(v)) {
					searches_met |= is_source_reachable(v);
					if (!searches_met) {
						dist[v] = blayer;
						bdeg += hg.degree(v);
						bqueue.push(v);
					} else {
						if (!n.isSource(v)) {
							dist[v] = meeting_dist;
						} else {
							// TODO actually n.isSource(source) will return false, since we overwrote its distance
							assert(v == source);
						}
					}
				}
			};
			
			auto finish_forward_layer = [&] {
				fqueue.finishNextLayer();
				flayer++;
			};
			
			auto finish_backward_layer = [&] {
				bqueue.finishNextLayer();
				blayer--;
			};
			
			for (auto& sp : cs.sourcePiercingNodes) {
				forward_visit(sp.node);
			}
			finish_forward_layer();
			for (auto& sp: cs.targetPiercingNodes) {
				backward_visit(sp.node);
			}
			finish_backward_layer();
			
			while (!searches_met && !fqueue.empty() && !bqueue.empty()) {
				if (fdeg < bdeg && !fqueue.empty()) {
					// advance forward search
					fdeg = 0;
					while (!fqueue.currentLayerEmpty()) {
						const Node u = fqueue.pop();
						for (InHe& inc_u : hg.hyperedgesOf(u)) {
							const Hyperedge e = inc_u.e;
							if (!h.areAllPinsSourceReachable(e)) {
								if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
									// u -> in-node(e) -> out-node(e) -> all pins | or | u -> out-node(e) -> all pins
									outDist[e] = flayer;
									for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
										forward_visit(pv.pin);
									}
								}
								if (!h.areFlowSendingPinsSourceReachable(e)) {
									// u -> in-node(e) -> all pins sending flow into e
									inDist[e] = flayer;
									for (const Pin& pv : hg.pinsSendingFlowInto(e)) {
										forward_visit(pv.pin);
									}
								}
							}
						}
					}
					finish_forward_layer();
				} else {
					// advance backward search
					bdeg = 0;
					while (!bqueue.currentLayerEmpty()) {
						const Node u = bqueue.pop();
						for (InHe& inc_u : hg.hyperedgesOf(u)) {
							const Hyperedge e = inc_u.e;
							if (!hg.isSaturated(e) || hg.flowSent(inc_u) > 0) {
								// u <- in-node(e) <- out-node(e) <- all pins | or | u <- out-node(e) <- all pins
								
							}
						}
					}
					finish_backward_layer();
					
				}
			}
			
			
			LOGGER_WN << V(searches_met) << "#BFS layers =" << (n.s.upper_bound - n.s.base);
			return searches_met;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			assert(cs.currentViewDirection() == 0);
			auto& n = cs.n;
			auto& h = cs.h;
			Flow f = 0;
			
			for (Node u : hg.nodeIDs()) {
				current_hyperedge[u] = hg.beginIndexHyperedges(u);
			}
			for (Hyperedge e : hg.hyperedgeIDs()) {
				current_pin[e] = hg.pinsNotSendingFlowIndices(e).begin();
				current_flow_sending_pin[e] = hg.pinsSendingFlowIndices(e).begin();
			}
			
			for (auto& sp : cs.sourcePiercingNodes) {
				assert(stack.empty());
				stack.push({ sp.node, InHeIndex::Invalid() });
				
				while (!stack.empty()) {
					const Node u = stack.top().u;
					Node v = invalidNode;
					InHeIndex inc_v_it = InHeIndex::Invalid();
					DistanceT req_dist = n.distance[u] + 1;
					assert(!n.isDistanceStale(u));
					assert(stack.size() + n.sourceBaseDistance() == req_dist);
					InHeIndex& he_it = current_hyperedge[u];
					for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
						InHe& inc_u = hg.getInHe(he_it);
						const Hyperedge e = inc_u.e;
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
						assert((residual > 0) == (!hg.isSaturated(e) || hg.absoluteFlowReceived(inc_u) > 0));
						const bool scanAll = req_dist == h.outDistance[e] && residual > 0;
						const bool scanFlowSending = req_dist == h.inDistance[e];
						
						if (scanFlowSending) {
							for (const PinIndex firstInvalid = hg.pinsSendingFlowIndices(e).end(); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
								const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
								if (residual + hg.absoluteFlowSent(pv) > 0 && (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist)) {
									v = pv.pin;
									inc_v_it = pv.he_inc_iter;
									break;
								}
							}
						}
						
						if (scanAll && v == invalidNode) {
							for (const PinIndex firstInvalid = hg.pinsNotSendingFlowIndices(e).end(); current_pin[e] < firstInvalid; current_pin[e]++) {
								const Pin& pv = hg.getPin(current_pin[e]);
								if (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist) {
									v = pv.pin;
									inc_v_it = pv.he_inc_iter;
									break;
								}
							}
						}
						
						if (v != invalidNode)
							break;		//don't advance hyperedge iterator
					}
					
					if (v == invalidNode) {
						assert(current_hyperedge[u] == hg.endIndexHyperedges(u));
						stack.pop();
						// Note: the iteration of u's predecessor on the stack still points to u. setting the distance to unreachable prevents the search from pushing u again.
						// It is fine to destroy the reachability datastructures, since we know that this function increases the flow.
						// An alternative method would be to advance the iteration manually, which would be hacky.
						n.distance[u] = ReachableNodes::unreachableDistance;
					}
					else {
						if (n.isTarget(v))
							f += augmentFromTarget(inc_v_it);
						else
							stack.push( { v, inc_v_it } );
					}
					
				}
			}
			assert(f > 0);
			return f;
		}
		
		
		
		Flow augmentFromTarget(InHeIndex inc_target_it) {
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
			assert(bottleneckCapacity > 0);
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