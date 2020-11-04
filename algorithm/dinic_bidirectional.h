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
		using DistanceT = ReachableNodes::DistanceT;
		
		
		FlowHypergraph& hg;
		LayeredQueue<Node> fqueue, bqueue;
		FixedCapacityStack<Pin> stack;
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
				hasCut = !buildLayeredNetwork(cs, true, false);
				LOGGER << V(hasCut);
				if (hasCut || cs.flowValue >= upperFlowBound) {
					break;
				}
				else {
					cs.flowValue += augmentFlowInLayeredNetwork(cs);
				}
			}
			finish(cs);
			LOGGER << V(cs.flowValue);
			return hasCut;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			prepare(cs);
			Flow f = 0;
			if (buildLayeredNetwork(cs, true, false))
				f += augmentFlowInLayeredNetwork(cs);
			finish(cs);
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			prepare(cs);
			bool found_target = buildLayeredNetwork(cs, false, true);
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
		
		
		enum class STATE {
			NONE, FLOW_SENDING, ALL,
		};
		std::vector<STATE> scan_for_backward;
		
		bool buildLayeredNetwork(CutterState<Type>& cs, const bool augment_flow, const bool only_forward) {
			unused(augment_flow);	//for debug builds only
			//cs.clearForSearch();
			assert(cs.currentViewDirection() == 0);
			
			if (scan_for_backward.size() != hg.numHyperedges()) {
				scan_for_backward.resize(hg.numHyperedges(), STATE::NONE);
			}
			std::fill(scan_for_backward.begin(), scan_for_backward.end(), STATE::NONE);
			
			auto& n = cs.n;
			auto& dist = cs.n.distance;
			auto& h = cs.h;
			auto& inDist = cs.h.inDistance; auto& outDist = cs.h.outDistance;
			fqueue.clear(); bqueue.clear();
			bool searches_met = false;
			
			// base values can occur in the dist labels, upper bound values not!
			DistanceT flayer = n.s.upper_bound + 1, blayer = n.t.base - 1;
			const DistanceT f_lb = n.s.upper_bound, meeting_dist = f_lb, b_ub = blayer;
			size_t fdeg = 0, bdeg = 0;
			
			LOGGER << V(flayer) << V(f_lb) << V(meeting_dist) << V(blayer) << V(b_ub);
			assert(std::none_of(dist.begin(), dist.end(), [&](auto x) { return x == meeting_dist || (x > flayer && x < blayer); }));
			
			// make piercing nodes single entries. we've stopped piercing multiple nodes a long time ago
			Node source = cs.sourcePiercingNodes.front().node, target = cs.targetPiercingNodes.front().node;
			
			// for prototyping purposes implemented here for now, since we would have to maintain flayer in the ReachableHyperedges object
			// in the end we'll put it all into the ReachableNodes/Hyperedges object and make it generic so we can run from both sides
			// this mostly benefits the DFS since we can leverage smaller degrees of piercing nodes in later iterations if only one side grows
			// should probably also merge ReachableNodes and ReachableHyperedges into one class
			auto is_source_reachable = [&](Node v) {
				return (dist[v] >= f_lb && dist[v] <= flayer) || n.isSource(v);	// meeting_dist == f_lb --> save one comparison
			};
			auto is_target_reachable = [&](Node v) {
				return (dist[v] <= b_ub && dist[v] >= blayer) || dist[v] == meeting_dist || n.isTarget(v);
			};
			auto are_all_pins_target_reachable = [&](Hyperedge e) {
				return inDist[e] == h.targetSettledDistance || (inDist[e] <= b_ub && inDist[e] >= blayer);
			};
			auto are_flow_receiving_pins_target_reachable = [&](Hyperedge e) {
				return outDist[e] == h.targetSettledDistance || (outDist[e] <= b_ub && outDist[e] >= blayer);
			};
			auto are_all_pins_source_reachable = [&](Hyperedge e) {
				assert(outDist[e] != meeting_dist);
				return outDist[e] == h.sourceSettledDistance || (outDist[e] >= f_lb && outDist[e] <= flayer);
			};
			auto are_flow_sending_pins_source_reachable = [&](Hyperedge e) {
				assert(inDist[e] != meeting_dist);
				return inDist[e] == h.sourceSettledDistance || (inDist[e] >= f_lb && inDist[e] <= flayer);
			};
			
			size_t intersection_size = 0;
			
			auto forward_visit = [&](Node v) {
				assert(!n.isTarget(v));		// we overwrote dist[target] and assume we can only find piercing nodes of the opposite side
				if (!is_source_reachable(v)) {
					searches_met |= is_target_reachable(v);
					if (!searches_met) {
						dist[v] = flayer;
						fdeg += hg.degree(v);
						fqueue.push(v);
						// set current_hyperedge[v] = hg.beginIndexHyperedges(v);
					} else if (is_target_reachable(v)) {
						dist[v] = meeting_dist;
						intersection_size++;
					}
				}
			};
			
			auto backward_visit = [&](Node v) {
				assert(!n.isSource(v)); 	// we overwrote dist[source] and assume we can only find piercing nodes of the opposite side
				if (!is_target_reachable(v)) {
					searches_met |= is_source_reachable(v);
					if (!searches_met) {
						dist[v] = blayer;
						bdeg += hg.degree(v);
						bqueue.push(v);
						// set current_hyperedge[v] = hg.beginIndexHyperedges(v);
					} else if (is_source_reachable(v)) {
						dist[v] = meeting_dist;
						intersection_size++;
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
				dist[sp.node] = ReachableNodes::unreachableDistance;
				forward_visit(sp.node);
			}
			finish_forward_layer();
			
			if (!only_forward) {
				for (auto& tp: cs.targetPiercingNodes) {
					dist[tp.node] = ReachableNodes::unreachableDistance;
					backward_visit(tp.node);
				}
				finish_backward_layer();
			}
			
			LOGGER << V(fdeg) << V(bdeg) << V(source) << V(dist[source]) << V(target) << V(dist[target]);
			
			while (!searches_met && !fqueue.empty() && !bqueue.empty()) {		// TODO this condition is not yet compatible with the only_forward parameter --> fix it later
				if (fdeg < bdeg && !fqueue.empty()) {
					// advance forward search
					LOGGER << "fsearch" << V(fqueue.currentLayerSize()) << V(flayer) << V(fdeg) << V(blayer) << V(bdeg);
					fdeg = 0;
					while (!fqueue.currentLayerEmpty()) {
						const Node u = fqueue.pop();
						for (InHe& inc_u : hg.hyperedgesOf(u)) {
							const Hyperedge e = inc_u.e;
							if (!are_all_pins_source_reachable(e)) {
								if (are_all_pins_target_reachable(e)) {
									LOGGER << V(u) << V(e) << V(outDist[e]) << V(inDist[e]);
									for (Pin pv : hg.pinsOf(e)) {
										LOGGER << V(pv.pin) << V(hg.flowSent(pv)) << V(dist[pv.pin]);
									}
								}
								assert(!are_all_pins_target_reachable(e));
								if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
									// u -> in-node(e) -> out-node(e) -> all pins | or | u -> out-node(e) -> all pins
									outDist[e] = flayer;
									for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
										forward_visit(pv.pin);
									}
								}
								if (!are_flow_sending_pins_source_reachable(e)) {
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
					LOGGER << "bsearch" << V(bqueue.currentLayerSize()) << V(flayer) << V(fdeg) << V(blayer) << V(bdeg);
					bdeg = 0;
					while (!bqueue.currentLayerEmpty()) {
						const Node u = bqueue.pop();
						for (InHe& inc_u : hg.hyperedgesOf(u)) {
							const Hyperedge e = inc_u.e;
							if (!are_all_pins_target_reachable(e)) {
								assert(!are_all_pins_source_reachable(e));
								// version 1 --> emulate what would happen if we flipped view direction. so this is compatible with growing target_reachable
								if (!hg.isSaturated(e) || hg.flowSent(inc_u) > 0) {
									// u <- in-node(e) <- out-node(e) <- all pins | or | u <- out-node(e) <- all pins
									inDist[e] = blayer;
									// in the DFS (forward search):
									// scan all pins if !hg.isSaturated(e)
									// scan flow sending pins if flow_sent(inc_u) > 0
									if (!hg.isSaturated(e)) {
										scan_for_backward[e] = STATE::ALL;
									} else if (scan_for_backward[e] != STATE::ALL && hg.flowSent(inc_u) > 0) {
										scan_for_backward[e] = STATE::FLOW_SENDING;
									}
									
									for (const Pin& pv : hg.pinsNotReceivingFlowFrom(e)) {
										backward_visit(pv.pin);
									}
									
								}
								if (!are_flow_receiving_pins_target_reachable(e)) {
									// u <- out-node(e) <- all pins receiving flow from e
									outDist[e] = blayer;
									// in the DFS: scan all pins
									scan_for_backward[e] = STATE::ALL;
									for (const Pin& pv : hg.pinsReceivingFlowFrom(e)) {
										backward_visit(pv.pin);
									}
									
									
								}
							}
							
							
							
						}
					}
					finish_backward_layer();
					
				}
			}
			
			n.s.base = f_lb + 1;			// f_lb was set to one lower ( == meeting_dist), so we could save a comparison for the meeting_dist
			n.s.upper_bound = flayer;		// value not used
			n.t.base = blayer;				// value not used
			n.t.upper_bound = b_ub;
			
			LOGGER << V(searches_met) << "#flayers =" << (n.s.upper_bound - n.s.base) << "#blayers =" << (n.t.upper_bound - n.t.base) << V(intersection_size);
			LOGGER << V(f_lb) << V(flayer) << V(b_ub) << V(blayer);
			return searches_met;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			assert(cs.currentViewDirection() == 0);
			auto& n = cs.n;
			auto& h = cs.h;
			
			const DistanceT meeting_dist = n.s.base - 1;
			const Node target = cs.targetPiercingNodes.front().node;
			Node source = cs.sourcePiercingNodes.front().node;
			
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
					const Node u = stack.top().pin;
					Pin next;
					InHeIndex& he_it = current_hyperedge[u];
					
					const bool at_meeting_node = n.distance[u] == meeting_dist;
					const bool in_forward_search = n.distance[u] < n.s.upper_bound;
					
					if (in_forward_search && !at_meeting_node) {
						const DistanceT req_dist_edge = n.distance[u] + 1;
						
						// -2 since we increment once after scanning the last layer of forward search, and because we want u one layer back
						const DistanceT req_dist_node = n.distance[u] == n.s.upper_bound - 2 ? meeting_dist : req_dist_edge;	// special case if next layer is intersection
						
						for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
							InHe& inc_u = hg.getInHe(he_it);
							const Hyperedge e = inc_u.e;
							const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
							
							if (req_dist_edge == h.inDistance[e]) {
								for (const PinIndex firstInvalid = hg.pinsSendingFlowIndices(e).end(); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
									const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
									// TODO why are we checking again for absoluteFlowSent(pv) > 0? should always be the case, right?
									// maybe if flow sending pins are stored at the end of the range and current_flow_sending_pin[e] < hg.flowSendingPinsIndices(e).begin()
									// because some flow was routed. in this case we could raise current_flow_sending_pin[e]
									assert(hg.absoluteFlowSent(pv) > 0 || current_flow_sending_pin[e] < hg.pinsSendingFlowIndices(e).begin());
									// TODO the lookup for flowSent from the pin is expensive since it's not a scan
									// --> either store flow at both pin and edge-incidence, or avoid lookup by raising iterator!
									// we can also sync the iterator from the routeFlow function, but that incurs more dense coupling --> bad software design?
									if (residual + hg.absoluteFlowSent(pv) > 0 && n.distance[pv.pin] == req_dist_node) {
										next = pv;
										break;
									}
								}
							}
							
							if (next.pin == invalidNode && residual > 0 && h.outDistance[e] == req_dist_edge) {
								for (const PinIndex firstInvalid = hg.pinsNotSendingFlowIndices(e).end(); current_pin[e] < firstInvalid; current_pin[e]++) {
									const Pin& pv = hg.getPin(current_pin[e]);
									if (n.distance[pv.pin] == req_dist_node) {
										next = pv;
										break;
									}
								}
							}
							
							if (next.pin != invalidNode)
								break;		//don't advance hyperedge iterator
						}
						
					} else {
						// meeting node was visited in layer n.t.base+1
						const DistanceT req_dist_node = at_meeting_node ? n.t.base + 2 : n.distance[u] + 1;
						const DistanceT req_dist_edge = req_dist_node - 1;
						
						for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
							InHe& inc_u = hg.getInHe(he_it);
							const Hyperedge e = inc_u.e;
							const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
							
							if (	req_dist_edge == h.inDistance[e]
									&& (scan_for_backward[e] == STATE::FLOW_SENDING || scan_for_backward[e] == STATE::ALL)
								)
							{
								for (const PinIndex firstInvalid = hg.pinsSendingFlowIndices(e).end(); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
									const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
									if (residual + hg.absoluteFlowSent(pv) > 0 && n.distance[pv.pin] == req_dist_node) {
										next = pv;
										break;
									}
								}
							}
							
							if (next.pin == invalidNode
								&& residual > 0
								&& (req_dist_edge == h.inDistance[e] || req_dist_edge == h.outDistance[e])
								&& scan_for_backward[e] == STATE::ALL
								)
							{
								for (const PinIndex firstInvalid = hg.pinsNotSendingFlowIndices(e).end(); current_pin[e] < firstInvalid; current_pin[e]++) {
									const Pin& pv = hg.getPin(current_pin[e]);
									if (n.distance[pv.pin] == req_dist_node) {
										next = pv;
										break;
									}
								}
							}
							
							if (next.pin != invalidNode)
								break;		//don't advance hyperedge iterator
						}
						
					}
					
					
					if (next.pin == invalidNode) {
						assert(current_hyperedge[u] == hg.endIndexHyperedges(u));
						stack.pop();
						// Note: the iteration of u's predecessor on the stack still points to u. setting the distance to unreachable prevents the search from pushing u again.
						// It is fine to destroy the reachability datastructures, since we know that this function increases the flow.
						// An alternative method would be to advance the iteration manually, which would be hacky.
						n.distance[u] = ReachableNodes::unreachableDistance;
					}
					else {
						//if (n.isTarget(v))
						if (next.pin == target)
							f += augmentFromTarget(next.he_inc_iter);
						else
							stack.push(next);
					}
					
				}
			}
			LOGGER << V(f);
			assert(f > 0);
			return f;
		}
		
		
		
		Flow augmentFromTarget(InHeIndex inc_target_it) {
			Flow bottleneckCapacity = maxFlow;
			int64_t lowest_bottleneck = std::numeric_limits<int64_t>::max();
			InHeIndex inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const Pin& t = stack.at(stack_pointer);
				const Flow residual = hg.residualCapacity(hg.getInHe(current_hyperedge[t.pin]), hg.getInHe(inc_v_it));
				if (residual <= bottleneckCapacity) {
					bottleneckCapacity = residual;
					lowest_bottleneck = stack_pointer;
				}
				inc_v_it = t.he_inc_iter;	// parent
			}
			assert(bottleneckCapacity > 0);
			inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const Pin& t = stack.at(stack_pointer);
				hg.routeFlow(hg.getInHe(current_hyperedge[t.pin]), hg.getInHe(inc_v_it), bottleneckCapacity);
				inc_v_it = t.he_inc_iter;
			}
			stack.popDownTo(lowest_bottleneck);
			return bottleneckCapacity;
		}
		
	};
}