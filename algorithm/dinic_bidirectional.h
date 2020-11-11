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

		static constexpr bool lazy_iterators = false;

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
				if (hasCut || cs.flowValue >= upperFlowBound) {
					break;
				}
				else {
					cs.flowValue += augmentFlowInLayeredNetwork(cs);
				}
			}
			finish(cs);
			timer.report(std::cout);
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
		
		
		bool buildLayeredNetwork(CutterState<Type>& cs, const bool augment_flow, const bool only_forward) {
			unused(augment_flow);	//for debug builds only
			//cs.clearForSearch();
			assert(cs.currentViewDirection() == 0);

			timer.start("BFS");

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

			size_t fvis_nodes = 0, bvis_nodes = 0, fvis_edges = 0, bvis_edges = 0, fvis_edges_fs = 0, bvis_edges_fr = 0, fvis_pins = 0, bvis_pins = 0;
			size_t intersection_size = 0;

			assert(std::none_of(dist.begin(), dist.end(), [&](auto x) { return x == meeting_dist || (x > flayer && x < blayer); }));
			
			// make piercing nodes single entries. we've stopped piercing multiple nodes a long time ago
			//Node source = cs.sourcePiercingNodes.front().node, target = cs.targetPiercingNodes.front().node;
			
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

			auto forward_visit = [&](Node v) {
				assert(!n.isTarget(v));		// we overwrote dist[target] and assume we can only find piercing nodes of the opposite side
				if (!is_source_reachable(v)) {
					if (is_target_reachable(v)) {
						searches_met = true;
						dist[v] = meeting_dist;
						intersection_size++;
					} else if (!searches_met) {
						fvis_nodes++;
						dist[v] = flayer;
						fdeg += hg.degree(v);
						fqueue.push(v);
						if constexpr (lazy_iterators) { current_hyperedge[v] = hg.beginIndexHyperedges(v); }
					}
				}
			};
			
			auto backward_visit = [&](Node v) {
				assert(!n.isSource(v)); 	// we overwrote dist[source] and assume we can only find piercing nodes of the opposite side
				if (!is_target_reachable(v)) {
					if (is_source_reachable(v)) {
						searches_met = true;
						dist[v] = meeting_dist;
						intersection_size++;
					} else if (!searches_met) {
						bvis_nodes++;
						dist[v] = blayer;
						bdeg += hg.degree(v);
						bqueue.push(v);
						if constexpr (lazy_iterators) { current_hyperedge[v] = hg.beginIndexHyperedges(v); }
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

			while (!searches_met && !fqueue.empty() && !bqueue.empty()) {
				if (fdeg < bdeg && !fqueue.empty()) {
					// advance forward search
					fdeg = 0;
					while (!fqueue.currentLayerEmpty()) {
						const Node u = fqueue.pop();
						for (InHe& inc_u : hg.hyperedgesOf(u)) {
							const Hyperedge e = inc_u.e;
							if (!are_all_pins_source_reachable(e)) {
								assert(!are_all_pins_target_reachable(e));
								if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
									// u -> in-node(e) -> out-node(e) -> all pins | or | u -> out-node(e) -> all pins
									outDist[e] = flayer;
									if constexpr (lazy_iterators) { current_pin[e] = hg.beginIndexPinsNotSendingFlow(e); }
									fvis_edges++;
									fvis_pins += hg.pinsNotSendingFlowInto(e).size();
									for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
										forward_visit(pv.pin);
									}
								}
								if (!are_flow_sending_pins_source_reachable(e)) {
									// u -> in-node(e) -> all pins sending flow into e
									inDist[e] = flayer;
									if constexpr (lazy_iterators) { current_flow_sending_pin[e] = hg.beginIndexPinsSendingFlow(e); }
									fvis_edges_fs++;
									fvis_pins += hg.pinsSendingFlowInto(e).size();
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
							if (!are_all_pins_target_reachable(e)) {
								assert(!are_all_pins_source_reachable(e));
								if (!hg.isSaturated(e) || hg.flowSent(inc_u) > 0) {
									// u <- in-node(e) <- out-node(e) <- all pins | or | u <- out-node(e) <- all pins
									inDist[e] = blayer;
									bvis_edges++;
									bvis_pins += hg.pinsNotReceivingFlowFrom(e).size();
									if constexpr (lazy_iterators) { if (hg.isSaturated(e)) { current_flow_sending_pin[e] = hg.beginIndexPinsSendingFlow(e); } }
									// in the DFS:
									// scan all pins if !hg.isSaturated(e)
									// scan flow sending pins if flow_sent(inc_u) > 0
									for (const Pin& pv : hg.pinsNotReceivingFlowFrom(e)) {
										backward_visit(pv.pin);
									}
									
								}
								if (!are_flow_receiving_pins_target_reachable(e)) {
									// u <- out-node(e) <- all pins receiving flow from e
									outDist[e] = blayer;
									bvis_edges_fr++;
									if constexpr (lazy_iterators) { current_pin[e] = hg.beginIndexPinsSendingFlow(e); }
									bvis_pins += hg.pinsReceivingFlowFrom(e).size();
									// in the DFS: scan all pins
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

			timer.stop("BFS");

			//LOGGER << V(searches_met) << "#flayers =" << (n.s.upper_bound - n.s.base) << "#blayers =" << (n.t.upper_bound - n.t.base) << V(intersection_size);
			//LOGGER << V(f_lb) << V(flayer) << V(b_ub) << V(blayer);
			LOGGER 	<< "#flayers=" << (n.s.upper_bound - n.s.base) << "#blayers=" << (n.t.upper_bound - n.t.base)
					<< V(intersection_size) << V(fvis_pins) << V(bvis_pins)
					<< V(fvis_nodes) << V(fvis_edges) << V(fvis_edges_fs)
					<< V(bvis_nodes) << V(bvis_edges) << V(bvis_edges_fr);
			return searches_met;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			assert(cs.currentViewDirection() == 0);
			auto& n = cs.n;
			auto& h = cs.h;
			
			const DistanceT meeting_dist = n.s.base - 1;
			const Node target = cs.targetPiercingNodes.front().node;
			Flow f = 0;

			if constexpr (!lazy_iterators) {
				timer.start("Init Iterators");
				for (Node u : hg.nodeIDs()) {
					current_hyperedge[u] = hg.beginIndexHyperedges(u);
				}
				for (Hyperedge e : hg.hyperedgeIDs()) {
					current_pin[e] = hg.beginIndexPinsNotSendingFlow(e);
					current_flow_sending_pin[e] = hg.beginIndexPinsSendingFlow(e);
				}
				timer.stop("Init Iterators");
			}

			timer.start("DFS");

			for (auto& sp : cs.sourcePiercingNodes) {
				assert(stack.empty());
				stack.push({ sp.node, InHeIndex::Invalid() });

				while (!stack.empty()) {
					const Node u = stack.top().pin;
					Pin next;

					// TODO simplify this bit?
					const bool at_meeting_node = n.distance[u] == meeting_dist;
					const bool in_forward_search = n.distance[u] < n.s.upper_bound;
					DistanceT req_dist_edge, req_dist_node;

					if (in_forward_search && !at_meeting_node) {
						req_dist_edge = n.distance[u] + 1;
						// -2 since we increment once after scanning the last layer of forward search, and because we want u one layer back
						req_dist_node = n.distance[u] == n.s.upper_bound - 2 ? meeting_dist : req_dist_edge;	// special case if next layer is intersection
					} else {
						// meeting node was visited in layer n.t.base+1
						req_dist_node = at_meeting_node ? n.t.base + 2 : n.distance[u] + 1;
						req_dist_edge = req_dist_node - 1;
					}

					for ( ; current_hyperedge[u] < hg.endIndexHyperedges(u); current_hyperedge[u]++) {
						const InHe& inc_u = hg.getInHe(current_hyperedge[u]);
						const Hyperedge e = inc_u.e;
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);

						if (req_dist_edge == h.inDistance[e]) {
/*							not necessary, if we're only working in forward view :)
							if (current_flow_sending_pin[e] < hg.beginIndexPinsSendingFlow(e)) {
								current_flow_sending_pin[e] = hg.beginIndexPinsSendingFlow(e);
							}
*/
							for (const PinIndex firstInvalid = hg.endIndexPinsSendingFlow(e); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
								const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
								assert(hg.absoluteFlowSent(pv) > 0);
								if (n.distance[pv.pin] == req_dist_node) {
									next = pv;
									break;
								}
							}
						}

						if (next.pin == invalidNode && residual > 0 && h.outDistance[e] == req_dist_edge) {
							// eliminating the forwardView() checks helps --> rework flow_hypergraph.h for a version with just bidirectional search later, after testing other optimizations
							for (const PinIndex firstInvalid = hg.endIndexPinsNotSendingFlow(e); current_pin[e] < firstInvalid; current_pin[e]++) {
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

			timer.stop("DFS");

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

		TimeReporter timer;

	};
}