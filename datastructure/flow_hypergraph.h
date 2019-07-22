#pragma once

#include "hypergraph.h"

namespace whfc {

	//TODO for integration: fill and reserve functions, to avoid reallocations
	class FlowHypergraph {
	public:
		struct Pin {
			Node pin = invalidNode;
			InHeIndex he_inc_iter;
		};

		struct InHe {	//Hyperedge Incidence
			Hyperedge e = invalidHyperedge;
			Flow flow = Flow(0);
			PinIndex pin_iter;
		};

		struct HyperedgeData {
			PinIndex first_out = PinIndex(0);
			Flow flow = Flow(0);
			Flow capacity = Flow(0);
		};

		struct NodeData {
			HyperedgeIndex first_out = HyperedgeIndex(0);
			NodeWeight weight = NodeWeight(0);
		};

		using PinRange = mutable_range<std::vector<Pin>>;
		using PinIterator = PinRange::iterator;
		using PinIndexRange = mutable_index_range<PinIndex>;
		using InHeRange = mutable_range<std::vector<InHe>>;
		using InHeIterator = InHeRange::iterator;
		using HyperedgeIndexRange = mutable_index_range<HyperedgeIndex>;



		inline auto nodeIDs() const { return boost::irange<Node>(Node(0), Node::fromOtherValueType(numNodes() - 1)); }
		inline auto hyperedgeIDs() const { return boost::irange<Hyperedge>(Hyperedge(0), Hyperedge::fromOtherValueType(numHyperedges() - 1)); }
		inline auto pinIndices() const { return boost::irange<PinIndex>(PinIndex(0), numPins()); }

		FlowHypergraph(std::vector<NodeWeight>& node_weights, std::vector<HyperedgeWeight>& hyperedge_weights, std::vector<PinIndex>& hyperedge_sizes, std::vector<Node>& _pins) :
				nodes(node_weights.size() + 1),
				hyperedges(hyperedge_weights.size() + 1),
				pins(_pins.size()),
				incident_hyperedges(_pins.size()),
				pins_sending_flow(),
				pins_receiving_flow(),
				total_node_weight(boost::accumulate(node_weights, NodeWeight(0))),
				total_hyperedge_weight(boost::accumulate(hyperedge_weights, HyperedgeWeight(0)))

		{
			std::size_t i = 0;
			for (const Node p : _pins) {
				pins[i++].pin = p;					//copy pins
				nodes[p + 1].first_out++;			//bucket sizes
			}
			HyperedgeIndex running_hyperedge_index(0); i = 0;
			for (NodeData& u : nodes) {
				u.first_out += running_hyperedge_index;			//prefix sum
				running_hyperedge_index += u.first_out;			//prefix sum
				u.weight = node_weights[i++];					//copy node weights
			}
			for (Hyperedge e : hyperedgeIDs()) {
				hyperedges[e].capacity = hyperedge_weights[e];
				hyperedges[e+1].first_out = hyperedges[e].first_out + hyperedge_sizes[e];		//prefix sum
				for (auto pin_it = beginIndexPins(e); pin_it != endIndexPins(e); pin_it++) {
					Pin& p = pins[pin_it];
					InHeIndex ind_he = InHeIndex::fromOtherValueType(nodes[p.pin].first_out++);							//destroy first_out temporarily and reset later
					InHe& inc_he = incident_hyperedges[ind_he];
					inc_he.e = e;
					inc_he.pin_iter = pin_it;				//set iterator for pin -> its position in the pins of the hyperedge
					p.he_inc_iter = ind_he;					//set iterator for incident hyperedge -> its position in incident_hyperedges of the node
				}
			}
			for (NodeIndex u(numNodes()-1); u > 0; u--) { nodes[u].first_out = nodes[u-1].first_out; }	//reset temporarily destroyed first_out
			nodes[0].first_out = HyperedgeIndex(0);
			PinIndex x = PinIndex(0);
			for (auto e : hyperedgeIDs()) {
				pins_sending_flow[e] = PinIndexRange(x, x);	//empty range starting at the first pin of e
				x += pinCount(e);
				pins_receiving_flow[e] = PinIndexRange(x, x);	//empty range starting at one past the last pin of e
			}
		}

		inline NodeIndex numNodes() const { return NodeIndex::fromOtherValueType(nodes.size()); }
		inline HyperedgeIndex numHyperedges() const { return HyperedgeIndex::fromOtherValueType(hyperedges.size()); }
		inline PinIndex numPins() const { return PinIndex::fromOtherValueType(pins.size()); }
		inline PinIndex pinCount(const Hyperedge e) const { return hyperedges[e+1].first_out - hyperedges[e].first_out; }
		inline HyperedgeIndex degree(const Node u) const { return nodes[u+1].first_out - nodes[u].first_out; }
		inline NodeWeight totalNodeWeight() const { return total_node_weight; }
		inline NodeWeight nodeWeight(const Node u) const { return nodes[u].weight; }
		inline HyperedgeWeight totalHyperedgeWeight() const { return total_hyperedge_weight; }

		inline HyperedgeIndex beginIndexHyperedges(Node u) const { return nodes[u].first_out; }
		inline HyperedgeIndex endIndexHyperedges(Node u) const { return nodes[u+1].first_out; }
		//interface for irange is front(), drop_front(), empty(), size(), begin(), end()
		inline decltype(auto) incidentHyperedgeIndices(const Node u) const { return boost::irange<HyperedgeIndex>(beginIndexHyperedges(u), endIndexHyperedges(u)); }
		inline InHe& getInHe(const InHeIndex ind_e) { return incident_hyperedges[ind_e]; }
		inline InHe& getInHe(const Pin& pin) { return getInHe(pin.he_inc_iter); }
		inline const InHe& getInHe(const InHeIndex ind_e) const { return incident_hyperedges[ind_e]; }
		inline const InHe& getInHe(const Pin& pin) const { return getInHe(pin.he_inc_iter); }

		inline PinIndex beginIndexPins(const Hyperedge e) const { return hyperedges[e].first_out; }
		inline PinIndex endIndexPins(const Hyperedge e) const { return hyperedges[e+1].first_out; }
		inline decltype(auto) pinIndices(const Hyperedge e) const { return boost::irange<PinIndex>(beginIndexPins(e), endIndexPins(e)); }
		inline Pin& getPin(const PinIndex ind_p) { return pins[ind_p]; }
		inline Pin& getPin(const InHe& inc_p) { return getPin(inc_p.pin_iter); }
		inline const Pin& getPin(const PinIndex ind_p) const { return pins[ind_p]; }
		inline const Pin& getPin(const InHe& inc_p) const { return getPin(inc_p.pin_iter); }


		inline InHeIterator beginHyperedges() { return incident_hyperedges.begin(); }
		inline InHeIterator endHyperedges() { return incident_hyperedges.end(); }
		inline InHeIterator beginHyperedges(const Node u) { return incident_hyperedges.begin() + nodes[u].first_out; }
		inline InHeIterator endHyperedges(const Node u) { return incident_hyperedges.begin() + nodes[u+1].first_out; }
		InHeRange hyperedgesOf(const Node u) { return InHeRange(beginHyperedges(u), endHyperedges(u)); }
		InHeRange hyperedgesInRange(const HyperedgeIndexRange hir) { return InHeRange(beginHyperedges() + hir.begin(), beginHyperedges() + hir.end()); }

		inline PinIterator beginPins() { return pins.begin(); }
		inline PinIterator endPins() { return pins.end(); }
		inline PinIterator beginPins(const Hyperedge e) { return pins.begin() + hyperedges[e].first_out; }
		inline PinIterator endPins(const Hyperedge e) { return pins.begin() + hyperedges[e+1].first_out; }
		PinRange pinsOf(const Hyperedge e) { return PinRange(beginPins(e), endPins(e)); }
		PinRange pinsInRange(const PinIndexRange pir) { return PinRange(beginPins() + pir.begin(), beginPins() + pir.end()); }
		PinRange pinsSendingFlowInto(const Hyperedge e) { return pinsInRange(pins_sending_flow[e]); }
		PinRange pinsReceivingFlowFrom(const Hyperedge e) { return pinsInRange(pins_receiving_flow[e]); }

		inline bool forwardView() const { return sends_multiplier == 1; }
		void flipViewDirection() { std::swap(pins_sending_flow, pins_receiving_flow); std::swap(sends_multiplier, receives_multiplier); }


		inline Flow capacity(const Hyperedge e) const { return hyperedges[e].capacity; }
		inline Flow flow(const Hyperedge e) const { return hyperedges[e].flow; }
		inline Flow& flow(const Hyperedge e) { return hyperedges[e].flow; }
		inline Flow residualCapacity(const Hyperedge e) const { return capacity(e) - flow(e); }
		inline bool isSaturated(const Hyperedge e) const { assert(flow(e) <= capacity(e)); return flow(e) == capacity(e); }

		inline Flow flowSent(const Flow f) const { return f * sends_multiplier; }
		//flow sent from u = getPin(inc_u.pin_iter).pin into e = inc_u.e
		inline Flow flowSent(const InHe& inc_u) const { return flowSent(inc_u.flow); }
		inline Flow absoluteFlowSent(const InHe& inc_u) const { return std::max(0, flowSent(inc_u)); }
		//inline Flow flowSent(const Pin& pin) const { return flowSent(getInHe(pin)); }

		inline Flow flowReceived(const Flow f) const { return f * receives_multiplier; }
		//flow that u = getPin(inc_u.pin_iter).pin receives from e = inc_u.e
		inline Flow flowReceived(const InHe& inc_u) const { return flowReceived(inc_u.flow); }
		inline Flow absoluteFlowReceived(const InHe& inc_u) const { return std::max(0, flowReceived(inc_u)); }
		//inline Flow flowReceived(const Pin& pin) const { return flowReceived(getInHe(pin)); }

		inline Flow residualCapacity(const InHe& inc_u, InHe& inc_v) const {
			return absoluteFlowReceived(inc_u) + absoluteFlowSent(inc_v) + residualCapacity(inc_u.e);
		}
		//inline Flow residualCapacityIfFromNodeDoesNotReceiveFlow(const Flow residualE, const Pin& v) const { return std::max(0, flowSent(v)) + residualE; }
		//inline Flow residualCapacityIfFromNodeReceivesFlow(const Flow flowUReceives, const Flow residualE, const Pin& v) const { return std::max(0, flowSent(v)) + flowUReceives + residualE; }
		//need this for scaling


		void routeFlow(InHe& inc_u, InHe& inc_v, const Flow _flow) {
			const Hyperedge e = inc_u.e;
			AssertMsg(inc_u.e == inc_v.e, "Routing flow but incident hyperedges are not the same");
			AssertMsg(_flow > 0, "Routing <= 0 flow.");
			AssertMsg(_flow <= residualCapacity(inc_u, inc_v), "Routing more flow than residual capacity");
			AssertMsg(flow(e) <= capacity(e), "Routing more flow than capacity");
			AssertMsg(std::abs(inc_u.flow) <= capacity(e), "Pin capacity violated (u)");
			AssertMsg(std::abs(inc_v.flow) <= capacity(e), "Pin capacity violated (v)");

			hyperedges[e].flow += (_flow - absoluteFlowReceived(inc_u) - absoluteFlowSent(inc_v));
			const Flow prevFlowU = inc_u.flow;
			const Flow prevFlowV = inc_v.flow;
			inc_u.flow += flowSent(_flow);
			inc_v.flow += flowReceived(_flow);

			if (flowReceived(prevFlowU) > 0 && flowSent(inc_u.flow) >= 0)	//u previously received flow and now either has none, or sends flow.
				removePinFromFlowPins(inc_u, true);
			if (flowSent(inc_u.flow) > 0 && flowSent(prevFlowU) <= 0) //u now sends flow and did not previously, thus must be inserted into pins_sending_flow
				insertPinIntoFlowPins(inc_u, false);

			if (flowSent(prevFlowV) > 0 && flowReceived(inc_v.flow) >= 0) //v previously sent flow and now either has none, or receives flow.
				removePinFromFlowPins(inc_v, false);
			if (flowReceived(inc_v.flow) > 0 && flowReceived(prevFlowV) <= 0)  //v now receives flow and did not previously, thus must be inserted into pins_receiving_flow
				insertPinIntoFlowPins(inc_v, true);

			AssertMsg(pin_is_categorized_correctly(inc_u), "Pin categorized incorrectly");
			AssertMsg(pin_is_categorized_correctly(inc_v), "Pin categorized incorrectly");
		}

	private:
		std::vector<NodeData> nodes;
		std::vector<HyperedgeData> hyperedges;
		std::vector<Pin> pins;
		std::vector<InHe> incident_hyperedges;

		//TODO get rid of the range and just store one index, if this turns out to be cache inefficient later on
		std::vector<PinIndexRange> pins_sending_flow;	//indexed by hyperedge id. gives range of pin ids/iterators sending flow to that hyperedge. grows right if forwardView = true
		std::vector<PinIndexRange> pins_receiving_flow;	//indexed by hyperedge id. gives range of pin ids/iterators receiving flow from that hyperedge. grows left if forwardView = true

		NodeWeight total_node_weight = NodeWeight(0);
		HyperedgeWeight total_hyperedge_weight = HyperedgeWeight(0);
		int sends_multiplier = 1;						//if forwardView = true, flow entering hyperedge e should be positive and flow exiting e should be negative. reverse, if forwardView = false.
		int receives_multiplier = -1;

		//this nasty stuff is implemented directly here and not in range.h, since the ranges would have to fiddle with memory ranges they don't manage.
		//this is a particular special case, which should not be accessible via a public interface
		inline PinIndexRange pins_without_flow(const Hyperedge e) const {
			return forwardView() ? PinIndexRange(pins_sending_flow[e].end(), pins_receiving_flow[e].begin()) : PinIndexRange(pins_receiving_flow[e].end(), pins_sending_flow[e].begin());
		}

		PinIndex removePinFromFlowPins(InHe& inc_u, bool flow_receiving_pins) {
			const Hyperedge e = inc_u.e;
			PinIndex it_u = inc_u.pin_iter;
			PinIndexRange& flow_pins = flow_receiving_pins ? pins_receiving_flow[e] : pins_sending_flow[e];
			Assert(!flow_pins.empty());
			Assert(flow_pins.contains(it_u));

			PinIndex it_o = (forwardView() == flow_receiving_pins) ? flow_pins.begin() : PinIndex(flow_pins.end() - 1);
			InHe& inc_o = getInHe(getPin(it_o));
			Assert(flow_receiving_pins ? flowReceived(inc_o) > 0 : flowSent(inc_o) > 0);	//ensure it_o, taken from flow_pins, actually receives or sends flow, as appropriate
			if (forwardView() == flow_receiving_pins)
				flow_pins.advance_begin();
			else
				flow_pins.retreat_end();
			std::swap(inc_u.pin_iter, inc_o.pin_iter);
			std::swap(pins[it_u], pins[it_o]);
			Assert(pins_without_flow(e).contains(it_o));
			return it_o;
		}

		PinIndex insertPinIntoFlowPins(InHe& inc_u, bool flow_receiving_pins) {
			const Hyperedge e = inc_u.e;
			PinIndex it_u = inc_u.pin_iter;
			PinIndexRange& flow_pins = flow_receiving_pins ? pins_receiving_flow[e] : pins_sending_flow[e];
			Assert(pins_without_flow(e).contains(it_u));
			PinIndex it_o = (forwardView() == flow_receiving_pins) ? PinIndex(flow_pins.begin() - 1) : flow_pins.end();
			InHe& inc_o = getInHe(getPin(it_o));
			if (forwardView() == flow_receiving_pins)
				flow_pins.retreat_begin();
			else
				flow_pins.advance_end();
			std::swap(inc_u.pin_iter, inc_o.pin_iter);
			std::swap(pins[it_u], pins[it_o]);
			Assert(flow_pins.contains(it_o));
			return it_o;
		}

		void assert_backpointers_correct(const InHe& inhe) const {
			const InHe& doubled = getInHe(getPin(inhe));
			AssertMsg(doubled.pin_iter == inhe.pin_iter, "Backpointer Pin Iter inconsistent");
			AssertMsg(doubled.e == inhe.e, "Backpointer hyperedge ID inconsistent");
			AssertMsg(doubled.flow == inhe.flow, "Backpointer Pin Iter inconsistent");
		}

		void assert_backpointers_correct(const Pin& pin) const {
			const Pin& doubled = getPin(getInHe(pin));
			AssertMsg(doubled.he_inc_iter== pin.he_inc_iter, "Backpointer HyperedgeIncidence iterator inconsistent");
			AssertMsg(doubled.pin  == pin.pin, "Backpointer Node ID inconsistent");
		}

		void sanity_check_pin_ranges(const Hyperedge e) const {
			//check left / right end of pin ranges agree with first_out
			const PinIndexRange& s = forwardView() ? pins_sending_flow[e] : pins_receiving_flow[e];
			const PinIndexRange& l = !forwardView() ? pins_sending_flow[e] : pins_receiving_flow[e];
			Assert(hyperedges[e].first_out == s.begin());
			Assert(hyperedges[e+1].first_out == l.end());

			Assert(s.begin() <= s.end());
			Assert(s.end() < l.begin());
			Assert(l.begin() <= l.end());
		}

		bool pin_is_categorized_correctly(const InHe& inc_u) {
			const Hyperedge e = inc_u.e;
			const PinIndex it_u = inc_u.pin_iter;
			bool sends = flowSent(inc_u) > 0 && pins_sending_flow[e].contains(it_u);
			bool receives = flowReceived(inc_u.flow) > 0 && pins_receiving_flow[e].contains(it_u);
			bool no_flow = inc_u.flow == 0 && pins_without_flow(e).contains(it_u);
			return 		(sends && !receives && !no_flow)
					|| 	(!sends && receives && !no_flow)
					|| 	(!sends && !receives && no_flow);
		}
	};
}