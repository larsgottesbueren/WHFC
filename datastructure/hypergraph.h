#pragma once

#include "../definitions.h"
#include <boost/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>


namespace whfc {
	class WeightedHypergraph {
	private:
		std::vector<NodeWeight> node_weights;
		std::vector<HyperedgeWeight> hyperedge_weights;

		std::vector<HyperedgeIndex> node_first_out;
		std::vector<Hyperedge> incident_hyperedges;
		std::vector<PinIndex> hyperedge_first_out;
		std::vector<Node> pins;

		NodeWeight __total_node_weight = NodeWeight(0);
		HyperedgeWeight __total_hyperedge_weight = HyperedgeWeight(0);

		void init(std::vector<PinIndex>& hyperedge_sizes) {
			for (Node pin : pins) { node_first_out[pin + 1]++; }
			std::partial_sum(node_first_out.begin(), node_first_out.end(), node_first_out.begin());
			for (HyperedgeIndex e(1); e < numHyperedges() + 1; e++) {
				hyperedge_first_out[e] = hyperedge_first_out[e-1] + hyperedge_sizes[e-1];
				for (Node pin : pinsOf(Hyperedge(e-1)))
					incident_hyperedges[node_first_out[pin]++] = Hyperedge(e - 1);
			}
			for (size_t u = numNodes()-1; u > 0; u--) { node_first_out[u] = node_first_out[u-1]; }
			node_first_out[0] = HyperedgeIndex(0);
		}

	public:
		using node_range = const_range<std::vector<Node>>;
		using node_iterator = node_range::const_iterator;
		using hyperedge_range = const_range<std::vector<Hyperedge>>;
		using hyperedge_iterator = hyperedge_range::const_iterator;

		WeightedHypergraph() = default;

		WeightedHypergraph(std::vector<NodeWeight> _node_weights, std::vector<HyperedgeWeight> _net_weights, std::vector<PinIndex>& hyperedge_sizes, std::vector<Node> _pins) :
				node_weights(std::move(_node_weights)),
				hyperedge_weights(std::move(_net_weights)),
				node_first_out(numNodes() + 1, HyperedgeIndex(0)),
			 	incident_hyperedges(_pins.size()),
				hyperedge_first_out(hyperedge_sizes.size() + 1, PinIndex(0)),
				pins(std::move(_pins)),
				__total_node_weight(boost::accumulate(node_weights, NodeWeight(0))),
				__total_hyperedge_weight(boost::accumulate(hyperedge_weights, HyperedgeWeight(0)))
		{
			init(hyperedge_sizes);
		}

		inline NodeIndex numNodes() const { return NodeIndex::fromOtherValueType(node_weights.size()); }
		inline HyperedgeIndex numHyperedges() const { return HyperedgeIndex::fromOtherValueType(hyperedge_weights.size()); }
		inline PinIndex numPins() const { return PinIndex::fromOtherValueType(pins.size()); }

		inline NodeWeight totalNodeWeight() const { return __total_node_weight; }
		inline HyperedgeWeight totalHyperedgeWeight() const { return __total_hyperedge_weight; }

		inline PinIndex pinCount(const Hyperedge e) const { return hyperedge_first_out[e+1] - hyperedge_first_out[e]; }
		inline HyperedgeIndex degree(const Node u) const { return node_first_out[u+1] - node_first_out[u]; }

		inline hyperedge_iterator beginHyperedges(const Node u) const { return incident_hyperedges.begin() + node_first_out[u]; }
		inline hyperedge_iterator endHyperedges(const Node u) const { return incident_hyperedges.begin() + node_first_out[u+1]; }
		inline hyperedge_iterator endHyperedges() const { return incident_hyperedges.cend(); }

		inline node_iterator beginPins(const Hyperedge e) const { return pins.cbegin() + hyperedge_first_out[e]; }
		inline node_iterator endPins(const Hyperedge e) const { return pins.cbegin() + hyperedge_first_out[e+1]; }
		inline node_iterator endPins() const { return pins.cend(); }

		inline PinIndex beginIndexPins(const Hyperedge e) const { return hyperedge_first_out[e]; }
		inline PinIndex endIndexPins(const Hyperedge e) const { return hyperedge_first_out[e+1]; }
		inline HyperedgeIndex beginIndexHyperedges(Node u) const { return node_first_out[u]; }
		inline HyperedgeIndex endIndexHyperedges(Node u) const { return node_first_out[u+1]; }

		const std::vector<Node>& getPins() const { return pins; }
		const std::vector<PinIndex>& getHyperedgeFirstOut() const { return hyperedge_first_out; }

		inline hyperedge_range hyperedgesOf(const Node u) const { return hyperedge_range(beginHyperedges(u) , endHyperedges(u)); }
		inline node_range pinsOf(const Hyperedge e) { return node_range(beginPins(e), endPins(e)); }

		inline NodeWeight nodeWeight(const Node u) const { return node_weights[u]; }
		inline HyperedgeWeight hyperedgeWeight(const Hyperedge e) const { return hyperedge_weights[e]; }

		inline Node get(const PinIndex ind) const { return pins[ind]; }
		inline Hyperedge get(const HyperedgeIndex ind) const { return incident_hyperedges[ind]; }

		inline decltype(auto) nodes() const { return boost::irange<Node>(Node(0), Node::fromOtherValueType(numNodes() - 1)); }
		inline decltype(auto) hyperedges() const { return boost::irange<Hyperedge>(Hyperedge(0), Hyperedge::fromOtherValueType(numHyperedges() - 1)); }
	};

	using Hypergraph = WeightedHypergraph;
}//namespace hyper
