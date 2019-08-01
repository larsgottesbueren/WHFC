#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <random>
#include <chrono>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include "logger.h"
#include "util/timer.h"
#include "util/range.h"
#include "util/tagged_integer.h"

template<typename int1, typename int2>
int1 ceil_div(int1 numerator, int2 denominator) {
	return numerator/denominator + (numerator % denominator == 0 ? 0 : 1);
}

namespace aux {
	template<typename T> std::pair<T,T> minmax(T a, T b) {
		return a < b ? std::make_pair(a,b) : std::make_pair(b,a);
	}

	template<typename T>
	void min_to(T& a, T& b) { a = std::min(a,b); }
}

namespace whfc {
							//every tagged integer takes as first template argument an integer tag, which is unique to its type
	using Node = TaggedInteger<0, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	static constexpr Node invalidNode = Node::Invalid();
	using Hyperedge = TaggedInteger<1, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	static constexpr Hyperedge invalidHyperedge = Hyperedge::Invalid();
	using Flow = int32_t;
	static constexpr Flow maxFlow = std::numeric_limits<int32_t>::max();
	using NodeWeight = TaggedInteger<3, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	using HyperedgeWeight = TaggedInteger<4, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;

	using NodeIndex = TaggedInteger<5, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	using HyperedgeIndex = TaggedInteger<6, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	using PinIndex = TaggedInteger<7, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	using HyperedgeIncidenceIndex = TaggedInteger<8, uint32_t, std::numeric_limits<uint32_t>::max(), 0>;
	using InHeIndex = HyperedgeIncidenceIndex;

	using HopDistance = uint32_t;
	static constexpr HopDistance maxHopDistance = std::numeric_limits<uint32_t>::max();

	using Index = uint32_t;
	static constexpr Index invalidIndex = std::numeric_limits<uint32_t>::max();

	struct STPair {
		std::vector<Node> s,t;
	};

	class Metrics {
	public:
		static size_t largerBlockSize(const size_t numNodes, const double eps) {
			size_t perf_balance = numNodes/2 + (numNodes % 2 == 0 ? 0 : 1);
			auto sm = static_cast<size_t>( std::lround(std::ceil((1.0+eps)/2.0 * numNodes)) );
			if (eps == 0.0) sm = perf_balance;
			return std::max(perf_balance, sm);
		}

		static size_t smallerBlockSize(const size_t numNodes, const double eps) {
			size_t perf_balance = numNodes/2;
			auto sm = static_cast<size_t>( std::lround(std::floor((1.0-eps)/2.0 * numNodes)) );
			if (eps == 0.0) sm = perf_balance;
			return std::min(perf_balance, sm);
		}

		static double imbalance(const size_t numNodes, const size_t _smallerBlockSize) {
			double achieved_eps = 1.0 - (2.0 * static_cast<double>(_smallerBlockSize) / static_cast<double>(numNodes));
			if (_smallerBlockSize == numNodes / 2) { achieved_eps = 0.0; }
			return achieved_eps;
		}
	};

}