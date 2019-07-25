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

	using PartitionBlock = uint32_t;
	static constexpr PartitionBlock invalidBlock = invalidNode;

	using HopDistance = uint32_t;

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

/*

namespace Random {
	uint32_t getSeed();
	void setSeed(uint32_t seed);
	std::mt19937& getRNG();

	bool coinToss() {
		std::uniform_int_distribution<int> dist(0,1);
		return static_cast<bool>(dist(getRNG()));
	}

	template<typename T> T randomNumber(T a = 0, T b=std::numeric_limits<T>::max())
	{
		std::uniform_int_distribution<T> dist(a,b);
		return dist(getRNG());
	}

	template<typename T> std::vector<T> randomSequence(std::size_t n, T a = 0, T b=std::numeric_limits<T>::max())
	{
		std::vector<T> res;
		std::uniform_int_distribution<T> dist(a,b);
		for (std::size_t i = 0; i < n; i++) res.push_back(dist(getRNG()));
		return res;
	}

	template<typename T> const T& sampleOne(const std::vector<T>& seq) {
		return seq[ randomNumber<std::size_t>(0, seq.size()-1) ]; }

	template<typename T> std::pair<T,T> sampleTwoDisjoint(std::vector<T>& seq) {
		std::vector<T> sample_out(2);
		std::sample(seq.begin(), seq.end(), sample_out.begin(), 2, std::mt19937{randomNumber<std::size_t>()});
		return std::make_pair(sample_out[0], sample_out[1]);
	}

}

*/
