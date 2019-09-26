#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <boost/dynamic_bitset.hpp>

#include "../util/functions.h"
#include "../util/filter.h"

namespace whfc {
	using BitVector = boost::dynamic_bitset<>;

	//A - B
	template<typename T>
	std::vector<T> setDifference(const std::vector<T>& A, const std::vector<T>& B, const size_t universeSize) {
		BitVector inB(universeSize);
		for (const T& b : B)
			inB.set(b);
		return util::remove_if_copy(A, Function::RandomAccessCallableWrapper(inB));
	}
}

namespace TestOptimization {
	template<typename id_t>
	class BoolVec {
	private:
		std::vector<bool> data;
	public:
		static constexpr std::size_t bits_per_block = 64;
		explicit BoolVec(const id_t __size) : data(__size, false) {}
		inline bool operator[](const id_t idx) const { return data[idx]; }
		inline void set(const id_t idx) { data[idx] = true; }
		inline void reset(const id_t idx) { data[idx] = false; }
		inline void reset() { data.assign(data.size(), false); }
		inline bool any() const { return std::any_of(data.begin(), data.end(), Function::Identity<bool>()); }
		inline id_t count() const {
			id_t res = 0;
			for (auto& x : data)
				res += static_cast<id_t>(x);
			return res;
		}
		inline std::size_t size() const { return data.size(); }

	};

}