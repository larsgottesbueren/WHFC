#pragma once

#include <vector>
#include <algorithm>

namespace whfc {
namespace util {
	/*
	 * Removes all elements which fulfill p and 1) copies to new vector 2) removes them in-place
	 */
	
	template<typename Container, typename Predicate>
	auto remove_if_copy(const Container& C, const Predicate& p) {
		std::vector<typename Container::value_type> out;
		for (const auto& x : C)
			if (!p(x))
				out.push_back(x);
		return out;
	}

	template<typename Container, typename Predicate>
	void remove_if_inplace(Container& C, const Predicate& p) {
		auto new_end = std::remove_if(C.begin(), C.end(), p);
		C.erase(new_end, C.end());
	}

}
}
