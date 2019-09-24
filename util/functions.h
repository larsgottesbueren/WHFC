#pragma once

namespace Function {
	template<typename T>
	struct Identity {
		const T& operator()(const T& x) const {
			return x;
		}
	};
	//included in the standard from C++20 onward

	template<typename RandomAccessObject>
	struct RandomAccessCallableWrapper {
		const RandomAccessObject& rao;
		RandomAccessCallableWrapper(const RandomAccessObject& rao) : rao(rao) { }
		const auto operator()(const size_t idx) const {
			return rao[idx];
		}
	};

	//use std::not_fn for negating a callable
}
