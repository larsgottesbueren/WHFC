#pragma once

namespace Function {
	template<typename T> struct Identity { const T& operator()(const T& x) const  { return x; } };
	//included in the standard from C++20 onward

	//use std::not_fn for negating a callable

}
