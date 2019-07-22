#pragma once

namespace Function {
	template<typename T> struct Identity { const T& operator()(const T& x) const  { return x; } };
}
