#pragma once

#include <vector>
#include "../util/custom_asserts.h"

template<typename T, typename index_t = uint32_t>
class FixedCapacityStack {
private:
	std::vector<T> stack;
	index_t __size;
public:
	explicit FixedCapacityStack(const index_t num_elements) : stack(num_elements), __size(0) {}
	inline void clear() { __size = 0; }
	inline bool empty() const { return __size == 0; }
	inline T pop() { Assert(!empty()); return stack[--__size]; }
	inline T top() { Assert(!empty()); return stack[__size - 1]; }
	inline void push(const T& x) { Assert(__size < stack.size()); stack[__size++] = x; }
	inline T elementAt(const index_t t) const { return stack[t]; }
	inline index_t size() { return __size; }

	inline index_t capacity() const { return static_cast<index_t>(stack.size()); }
	inline std::vector<T>& data() { return stack; }
};
