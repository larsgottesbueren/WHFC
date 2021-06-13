#pragma once

#include <vector>
#include "../util/range.h"
#include <boost/range/irange.hpp>

template<typename T>
class LayeredQueue {
private:
	std::vector<T> queue;
public:
	size_t layerfront = 0, layerend = 0, qfront = 0;
	explicit LayeredQueue(const size_t num_elements = 0) {
		queue.reserve(num_elements);
	}
	//Note. Use reinitialize() if you want to keep entries in the underlying vector intact, and ensure these won't be pushed again
	inline void reinitialize(size_t x) {
		layerfront = x; layerend = x; qfront = x;
		queue.erase(queue.begin() + x, queue.end());
	}
	inline void reinitialize() { reinitialize(queueEnd()); }
	inline void clear() { reinitialize(0); }
	inline bool empty() const { return qfront == queueEnd(); }
	inline bool currentLayerEmpty() const { return qfront == layerend; }
	inline T pop() { return queue[qfront++]; }
	inline T previousLayerPop() { return queue[layerfront++]; }
	inline void finishNextLayer() { layerend = queueEnd(); }
	inline void push(const T x) { queue.push_back(x); }
	inline bool previousLayerEmpty() const { return layerfront == layerend; }
	inline T capacity() const { return static_cast<T>(queue.size()); }
	inline std::vector<T>& data() { return queue; }
	template<typename Func> inline void forAllEverContainedElements(Func f) { for (size_t i = 0; i < queueEnd(); i++) { f(queue[i]); } }
	inline const_range<std::vector<T>> range(size_t __begin, size_t __end) { return { queue.begin() + __begin, queue.cbegin() + __end }; }
	inline const_range<std::vector<T>> currentLayer() { return range(qfront, layerend); }
	inline const_range<std::vector<T>> allElements() { return range(0, queueEnd()); }
	inline size_t queueEnd() const { return queue.size(); }
	inline T popBack() {
		T back = queue.back();
		queue.pop_back();
		return back;
	}

	inline decltype(auto) currentLayerIndices() { return boost::irange<size_t>(qfront, layerend); }

	inline T elementAt(const size_t pos) const { return queue[pos]; }
	inline void setTo(const size_t pos, T element) { queue[pos] = element; }

	inline T swapFrontToPositionAndPop(size_t pos) {
		std::swap(queue[pos], queue[qfront]);
		return pop();
	}

	template<typename URBG>
	void shuffleQueue(URBG& urbg, size_t a, size_t b) {
		std::shuffle(queue.begin() + a, queue.begin() + b, urbg);
	}

	template<typename URBG>
	void shuffleQueue(URBG& urbg) {
		shuffleQueue(urbg, qfront, queueEnd());
	}

	template<typename URBG>
	void shuffleCurrentLayer(URBG& urbg) {
		shuffleQueue(urbg, qfront, layerend);
	}
};
