#pragma once

#include <vector>
#include "../util/range.h"

template<typename T, typename queue_size_type=uint32_t>
class LayeredQueue {
private:
	std::vector<T> queue;
public:
	using size_type = queue_size_type;
	size_type layerfront, layerend, qfront, qend;
	explicit LayeredQueue(const size_type num_elements) : queue(num_elements), layerfront(0), layerend(0), qfront(0), qend(0) {}
	inline void clear() { layerfront = 0; layerend = 0; qfront = 0; qend = 0; }
	inline bool empty() const { return qfront == qend; }
	inline bool currentLayerEmpty() const { return qfront == layerend; }
	inline T pop() { return queue[qfront++]; }
	inline T previousLayerPop() { return queue[layerfront++]; }
	inline void finishNextLayer() { layerend = qend; }
	inline void push(const T x) { Assert(qend < queue.size()); queue[qend++] = x; }
	inline bool previousLayerEmpty() const { return layerfront == layerend; }
	inline T capacity() const { return static_cast<T>(queue.size()); }
	inline std::vector<T>& data() { return queue; }
	template<typename Func> inline void forAllEverContainedElements(Func f) { for (size_type i = 0; i < qend; i++) { f(queue[i]); } }
	inline const_range<T> range(size_type __begin, size_type __end) { return { queue.begin() + __begin, queue.cbegin() + __end }; }
	inline const_range<T> currentLayer() { return range(layerfront, qend); }
	inline size_type numberOfPushedElements() const { return qend; }

	inline T elementAt(const size_type pos) const { return queue[pos]; }
	inline void setTo(const size_type pos, T element) { queue[pos] = element; }

	inline T swapFrontToPositionAndPop(size_type pos) {
		std::swap(queue[pos], queue[qfront]);
		return pop();
	}

	std::vector<T> extract() { return std::move(queue); }

	template<typename URNG>
	void shuffleNextLayer(URNG&& urng) {
		std::shuffle(queue.begin() + layerfront, queue.begin() + layerend, urng);
	}

	template<bool resize = true>
	void inject(std::vector<T> external_queue, size_type num_elements) {
		queue = std::move(external_queue);
		layerfront = 0;
		layerend = 0;
		qfront = 0;
		qend = queue.size();
		if (resize)
			queue.resize(num_elements);
	}
};
