#pragma once

#include "../definitions.h"

namespace whfc {

class Scaling {
private:
	static constexpr Flow DefaultInitialCapacity = 1 << 24;
	Flow initialCapacity = DefaultInitialCapacity;
	Flow capacity = initialCapacity;
	Flow CutOff = 3; //NOTE choose sensibly
	bool enabled = true;
public:
	
	void reduceCapacity() {
		capacity /= 2;
	}
	
	void reset() {
		capacity = initialCapacity;
	}
	
	Flow getCapacity() const {
		return use() ? capacity : 1;
	}
	
	void initialize(Flow maxScalingCap) {
		maxScalingCap = std::min(DefaultInitialCapacity, maxScalingCap);
		initialCapacity = 1;
		while (2 * initialCapacity <= maxScalingCap) {
			initialCapacity *= 2;
		}
		capacity = initialCapacity;
	}
	
	void enable() {
		enabled = true;
	}
	
	void disable() {
		enabled = false;
	}
	
	bool use() const {
		return enabled && capacity > CutOff;
	}
};

}