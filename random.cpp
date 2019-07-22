#include "definitions.h"

namespace Random {
	static uint32_t seed;
	static std::mt19937 RNG;

	void setSeed(uint32_t _seed) {
		seed = _seed;
		RNG.seed(seed);
	}
	uint32_t getSeed() { return seed; }
	std::mt19937& getRNG() { return RNG; }

}