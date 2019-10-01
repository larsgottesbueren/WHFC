#pragma once


#include <algorithm>
#include <random>


//adapted from KaHyPar
namespace whfc {
	class Random {
	public:
		Random(const Random &) = delete;
		Random(Random &&) = delete;
		Random &operator=(const Random &) = delete;
		Random &operator=(Random &&) = delete;

		static Random &instance() {
			static Random instance;
			return instance;
		}

		static bool coinToss() {
			return static_cast<bool>(instance()._bool_dist(instance()._gen));
		}

		static uint32_t randomNumber() {
			return instance().uint_distribution(instance()._gen);
		}

		static void setSeed(int seed) {
			instance()._gen.seed(seed);
		}

		std::mt19937& getGenerator() {
			return _gen;
		}

	protected:
		Random() :
				_gen(),
				_bool_dist(0, 1),
				uint_distribution(0, std::numeric_limits<uint32_t>::max())
		{ }

		std::mt19937 _gen;
		std::uniform_int_distribution<int> _bool_dist;
		std::uniform_int_distribution<uint32_t> uint_distribution;
	};
}