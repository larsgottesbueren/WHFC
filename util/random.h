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
			return static_cast<bool>(instance().bool_dist(instance().gen));
		}

		static uint32_t randomNumber() {
			return instance().uint_dist(instance().gen);
		}

		static uint32_t randomNumber(const uint32_t a, const uint32_t b) {
			return instance().uint_dist(instance().gen, std::uniform_int_distribution<uint32_t>::param_type(a,b));
		}
		
		template<class T>
		static T selectRandomElement(std::vector<T>& range) {
			if (range.empty())
				return T::Invalid();
			uint32_t index = randomNumber(0, range.size() - 1);
			return range[index];
		}
		
		static void setSeed(int seed) {
			instance().gen.seed(seed);
		}

		std::mt19937& getGenerator() {
			return gen;
		}

	protected:
		Random() :
				gen(),
				bool_dist(0, 1),
				uint_dist(0, std::numeric_limits<uint32_t>::max())
		{ }

		std::mt19937 gen;
		std::uniform_int_distribution<int> bool_dist;
		std::uniform_int_distribution<uint32_t> uint_dist;
	};
}