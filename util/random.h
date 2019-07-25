#pragma once


#include <algorithm>
#include <random>


//stolen from KaHyPar
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

		static void setSeed(int seed) {
			instance()._gen.seed(seed);
		}

		/*
		// returns uniformly random int from the interval [low, high]
		int getRandomInt(int low, int high) {
			return _int_dist(_gen, std::uniform_int_distribution<int>::param_type(low, high));
		}

		// returns uniformly random float from the interval [low, high)
		float getRandomFloat(float low, float high) {
			return _float_dist(_gen, std::uniform_real_distribution<float>::param_type(low, high));
		}

		float getNormalDistributedFloat(float mean, float std_dev) {
			return _norm_dist(_gen, std::normal_distribution<float>::param_type(mean, std_dev));
		}
		*/
		std::mt19937 &getGenerator() {
			return _gen;
		}

	private:
		Random() :
				_gen(),
				_bool_dist(0, 1)
		//_int_dist(0, std::numeric_limits<int>::max()),
		//_float_dist(0, 1),
		//_norm_dist(0, 1)
		{}

		~Random() = default;

		std::mt19937 _gen;
		std::uniform_int_distribution<int> _bool_dist;
		//std::uniform_int_distribution<int> _int_dist;
		//std::uniform_real_distribution<float> _float_dist;
		//std::normal_distribution<float> _norm_dist;
	};
}