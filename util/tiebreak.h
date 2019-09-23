#pragma once

#include "random.h"

namespace whfc {
	class TieBreak {
	public:
		enum Option {
			RandomScore, Random, First, Last
		};
		Option option = Option::Random;

		bool acceptEqual() {
			//honestly. there is no harm in a nice switch statement here. no need for inheritance with things so local :D
			// For option parsing you would need an enum and switch statement anyways
			switch (option) {
				case Option::RandomScore : return replaceByRandomScore();
				case Option::Random : return Random::coinToss();
				case Option::First : return false;
				case Option::Last : return true;
				default : throw std::runtime_error("You apparently forgot to include an option in the switch statement.");
			}
		}

		uint32_t bestScore = 0;

		void reset() {
			bestScore = 0;
		}

		//replace by (unbiased) coin toss favors later elements, i.e. rejection is geometrically distributed in the order of appearance. score gives uniform distribution (assuming unique scores, which is reasonable enough)
		bool replaceByRandomScore() {
			uint32_t newScore = Random::randomNumber();
			if (bestScore < newScore) {
				bestScore = newScore;
				return true;
			}
			return false;
		}
	};
}