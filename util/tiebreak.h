#pragma once

#include "random.h"

namespace whfc {
	class TieBreak {
	public:
		enum Option {
			Random, First, Last
		};
		Option option = Option::Random;

		bool acceptEqual() {
			//honestly. there is no harm in a nice switch statement here. no need for inheritance with things so local :D
			// For option parsing you would need an enum and switch statement anyways
			switch (option) {
				case Option::Random : return Random::coinToss();
				case Option::First : return false;
				case Option::Last : return true;
				default : throw std::runtime_error("You idiot apparently forgot to include an option in the switch statement.");
			}
		}
	};
}