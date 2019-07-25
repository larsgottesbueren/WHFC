#pragma once

#include "../util/tiebreak.h"

namespace whfc {
	class Piercer {
	public:
		bool useDistancesFromCut = true;
		bool avoidAugmentingPaths = true;
		TieBreak tiebreak;


	};
}