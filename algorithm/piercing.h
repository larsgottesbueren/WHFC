#pragma once

#include "../definitions.h"
#include "../util/random.h"
#include "../util/comparison.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/border.h"

namespace whfc {
	class Piercer {
	public:
		bool useDistancesFromCut = false;
		bool avoidAugmentingPaths = true;
		std::vector<HopDistance> distanceFromCut;


		bool distancesFromCutAvailable() const { return useDistancesFromCut && distanceFromCut.size() == hg.numNodes(); }

		template<class ReachableNodes>
		const Node findPiercingNode(ReachableNodes& n, NodeBorder& border) {
			Score maxScore;
			for (const Node u : border.sourceSideBorder) {
				const Score score_u(doesNodeAvoidAugmentingPath(!n.isTargetReachable(u)), getHopDistanceFromCut(u), Random::randomNumber(), u);
				if (maxScore < score_u)
					maxScore = score_u;
			}
			return maxScore.candidate;
		}


	private:

		struct Score {
			using IHD = Util::InvertComparison<HopDistance>;
			bool avoidsAugmentingPaths = false;
			IHD hopDistance = IHD(maxHopDistance);
			uint32_t randomScore = 0;
			Node candidate = invalidNode;

			explicit Score() = default;
			Score(bool aap, HopDistance hd, uint32_t rs, Node cd) : avoidsAugmentingPaths(aap), hopDistance(hd), randomScore(rs), candidate(cd) { }

			bool operator<(const Score& o) const {
				return std::tie(avoidsAugmentingPaths,hopDistance,randomScore,candidate) < std::tie(o.avoidsAugmentingPaths,hopDistance,randomScore,candidate);
			}
		};


		HopDistance getHopDistanceFromCut(const Node x) {
			return distancesFromCutAvailable() ? distanceFromCut[x] : 0;
		}

		bool doesNodeAvoidAugmentingPath(const bool does_it) {
			return avoidAugmentingPaths ? does_it : false;
		}


		FlowHypergraph& hg;
	};
}