#pragma once

#include "../definitions.h"
#include "../util/random.h"
#include "../util/comparison.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/border.h"

namespace whfc {
	class Piercer {
	public:
		explicit Piercer(FlowHypergraph& _hg) : hg(_hg), distanceFromCut(hg.numNodes(), 0) { }

		bool useDistancesFromCut = false;
		bool avoidAugmentingPaths = true;
		int multiplier = -1;
		
		static constexpr bool log = true;
		
		void flipViewDirection() {
			multiplier *= -1;
		}
		
		void clear() {
			multiplier = -1;
		}
		
		template<class ReachableNodes>
		const Node findPiercingNode(ReachableNodes& n, const NodeBorder& border, const NodeWeight maxBlockWeight) {
			Score maxScore;
			for (const Node u : border.sourceSideBorder) {
				if (n.sourceWeight + hg.nodeWeight(u) <= maxBlockWeight) {
					const Score score_u(doesNodeAvoidAugmentingPath(!n.isTargetReachable(u)), getHopDistanceFromCut(u), Random::randomNumber(), u);
					if (maxScore < score_u)
						maxScore = score_u;
				}
			}
			
			return maxScore.candidate;
		}
		
	private:

		struct Score {
			bool avoidsAugmentingPaths = false;
			Util::InvertComparison<HopDistance> hopDistance = Util::InvertComparison<HopDistance>(maxHopDistance);
			uint32_t randomScore = 0;
			Node candidate = invalidNode;

			explicit Score() = default;
			Score(bool aap, HopDistance hd, uint32_t rs, Node cd) : avoidsAugmentingPaths(aap), hopDistance(hd), randomScore(rs), candidate(cd) { }

			bool operator<(const Score& o) const {
				auto a = std::tie(avoidsAugmentingPaths,hopDistance,randomScore);
				auto b = std::tie(o.avoidsAugmentingPaths,o.hopDistance,o.randomScore);
				return a < b || (a == b && candidate > o.candidate);
			}
			
			inline friend std::ostream& operator<<(std::ostream& out, const Score& score) noexcept {
				return out << "{ aap " << score.avoidsAugmentingPaths << " randomScore " << score.randomScore << " node " << score.candidate << " }";
			}
		};

		bool distancesFromCutAvailable() const { return useDistancesFromCut; }

		HopDistance getHopDistanceFromCut(const Node x) {
														//distances of vertices on opposite side are negative --> throw away
			return distancesFromCutAvailable() ? std::max(multiplier * distanceFromCut[x], 0) : 0;
		}

		bool doesNodeAvoidAugmentingPath(const bool does_it) {
			return avoidAugmentingPaths ? does_it : false;
		}

		FlowHypergraph& hg;
		
	public:
		//negative entries for original source-side, positive entries for original target-side. start counting at -1/1
		std::vector<HopDistance> distanceFromCut;
	};
}