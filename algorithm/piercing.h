#pragma once

#include "../definitions.h"
#include "../util/random.h"
#include "../util/comparison.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/border.h"

namespace whfc {
	class Piercer {
	public:
		explicit Piercer(FlowHypergraph& _hg) : hg(_hg) { }

		bool useDistancesFromCut = false;
		bool avoidAugmentingPaths = true;
		std::vector<HopDistance> distanceFromCut;

		static constexpr bool log = true;
		
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
		
		void clear() {
			//once distanceFromCut is implemented, this will contain code
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

		bool distancesFromCutAvailable() const { return useDistancesFromCut && distanceFromCut.size() == hg.numNodes(); }

		HopDistance getHopDistanceFromCut(const Node x) {
			//TODO multiply with -1 for nodes which were not on the current source side, in the original partition
			return distancesFromCutAvailable() ? distanceFromCut[x] : 0;
		}

		bool doesNodeAvoidAugmentingPath(const bool does_it) {
			return avoidAugmentingPaths ? does_it : false;
		}

		FlowHypergraph& hg;
	};
}