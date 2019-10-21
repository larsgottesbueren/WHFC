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
		int multiplier = -1;
		
		static constexpr bool log = true;
		
		void flipViewDirection() {
			multiplier *= -1;
		}
		
		void clear() {
			multiplier = -1;
		}
		
		template<class CutterState, class NodeRange>
		const Node findPiercingNode(CutterState& cs, NodeRange& candidates, const NodeWeight maxBlockWeight) {
			Score maxScore;
			for (const Node u : candidates) {
				if (cs.canBeSettled(u) && cs.n.sourceWeight + hg.nodeWeight(u) <= maxBlockWeight) {		//the check is necessary for the fallback
					const Score score_u(!cs.n.isTargetReachable(u), getHopDistanceFromCut(u), Random::randomNumber(), u);
					if (maxScore < score_u)
						maxScore = score_u;
				}
			}
			return maxScore.candidate;
		}
		
	private:

		struct Score {
			bool avoidsAugmentingPaths = false;
			/* Util::InvertComparison<HopDistance> hopDistance = Util::InvertComparison<HopDistance>(maxHopDistance); */
			HopDistance hopDistance = std::numeric_limits<HopDistance>::min();
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

		HopDistance getHopDistanceFromCut(const Node x) {
												//distances of vertices on opposite side are negative --> throw away
			return useDistancesFromCut ? std::max(multiplier * distanceFromCut[x], 0) : 0;
		}

		FlowHypergraph& hg;
		
	public:
		//negative entries for original source-side, positive entries for original target-side. start counting at -1/1
		std::vector<HopDistance> distanceFromCut;
	};
}