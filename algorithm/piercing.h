#pragma once

#include "../definitions.h"
#include "../util/comparison.h"
#include "cutter_state.h"

namespace whfc {
	
	template<class FlowAlgorithm>
	class Piercer {
	public:
		explicit Piercer(FlowHypergraph& hg, CutterState<FlowAlgorithm>& cs) : hg(hg), cs(cs), distanceFromCut(hg.numNodes(), 0) { }

		bool useDistancesFromCut = false;
		int multiplier = -1;
		
		static constexpr bool log = true;
		
		void clear() {
			multiplier = -1;
		}
		
		const Node findPiercingNode() {
			Assert(cs.hasCut);
			Assert(cs.n.sourceWeight == cs.n.sourceReachableWeight);
			Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
			LOGGER << "piercing";
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;
			multiplier = cs.currentViewDirection() == 0 ? -1 : 1;
			LOGGER << V(multiplier);
			if (!cs.mostBalancedCutMode) {
				cs.borderNodes.sourceSide.cleanUp([&](const Node& x) { return !cs.canBeSettled(x); });
				Score first_try = multiCriteriaCandidateCheck(cs.borderNodes.sourceSide.persistent_entries);
				if (first_try.candidate != invalidNode)
					return first_try.candidate;
				auto allNodes = hg.nodeIDs();
				Score second_try = multiCriteriaCandidateCheck(allNodes);
				return second_try.candidate;
			}
			else {
				while (!cs.borderNodes.sourceSide.empty()) {
					Node p = cs.borderNodes.sourceSide.popRandomEntryPreferringPersistent();
					LOGGER << "piercing node" << V(p) << V(cs.n.isTargetReachable(p)) << V(cs.canBeSettled(p)) << V(hg.nodeWeight(p));
					if (!cs.n.isTargetReachable(p) && isCandidate(p))
						return p;
				}
				LOGGER << "no piercing possible";
				return invalidNode;
			}
		}
		
	private:
		bool isCandidate(const Node u) const {
			return cs.canBeSettled(u) && cs.n.sourceWeight + hg.nodeWeight(u) <= cs.maxBlockWeight;
		}
		
		struct Score {
			bool avoidsAugmentingPaths = false;
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
			return useDistancesFromCut ? std::max(multiplier * distanceFromCut[x], 0) : 0; // distances of vertices on opposite side are negative --> throw away
		}
		
		
		template<class NodeRange>
		Score multiCriteriaCandidateCheck(NodeRange& candidates) {
			Score maxScore;
			for (const Node u : candidates) {
				if (isCandidate(u)) {		//the canBeSettled(u) check is necessary for the fallback
					const Score score_u(!cs.n.isTargetReachable(u), getHopDistanceFromCut(u), Random::randomNumber(), u);
					if (maxScore < score_u)
						maxScore = score_u;
				}
			}
			return maxScore;
		}

		FlowHypergraph& hg;
		CutterState<FlowAlgorithm>& cs;
		
	public:
		//negative entries for original source-side, positive entries for original target-side. start counting at -1/1
		std::vector<HopDistance> distanceFromCut;
	};
}