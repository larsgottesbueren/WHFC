#pragma once

#include "../definitions.h"
#include "../util/random.h"
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
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;
			
			multiplier = cs.currentViewDirection() == 0 ? -1 : 1;
			
			Assert(cs.n.sourceWeight == cs.n.sourceReachableWeight);
			Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
			cs.cleanUpCut();
			cs.verifyCutPostConditions();
			
			if (!cs.mostBalancedCutMode) {
				cs.cleanUpBorder();
				Score first_try = checkAllCandidates(cs.borderNodes.sourceSideBorder);
				if (first_try.candidate != invalidNode)
					return first_try.candidate;
				auto allNodes = hg.nodeIDs();
				Score second_try = checkAllCandidates(allNodes);
				return second_try.candidate;
			}
			else {
				return selectRandomAAPNodeAndRemoveNonAAPNodes();
			}
		}
		
	private:
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
		
		Node selectRandomAAPNodeAndRemoveNonAAPNodes() {
			std::vector<Node>& b = cs.borderNodes.sourceSideBorder;
			while (!b.empty()) {
				uint32_t index = Random::randomNumber(0, b.size() - 1);
				const Node u = b[index];
				
				// remove u
				b[index] = b.back();
				b.pop_back();
				
				if (cs.canBeSettled(u) && !cs.n.isTargetReachable(u))
					return u;
			}
			return invalidNode;
		}
		
		template<class NodeRange>
		Score checkAllCandidates(NodeRange& candidates) {
			Score maxScore;
			for (const Node u : candidates) {
				if (cs.canBeSettled(u) && cs.n.sourceWeight + hg.nodeWeight(u) <= cs.maxBlockWeight) {		//the canBeSettled(u) check is necessary for the fallback
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