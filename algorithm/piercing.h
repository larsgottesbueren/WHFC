#pragma once

#include "../definitions.h"
#include "cutter_state.h"

namespace whfc {
	
	template<class FlowAlgorithm>
	class Piercer {
	public:
		
		explicit Piercer(FlowHypergraph& hg, CutterState<FlowAlgorithm>& cs, TimeReporter& timer) : hg(hg), cs(cs), timer(timer) { }

		const Node findPiercingNode() {
			Assert(cs.hasCut);
			Assert(cs.n.sourceWeight == cs.n.sourceReachableWeight);
			Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
			
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;
			
			// first look for piercing node in the bucket pqs
			
			NodeBorder& border = *cs.borderNodes.sourceSide;
			
			// 0 = not target-reachable, 1 == target-reachable or inserted during most BalancedCutMode
			for (Index reachability_bucket_type(0); reachability_bucket_type < 2; ++reachability_bucket_type) {
				HopDistance& d = border.maxOccupiedBucket[reachability_bucket_type];
			
				for ( ; d >= border.minOccupiedBucket[reachability_bucket_type]; --d) {
					NodeBorder::Bucket& b = border.buckets[d][reachability_bucket_type];
					while (!b.empty()) {
						Node p = Random::selectAndRemoveRandomElement(b);
				
						if (cs.mostBalancedCutMode) {
							border.removed_during_most_balanced_cut_mode[reachability_bucket_type].push_back(p);
						}
				
						if (isCandidate(p)) {
							//Note: the first condition relies on not inserting target-reachable nodes during most balanced cut mode
							if (reachability_bucket_type != NodeBorder::not_target_reachable_bucket_index || !cs.n.isTargetReachable(p)) {
								return p;
							}
							
							if (!cs.mostBalancedCutMode) {
								border.insertIntoBucket(p, NodeBorder::target_reachable_bucket_index, d);
							}
						}
					}
				}
				
				border.clearBuckets(reachability_bucket_type);
			}
			
			std::cout << "Piercing Fallback" << std::endl;
			
			// didn't find one in the bucket PQs, so pick a random unsettled node
			Node p = invalidNode;
			uint32_t rndScore = 0;
			HopDistance d = 0;
			for (const Node u : hg.nodeIDs()) {
				if (isCandidate(u)) {
					const HopDistance dist_u = cs.borderNodes.distance.getHopDistanceFromCut(u);
					if (dist_u > d) {
						const uint32_t score_u = Random::randomNumber(1, max_random_score);
						if (score_u > rndScore) {
							rndScore = score_u;
							p = u;
							d = dist_u;
						}
					}
				}
			}
			
			return p;
		}
		
	private:

		bool isCandidate(const Node u) const {
			return cs.canBeSettled(u) && cs.n.sourceWeight + hg.nodeWeight(u) <= cs.maxBlockWeight;
		}
		
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm>& cs;
		TimeReporter& timer;
		
		static constexpr uint32_t max_random_score = 1 << 25;
		
	};
}