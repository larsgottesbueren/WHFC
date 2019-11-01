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
			
			return invalidNode;
		}
		
	private:

		bool isCandidate(const Node u) const {
			return cs.canBeSettled(u) && cs.n.sourceWeight + hg.nodeWeight(u) <= cs.maxBlockWeight;
		}
		
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm>& cs;
		TimeReporter& timer;
		
	};
}