#pragma once

#include "../definitions.h"
#include "cutter_state.h"

namespace whfc {

	template<class FlowAlgorithm>
	class Piercer {
	public:

		explicit Piercer(FlowHypergraph& hg, CutterState<FlowAlgorithm>& cs) : hg(hg), cs(cs) { }

		Node findPiercingNode() {
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;

			NodeBorder* border = cs.side_to_pierce == 0 ? &cs.borderNodes.sourceSide : &cs.borderNodes.targetSide;
			base_weight = cs.side_to_pierce == 0 ? cs.source_weight : cs.target_weight;

			// first look for piercing node in the bucket pqs
			// 0 = not target-reachable, 1 == target-reachable or inserted during most BalancedCutMode
			for (Index reachability_bucket_type(0); reachability_bucket_type < 2; ++reachability_bucket_type) {
				HopDistance& d = border->maxOccupiedBucket[reachability_bucket_type];

				for ( ; d >= border->minOccupiedBucket[reachability_bucket_type]; --d) {
					NodeBorder::Bucket& b = border->buckets[d][reachability_bucket_type];
					while (!b.empty()) {
						Node p = cs.rng.selectAndRemoveRandomElement(b);

						if (cs.mostBalancedCutMode) {
							border->removed_during_most_balanced_cut_mode[reachability_bucket_type].push_back(p);
						}

						if (isCandidate(p)) {
							//Note: the first condition relies on not inserting target-reachable nodes during most balanced cut mode
							if (reachability_bucket_type != NodeBorder::not_target_reachable_bucket_index || !cs.reachableFromOppositeSide(p)) {
								return p;
							}

							if (!cs.mostBalancedCutMode) {
								border->insertIntoBucket(p, NodeBorder::target_reachable_bucket_index, d);
							}
						}
					}
				}

				border->clearBuckets(reachability_bucket_type);
			}

			Node p = invalidNode;

			if (piercingFallbacks[cs.side_to_pierce]++ < piercingFallbackLimitPerSide) {
				// didn't find one in the bucket PQs, so pick a random unsettled node
				uint32_t rndScore = 0;
				HopDistance d = 0;
				for (const Node u : hg.nodeIDs()) {
					if (isCandidate(u)) {
						const HopDistance dist_u = border->getDistance(u);
						if (dist_u >= d) {
							const uint32_t score_u = cs.rng.randomNumber(1, max_random_score);
							if (dist_u > d || score_u > rndScore) {
								rndScore = score_u;
								p = u;
								d = dist_u;
							}
						}
					}
				}
			}

			return p;
		}

		void reset() {
			piercingFallbacks = { 0, 0 };
		}

	private:

		bool isCandidate(const Node u) const {
			return cs.canBeSettled(u) && base_weight + hg.nodeWeight(u) <= cs.maxBlockWeight(cs.side_to_pierce);
		}

		FlowHypergraph& hg;
		CutterState<FlowAlgorithm>& cs;
		NodeWeight base_weight = 0;

		static constexpr uint32_t max_random_score = 1 << 25;

		std::array<int, 2> piercingFallbacks = { 0, 0 };
		static constexpr int piercingFallbackLimitPerSide = 3;
	};
}
