#pragma once

#include "../definitions.h"
#include "cutter_state.h"

namespace whfc {

	template<class FlowAlgorithm>
	class Piercer {
	public:

		explicit Piercer(FlowHypergraph& hg, CutterState<FlowAlgorithm>& cs) : hg(hg), cs(cs) { }

		bool findPiercingNode() {
			if (cs.notSettledNodeWeight() == 0)
				return false;

			NodeBorder* border = cs.side_to_pierce == 0 ? &cs.borderNodes.sourceSide : &cs.borderNodes.targetSide;
			cs.clearPiercingNodes();
			size_t num_piercing_nodes = 0;
			const bool add_all_unreachables = cs.addingAllUnreachableNodesDoesNotChangeHeavierBlock() && !cs.mostBalancedCutMode;

			static constexpr bool log = true;

			for (Index i = 0; i < 2; ++i) {
				HopDistance& dist = border->maxOccupiedBucket[i];
				const size_t max_num_piercing_nodes = (i == 0 || cs.mostBalancedCutMode) ? 1 : estimateMaxNumPiercingNodes();

				for ( ; dist >= border->minOccupiedBucket[i]; --dist) {
					NodeBorder::Bucket& bucket = border->buckets[dist][i];

					if (false && i == NodeBorder::not_reachable_bucket_index && add_all_unreachables) {
						// add all unreachable border nodes to speed up the process (we're going to add all unreachable nodes anyway)
						for (Node candidate : bucket) {
							if (cs.isNonTerminal(candidate) && settlingDoesNotExceedMaxWeight(candidate)) {
								if (!cs.reachableFromSideNotToPierce(candidate)) {
									cs.addPiercingNode(candidate);
									num_piercing_nodes++;
								} else {
									border->insertIntoBucket(candidate, NodeBorder::reachable_bucket_index, dist);
								}
							}
						}
						bucket.clear();
					} else {
						// the old random, lazy-clear method. except we might do more than one node
						while (!bucket.empty()) {
							size_t it = cs.rng.randomIndex(0, bucket.size() - 1);
							Node candidate = bucket[it];
							bucket[it] = bucket.back();
							bucket.pop_back();

							if (cs.mostBalancedCutMode) {
								// track for reset
								border->removed_during_most_balanced_cut_mode[i].push_back(candidate);
							}

							if (cs.isNonTerminal(candidate) && settlingDoesNotExceedMaxWeight(candidate)) {
								if (i != NodeBorder::not_reachable_bucket_index || !cs.reachableFromSideNotToPierce(candidate)) {
									cs.addPiercingNode(candidate);
									if (++num_piercing_nodes >= max_num_piercing_nodes) {
										return true;
									}
									// TODO maybe restrict adding multiple nodes to one distance bucket at a time?
								} else if (!cs.mostBalancedCutMode) {
									// node got reachable --> move to other bucket. (no need to move if it can't be pierced in the future)
									border->insertIntoBucket(candidate, NodeBorder::reachable_bucket_index, dist);
								}
							}
						}
					}
				}

				border->clearBuckets(i);

				if (num_piercing_nodes > 0) {
					return true;
				} else if (i == NodeBorder::not_reachable_bucket_index && !cs.mostBalancedCutMode && cs.rejectPiercingIfAugmenting()) {
					// in mbc mode there can be unreachable nodes in the 2nd bucket
					return false;
				}
			}

			if (cs.rejectPiercingIfAugmenting()) {
				return false;
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

			if (p != invalidNode) {
				cs.addPiercingNode(p);
				return true;
			} else {
				return false;
			}
		}

		void reset() {
			piercingFallbacks = { 0, 0 };
		}

	private:

		size_t estimateMaxNumPiercingNodes() const {
			// TODO implement...
			return 1;
		}

		bool isCandidate(const Node u) const {
			return cs.isNonTerminal(u) && settlingDoesNotExceedMaxWeight(u);
		}

		bool settlingDoesNotExceedMaxWeight(const Node u) const {
			return (cs.side_to_pierce == 0 ? cs.source_weight : cs.target_weight) + hg.nodeWeight(u) <= cs.maxBlockWeight(cs.side_to_pierce);
		}

		FlowHypergraph& hg;
		CutterState<FlowAlgorithm>& cs;

		static constexpr uint32_t max_random_score = 1 << 25;

		std::array<int, 2> piercingFallbacks = { 0, 0 };
		static constexpr int piercingFallbackLimitPerSide = 3;
	};
}
