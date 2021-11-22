#pragma once

#include "../definitions.h"
#include "bitvector.h"

namespace whfc {

class NodeBorder {
public:
	NodeBorder(const size_t initialN, const std::vector<HopDistance>& dfc, const int multiplier) :
			was_added(initialN),
			buckets(10, {Bucket(), Bucket()}),
			maxOccupiedBucket({-1,-1}),
			minOccupiedBucket({0,0}),
			backupMaxOccupiedBucket({-1,-1}),
			backupMinOccupiedBucket({0,0}),
			removed_during_most_balanced_cut_mode({Bucket(), Bucket()}),
			distance(dfc),
			multiplier(multiplier)
	{

	}

	bool wasAdded(const Node u) const {
		return was_added[u];
	}

	void add(const Node u, bool is_tr) {
		assert(!mostBalancedCutMode || !is_tr);
		assert(!wasAdded(u));
		was_added.set(u);
		const HopDistance d = getDistance(u);
		is_tr |= mostBalancedCutMode;				//reuse target_reachable_bucket_index buckets for nodes inserted during mbc
		const auto i = static_cast<Index>(is_tr);
		insertIntoBucket(u, i, d);
	}

	void insertIntoBucket(const Node u, const Index i, const HopDistance d) {
		buckets[d][i].push_back(u);
		maxOccupiedBucket[i] = std::max(maxOccupiedBucket[i], d);
		minOccupiedBucket[i] = std::min(minOccupiedBucket[i], d);
	}

	void reset(const size_t newN) {
		mostBalancedCutMode = false;
		was_added.resize(newN);
		was_added.reset(0, newN);

		for (Index i = 0; i < 2; ++i) {
			clearBuckets(i);
			assert(removed_during_most_balanced_cut_mode[i].empty());
		}
		verifyBucketsAreClean();

		HopDistance maxDistance = 0;
		for (Node i(0); i < newN; ++i) {
			maxDistance = std::max(maxDistance, getDistance(i));
		}
		if (static_cast<size_t>(maxDistance) >= buckets.size()) {
			buckets.resize(static_cast<size_t>(maxDistance + 1));
		}
	}

	void resetForMostBalancedCut() {
		// remove everything that was added during most balanced cut and is still in the buckets
		for (HopDistance d = minOccupiedBucket[most_balanced_cut_bucket_index]; d <= maxOccupiedBucket[most_balanced_cut_bucket_index]; ++d) {
			for (Node u : buckets[d][most_balanced_cut_bucket_index]) {
				was_added.reset(u);
			}
			buckets[d][most_balanced_cut_bucket_index].clear();
		}

		// reinsert the non-target-reachable nodes that were removed during most balanced cut
		for (Node u : removed_during_most_balanced_cut_mode[not_reachable_bucket_index]) {
			buckets[getDistance(u)][not_reachable_bucket_index].push_back(u);
		}

		for (Node u : removed_during_most_balanced_cut_mode[most_balanced_cut_bucket_index]) {
			was_added.reset(u);
		}

		removed_during_most_balanced_cut_mode[not_reachable_bucket_index].clear();
		removed_during_most_balanced_cut_mode[most_balanced_cut_bucket_index].clear();

		maxOccupiedBucket = backupMaxOccupiedBucket;
		minOccupiedBucket = backupMinOccupiedBucket;
	}

	using Bucket = std::vector<Node>;


	void clearBuckets(const Index i) {
		for (HopDistance d = minOccupiedBucket[i]; d <= maxOccupiedBucket[i]; ++d) {
			buckets[d][i].clear();
		}
		minOccupiedBucket[i] = 0;
		maxOccupiedBucket[i] = -1;
	}

	void enterMostBalancedCutMode () {
		mostBalancedCutMode = true;
		clearBuckets(reachable_bucket_index);
		// TODO could also filter non_reachable_bucket for already reachable nodes
		backupMaxOccupiedBucket = maxOccupiedBucket;
		backupMinOccupiedBucket = minOccupiedBucket;
	}

	HopDistance getDistance(const Node u) const {
		return std::max(multiplier * distance[u], 0); // distances of vertices on opposite side are negative --> throw away
	}

	BitVector was_added;

	static constexpr Index not_reachable_bucket_index = 0, reachable_bucket_index = 1, most_balanced_cut_bucket_index = 1;
	std::vector< std::array<Bucket, 2> > buckets;

	std::array<HopDistance, 2> maxOccupiedBucket, minOccupiedBucket, backupMaxOccupiedBucket, backupMinOccupiedBucket;

	std::array<Bucket, 2> removed_during_most_balanced_cut_mode;

	const std::vector<HopDistance>& distance;

	int multiplier;
	bool mostBalancedCutMode = false;

private:

	void verifyBucketsAreClean() {
#ifndef NDEBUG
		for (auto& bb : buckets) {
			for (Bucket& b : bb) {
				assert(b.empty());
			}
		}
#endif
	}

};

class NodeBorders {
public:
	NodeBorders(const size_t initialN) : distance(initialN, 0),
										 sourceSide(initialN, distance, -1),
										 targetSide(initialN, distance, 1) { }

	void reset(const size_t newN) {
		distance.resize(newN, 0);		// resize here in case distances are not used. however, users have to resize themselves at construction time
		sourceSide.reset(newN);
		targetSide.reset(newN);
	}

	void enterMostBalancedCutMode() {
		sourceSide.enterMostBalancedCutMode();
		targetSide.enterMostBalancedCutMode();
	}

	void resetForMostBalancedCut() {
		sourceSide.resetForMostBalancedCut();
		targetSide.resetForMostBalancedCut();
	}

	std::vector<HopDistance> distance;
	NodeBorder sourceSide, targetSide;
};

}
