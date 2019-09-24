#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"
#include "../util/gcd.h"

namespace whfc {
	class IsolatedNodes {
	private:
		FlowHypergraph& hg;
	public:
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<HyperedgeIndex> mixedIncidentHyperedges;

		struct SummableRange {
			NodeWeight from, to;
			SummableRange(NodeWeight _from, NodeWeight _to) : from(_from), to(_to) { AssertMsg(_from != NodeWeight::Invalid() && _to != NodeWeight::Invalid(), "Invalid range"); }
			bool inRange(const NodeWeight w) const { return from <= w && w <= to; }
			bool operator<(const SummableRange& o) const { return std::tie(from, to) < std::tie(o.from, o.to); }
			bool operator==(const SummableRange& o) const { return from == o.from && to == o.to; }
			bool operator!=(const SummableRange& o) const { return !operator==(o); }
		};

	private:

		//If the subset sum approach gets too slow, there are two alternatives:
		//	1) try scaling the weights of all nodes in the flow network --> more locality in the DP table. if GCD(weights) > 1 then using ranges for balance checking doesn't help at all
				//as the ranges with rescaled weights would consist of single elements
				//we would need to implement two versions depending on GCD(weights) { == 1, > 1 }
		//	2) just assign the nodes ad-hoc, after growAssimilated. shouldn't be too bad. works well enough with AAP
		NodeWeight maxSubsetSumWeight = NodeWeight(0);


		struct TableEntry {
			Node node;
			Index sumsIndex;
			TableEntry() : node(invalidNode), sumsIndex(invalidIndex) { }
			bool summable() const {
						//sumsIndex = 0 is reserved for the first range.
				AssertMsg(sumsIndex == invalidIndex || node != invalidNode || sumsIndex == 0, "Summable but no node");
				return sumsIndex != invalidIndex;
			}
		};

		std::vector<TableEntry> DPTable;

		std::vector<SummableRange> sumRanges;
		std::vector<SummableRange> nextSumRanges;
		bool newSumAvailable = true;

		std::vector<Node> nodesNotInTheDPTable;


		//This variant assumes that DPTable[sumRanges[i].from] == DPTable[sumRanges[i].to] == i
		//Pointers in between are considered stale.
		//Updating DP Table entries requires computing all sums available with node u.

		void updateDPTableWithSumRanges() {
			for (const Node u : nodesNotInTheDPTable) {

				if (newSumAvailable)
					nextSumRanges = sumRanges;        //find way to avoid copies. must be possible. this is so horrible.
				AssertMsg(nextSumRanges == sumRanges, "nextSumRanges not up to date");
				newSumAvailable = false;

				const NodeWeight wu = hg.nodeWeight(u);
				AssertMsg(wu > 0, "Node has zero weight");

				for (const SummableRange &sr : sumRanges) {
					for (NodeWeight new_sum = sr.from + wu, _end = std::min(sr.to + wu, maxSubsetSumWeight);
						 new_sum <= _end; ++new_sum) {
						if (!isSummable(new_sum)) {
							newSumAvailable = true;
							DPTable[new_sum].node = u;

							NodeWeight left(new_sum - 1), right(new_sum + 1);
							Index leftIndex = DPTable[left].sumsIndex, rightIndex = DPTable[right].sumsIndex;
							bool hasLeft = DPTable[left].summable();
							bool hasRight = DPTable[right].summable();    //new_sum + 1 is valid because of right-ward sentinel.

							if (hasLeft &&
								hasRight) { //merge ranges. keep left range, and extend it to cover the right range
								AssertMsg(nextSumRanges[DPTable[left].sumsIndex].to == left,
										  "hasLeft && hasRight: left range does not extend to new_sum-1");
								AssertMsg(nextSumRanges[DPTable[right].sumsIndex].from == right,
										  "hasLeft && hasRight: right range does not start at new_sum+1");
								DPTable[new_sum].sumsIndex = leftIndex; //bridging cell. the index is stale, but we're setting it anyway to mark it as summable

								SummableRange &leftRange = nextSumRanges[leftIndex], rightRange = nextSumRanges[rightIndex];

								//extend leftRange to cover rightRange
								DPTable[rightRange.to].sumsIndex = leftIndex;
								leftRange.to = rightRange.to;

								//delete rightRange
								SummableRange back = nextSumRanges.back();        //make copy! thus immediate pop_back is safe
								nextSumRanges.pop_back();
								nextSumRanges[rightIndex] = back;

								//update sumsIndex of the range that was swapped to rightIndex
								DPTable[back.from].sumsIndex = rightIndex;
								DPTable[back.to].sumsIndex = rightIndex;
							} else if (hasLeft) {
								//extend left range's .to by +1
								AssertMsg(nextSumRanges[DPTable[left].sumsIndex].to == left,
										  "hasLeft: left range does not extend to new_sum-1");
								nextSumRanges[DPTable[left].sumsIndex].to = new_sum;
								DPTable[new_sum].sumsIndex = leftIndex;
							} else if (hasRight) {
								//extend right range's .from by -1
								AssertMsg(nextSumRanges[DPTable[right].sumsIndex].from == right,
										  "hasRight: right range does not start at new_sum+1");
								nextSumRanges[DPTable[right].sumsIndex].from = new_sum;
								DPTable[new_sum].sumsIndex = rightIndex;
							} else {
								//start new range
								DPTable[new_sum].sumsIndex = static_cast<Index>(nextSumRanges.size());
								nextSumRanges.emplace_back(new_sum, new_sum);
							}
						}
					}
				}

				std::swap(nextSumRanges, sumRanges);
			}
		}

		//This variant assumes that DPTable[sumRanges[i].from] == i and DPTable[sumRanges[i].to] == sumRanges[i].from
		//Pointers in between point to -- potentially former -- left ends (sr.from) of the range sr. Using path compression we can keep "hop distance" to sr.from "short"
		//When a newly computed sum equals an already computed sum, we can jump to the left end of the range, get the corresponding index into sumRanges,
		//thus the entire range, and prune sums that have been previously computed
		void updateDPTableWithSumRangesAndRangePruning() {
			//Implement me if the other approach is too slow.
		}

	public:
		explicit IsolatedNodes(FlowHypergraph& hg, NodeWeight maxBlockWeight) :
				hg(hg),
				mixedIncidentHyperedges(hg.numNodes(), HyperedgeIndex(0)),
				maxSubsetSumWeight(maxBlockWeight),
				DPTable(maxBlockWeight + 2, TableEntry())
				//We could use maxBlocKWeight+2 for a right-ward sentinel. With +3 even a left-ward sentinel, but that would required reindexing everything by +1 which is nasty
				//So, only rightward sentinel
		{
			sumRanges.emplace_back(NodeWeight(0), NodeWeight(0));
			DPTable[0].sumsIndex = 0;
		}

		const std::vector<SummableRange>& getSumRanges() const {
			return sumRanges;
		}

		bool isSummable(const NodeWeight w) const {
			Assert(w < DPTable.size());
			return DPTable[w].summable();
		}

		void add(const Node u) {
			nodes.push_back(u);
			nodesNotInTheDPTable.push_back(u);
			weight += hg.nodeWeight(u);
		}

		bool isDPTableUpToDate() const {
			return nodesNotInTheDPTable.empty();
		}

		void updateDPTable() {
			if (weight + 1 > DPTable.size()) {
				//at the moment we're allocating for @maxBlocKWeight elements in DPTable, so this branch is never executed. in the future we might care about saving some memory
				//if memory actually becomes critical, we can also compress the DPTable by a two-level indexing approach, where ranges only consume one entry. compression called infrequently
				DPTable.resize(std::min((size_t)2*(weight+1), (size_t)maxSubsetSumWeight + 1), TableEntry());
			}
			updateDPTableWithSumRanges();
			//Internal_UpdateDPTableWithSumRangesAndRangePruning();		NOT IMPLEMENTED YET. Note sure if faster or necessary
			nodesNotInTheDPTable.clear();
		}

		std::vector<Node> extractSubset(NodeWeight sum) const {
			AssertMsg(isSummable(sum),
					  std::string("Trying to extract subset for not achieved subset sum. ")
					  + (isDPTableUpToDate() ? "Call updateDPTable() before calling this method. There are nodes that are not included in the table yet." : " ")
			);

			std::vector<Node> result;
			while (sum > 0) {
				const Node u = DPTable[sum].node;
				result.push_back(u);
				sum -= hg.nodeWeight(u);
			}
			return result;
		}

		std::pair<std::vector<Node>, std::vector<Node>> extractBipartition(NodeWeight sum) const {
			auto first = extractSubset(sum);
			auto second = setDifference(nodes, first, hg.numNodes());
			return std::make_pair(std::move(first), std::move(second));
		}


		bool isCandidate(const Node u) const {
			return mixedIncidentHyperedges[u] == hg.degree(u);
		}

	};
}