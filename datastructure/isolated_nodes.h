#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"
#include "../util/gcd.h"

//#include <tlx/container.hpp>
//one optimization for speeding up the balance constraint check might require sorted sumRanges

namespace whfc {
	class IsolatedNodes {
	private:
		FlowHypergraph& hg;
	public:
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<HyperedgeIndex> mixedIncidentHyperedges;


	private:

		//If the subset sum approach gets too slow, there are two alternatives:
		//	1) try scaling the weights of all nodes in the flow network --> more locality in the DP table
		//	2) just assign the nodes ad-hoc, after growAssimilated. shouldn't be too bad. works well enough with AAP
		NodeWeight maxSubsetSumWeight = NodeWeight(0);
		NodeWeight weightScaling = NodeWeight(1);


		struct TableEntry {
			Node node;
			Index sumsIndex;
			TableEntry() : node(invalidNode), sumsIndex(invalidIndex) { }
			bool summable() const {
#ifndef NDEBUG
				if (sumsIndex != invalidIndex)						//sumsIndex = 0 is reserved for the first range.
					AssertMsg(node != invalidNode || sumsIndex == 0, "Summable but no node");
#endif
				return sumsIndex != invalidIndex;
			}
		};

		std::vector<TableEntry> DPTable;

		struct SummableRange {
			NodeWeight from, to;
			SummableRange(NodeWeight _from, NodeWeight _to) : from(_from), to(_to) { }
			bool inRange(const NodeWeight w) const { return from <= w && w <= to; }
		};
		std::vector<SummableRange> sumRanges;
		std::vector<SummableRange> nextSumRanges;
		bool newSumAvailable = true;

		std::vector<Node> nodesNotInTheDPTable;


		//This variant assumes that DPTable[sumRanges[i].from] == DPTable[sumRanges[i].to] == i
		//Pointers in between are considered stale.
		//Updating DP Table entries requires computing all sums available with node u.
		void Internal_UpdateDPTableWithSumRanges() {
			for (const Node u : nodesNotInTheDPTable) {

				if (newSumAvailable)
					nextSumRanges = sumRanges;
				newSumAvailable = false;

				const NodeWeight wu = hg.nodeWeight(u);
				AssertMsg(wu > 0, "Node has zero weight");

				for (const SummableRange& sr : sumRanges) {
					for (NodeWeight new_sum = sr.from + wu; new_sum <= std::min(sr.to + wu, maxSubsetSumWeight); ++new_sum) {
						if (!isSummable(new_sum)) {
							newSumAvailable = true;
							DPTable[new_sum].node = u;

							NodeWeight left(new_sum - 1), right(new_sum + 1);
							Index leftIndex = DPTable[left].sumsIndex, rightIndex = DPTable[right].sumsIndex;
							bool hasLeft = DPTable[left].summable();
							bool hasRight = DPTable[right].summable();	//new_sum + 1 is valid because of right-ward sentinel.
							//If no sentinels. bool hasRight = new_sum + 1 <= maxSubsetSumWeight && DPTable[new_sum+1].summable();

							if (hasLeft && hasRight) { //merge ranges. keep left range, and extend it to cover the right range
								AssertMsg(nextSumRanges[DPTable[left].sumsIndex].to == left, "hasLeft && hasRight: left range does not extend to new_sum-1");
								AssertMsg(nextSumRanges[DPTable[right].sumsIndex].from == right, "hasLeft && hasRight: right range does not start at new_sum+1");
								DPTable[new_sum].sumsIndex = leftIndex; //bridging cell. the index is stale, but we're setting it anyway to mark it as summable

								SummableRange& leftRange = nextSumRanges[leftIndex], rightRange = nextSumRanges[rightIndex];

								//extend leftRange to cover rightRange
								DPTable[rightRange.to].sumsIndex = leftIndex;
								leftRange.to = rightRange.to;

								//delete rightRange
								SummableRange back = nextSumRanges.back();		//copy! immediate pop_back is safe
								nextSumRanges.pop_back();
								nextSumRanges[rightIndex] = back;

								//update sumsIndex of the range that was swapped to rightIndex
								DPTable[back.from].sumsIndex = rightIndex;
								DPTable[back.to].sumsIndex = rightIndex;
							}
							else if (hasLeft) {
								//extend left range's .to by +1
								AssertMsg(nextSumRanges[DPTable[left].sumsIndex].to == left, "hasLeft: left range does not extend to new_sum-1");
								nextSumRanges[DPTable[left].sumsIndex].to = new_sum;
								DPTable[new_sum].sumsIndex = leftIndex;
							}
							else if (hasRight) {
								//extend right range's .from by -1
								AssertMsg(nextSumRanges[DPTable[right].sumsIndex].from == right, "hasRight: right range does not start at new_sum+1");
								nextSumRanges[DPTable[right].sumsIndex].from = new_sum;
								DPTable[new_sum].sumsIndex = rightIndex;
							}
							else {
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
		void Internal_UpdateDPTableWithSumRangesAndRangePruning() {
			//Implement me.
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

			std::vector<NodeWeight> weights;
			for (Node u : hg.nodeIDs())
				weights.push_back(hg.nodeWeight(u));
			weightScaling = GreatestCommonDivisor::compute(weights);
			//TODO scale everything exposed to the outside with this factor
			//also maxSubsetSumWeight etc.
		}

		bool isSummable(const NodeWeight w) const {
			assert(w < DPTable.size());
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
				//if memory actually becomes critical, we can also compress the DPTable by a two-level indexing approach, where ranges only consume one entry
				DPTable.resize(std::min((size_t)2*(weight+1), (size_t)maxSubsetSumWeight + 1), TableEntry());
			}
			Internal_UpdateDPTableWithSumRanges();
			//Internal_UpdateDPTableWithSumRangesAndRangePruning();		NOT IMPLEMENTED YET. Note sure if faster or necessary
			nodesNotInTheDPTable.clear();
		}

		std::vector<Node> extractSubset(NodeWeight sum) {
			AssertMsg(sum == 0 || isSummable(sum), "Trying to extract subset for not achieved subset sum");
			AssertMsg(isDPTableUpToDate(), "There are still nodes that are not included in the DP table. Make sure to call updateDPTable() before calling this method.");
			std::vector<Node> result;
			while (sum > 0) {
				const Node u = DPTable[sum].node;
				result.push_back(u);
				sum -= hg.nodeWeight(u);
			}
			return result;
		}

		template<typename PredicateIsIsolated>
		void accommodateNewlyMixedHyperedge(const Hyperedge e, PredicateIsIsolated isIsolated) {
			for (const auto& px : hg.pinsOf(e)) {
				const Node p = px.pin;
				mixedIncidentHyperedges[p]++;
				/*
				 * Previously, we identified candidates, via mixedIncidentHyperedges[p] == hg.degree(p), and later accepted them as isolated, only if they were not settled.
				 * This was particularly necessary for piercing hyperedges.
				 * However, it is suboptimal since the later settled candidates could just as well be moved between the blocks without modifying the cut.
				 * Immediately isolating them can yield better final cuts and is substantially less code.
				 * TODO this does not seem quite done yet. think about it and then come back.
				 * However, we have to move the ones that were marked reachable out of that set. Also make sure that it is IMPOSSIBLE to actually ever reach an isolated node.
				 */
				if (isIsolated(p))
					add(p);
			}
		}

		bool isCandidate(const Node u) const {
			return mixedIncidentHyperedges[u] == hg.degree(u);
		}

		void settleNode(const Node u) {
			maxSubsetSumWeight -= hg.nodeWeight(u);
		}

	};
}