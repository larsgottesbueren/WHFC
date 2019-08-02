#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"

//#include <tlx/container.hpp>
//one optimization for speeding up the balance constraint check might require sorted sums

namespace whfc {
	class IsolatedNodes {
	private:
		FlowHypergraph& hg;
	public:
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<HyperedgeIndex> mixedIncidentHyperedges;

		//If the subset sum approach gets too slow, there are two alternatives:
		//	1) try scaling the weights of all nodes in the flow network --> more locality in the DP table
		//	2) just assign the nodes ad-hoc, after growAssimilated. shouldn't be too bad. works well enough with AAP
	private:
		NodeWeight maxSubsetSumWeight = NodeWeight(0);
		std::vector<Node> DPTable;
		std::vector<NodeWeight> sums;

		std::vector<Node> nodesNotInTheDPTable;

	public:
		explicit IsolatedNodes(FlowHypergraph& hg) :
				hg(hg),
				mixedIncidentHyperedges(hg.numNodes(), HyperedgeIndex(0)),
				maxSubsetSumWeight(hg.totalNodeWeight()),
				DPTable(200, invalidNode)
		{
			sums.emplace_back(0);
		}

		bool isSummable(const NodeWeight w) const {
			assert(w < DPTable.size());
			return DPTable[w] != invalidNode;
		}

		void add(const Node u) {
			nodes.push_back(u);
			nodesNotInTheDPTable.push_back(u);
			weight += hg.nodeWeight(u);
		}

		bool isDPTableUpToDate() const {
			return nodesNotInTheDPTable.empty();
		}

		//TODO Represent entire ranges -- that can be summed -- more efficiently. (not less memory, I want to traverse them faster)

		void updateDPTable() {
			if (weight + 1 > DPTable.size())
				DPTable.resize(std::min((size_t)2*(weight+1), (size_t)maxSubsetSumWeight + 1), invalidNode);

			for (const Node u : nodesNotInTheDPTable) {
				const NodeWeight wu = hg.nodeWeight(u);
				for (size_t i = 0, it_end = sums.size(); i < it_end; ++i) {
					const NodeWeight now = sums[i] + wu;
					if (!isSummable(now)) {
						DPTable[now] = u;
						sums.push_back(now);
					}
				}
			}
			nodesNotInTheDPTable.clear();
		}

		std::vector<Node> extractSubset(NodeWeight sum) {
			AssertMsg(sum == 0 || isSummable(sum), "Trying to extract subset for not achieved subset sum");
			AssertMsg(isDPTableUpToDate(), "There are still nodes that are not included in the DP table. Make sure to call updateDPTable() before calling this method.");
			std::vector<Node> result;
			while (sum > 0) {
				const Node u = DPTable[sum];
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