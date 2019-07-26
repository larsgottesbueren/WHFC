#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"

namespace whfc {
	class IsolatedNodes {
	private:
		const FlowHypergraph& hg;
	public:
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<HyperedgeIndex> mixedIncidentHyperedges;

		explicit IsolatedNodes(const FlowHypergraph& hg) : hg(hg), mixedIncidentHyperedges(hg.numNodes(), HyperedgeIndex(0)) { }

		void add(const Node u) {
			nodes.push_back(u);
			weight += hg.nodeWeight(u);
		}

		bool isCandidate(const Node u) const {
			return mixedIncidentHyperedges[u] == hg.degree(u);
		}

	};
}