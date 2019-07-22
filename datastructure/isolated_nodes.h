#pragma once

#include "../definitions.h"

namespace whfc {
	class IsolatedNodes {
	public:
		explicit IsolatedNodes(size_t nNodes) : mixedIncidentHyperedges(nNodes, HyperedgeIndex(0)) { }
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<HyperedgeIndex> mixedIncidentHyperedges;
		void add(const Node u, const NodeWeight w) {
			nodes.push_back(u);
			weight += w;
		}
	};
}