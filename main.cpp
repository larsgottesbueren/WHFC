#include <iostream>
#include "datastructure/hypergraph.h"
#include "util/range.h"
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/grow_assimilated.h"
#include "algorithm/ford_fulkerson.h"


namespace whfc {
	void run() {
		NodeWeight mbw = NodeWeight(42);

		FlowHypergraph hg;
		std::cout << "init hg" << std::endl;
		std::cout << hg.numNodes() << std::endl;

		ScalingFordFulkerson fs(hg);
		std::cout << "init fs" << std::endl;
		std::cout << fs.nodes_to_scan.size() << std::endl;

		CutterState<ScalingFordFulkerson> cs(hg, mbw);
		std::cout << "init cs" << std::endl;
		std::cout << cs.sourcePiercingNodes.size() << std::endl;

		HyperFlowCutter<ScalingFordFulkerson> hfc(hg, mbw);
		std::cout << "init hfc" << std::endl;
		std::cout << hfc.upperFlowBound << std::endl;
	}
}

int main() {
	whfc::run();
	return 0;
}