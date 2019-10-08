#include <iostream>
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/ford_fulkerson.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "util/random.h"
#include "datastructure/flow_hypergraph_builder.h"


namespace whfc {
	void runSnapshotTester(const std::string& filename, std::string& interleaving) {
		
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		NodeWeight mbw = info.maxBlockWeight;
		std::cout << s << " " << t << " " << mbw << " " << info.upperFlowBound<< std::endl;
		
		FlowHypergraph hg = HMetisIO::readFlowHypergraph(filename);
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");
		
		HyperFlowCutter<ScalingFordFulkerson> hfc(hg, mbw);
		hfc.upperFlowBound = info.upperFlowBound;
		
		auto time = time_now();
		
		if (interleaving == "flowbased")
			hfc.runUntilBalancedOrFlowBoundExceeded(s, t);
		else if (interleaving == "cutbased")
			hfc.findCutsUntilBalancedOrFlowBoundExceeded(s, t);
		else
			throw std::runtime_error("Unknown interleaving option");
		
		std::cout << second_duration(time_now() - time).count() << " [s]" << std::endl;
	}
}

int main(int argc, const char* argv[]) {
	whfc::Random::setSeed(42);
	if (argc != 3)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile interleaving-style (flowbased or cutbased)");
	std::string hgfile = argv[1];
	std::string interleavingstyle = argv[2];
	whfc::runSnapshotTester(hgfile, interleavingstyle);
	return 0;
}