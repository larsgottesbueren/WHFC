#include <iostream>
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/ford_fulkerson.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "util/random.h"
#include "datastructure/flow_hypergraph_builder.h"


namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		NodeWeight mbw = info.maxBlockWeight;
		std::cout << s << " " << t << " " << mbw << " " << info.upperFlowBound<< std::endl;
		
		FlowHypergraph hg = HMetisIO::readFlowHypergraph(filename);
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");
		
		HyperFlowCutter<BasicEdmondsKarp> hfc(hg, mbw);
		hfc.upperFlowBound = info.upperFlowBound;
		
		auto time = time_now();
		hfc.initialize(s,t);
		hfc.runUntilBalanced();
		std::cout << second_duration(time_now() - time).count() << " [s]" << std::endl;
	}
}

int main(int argc, const char* argv[]) {
	whfc::Random::setSeed(42);
	if (argc != 2)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile");
	std::string hgfile = argv[1];
	whfc::runSnapshotTester(hgfile);
	return 0;
}