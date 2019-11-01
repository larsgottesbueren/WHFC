#include <iostream>
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/ford_fulkerson.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "util/random.h"
#include "datastructure/flow_hypergraph_builder.h"
#include "algorithm/dinic.h"


namespace whfc {
	void runSnapshotTester(const std::string& filename, std::string& interleaving) {
		
		using FlowAlgorithm = Dinic;
		//using FlowAlgorithm = ScalingDinic;
		
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		NodeWeight mbw = info.maxBlockWeight;
		std::cout << s << " " << t << " " << mbw << " " << info.upperFlowBound<< std::endl;
		
		FlowHypergraphBuilder hg = HMetisIO::readFlowHypergraphWithBuilder(filename);
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");
		
		int seed = 42;
		HyperFlowCutter<FlowAlgorithm> hfc(hg, mbw, seed);
		hfc.upperFlowBound = info.upperFlowBound;
		
		WHFC_IO::readRandomGeneratorState(filename);
		
		hfc.timer.start();
		if (interleaving == "flowbased")
			hfc.runUntilBalancedOrFlowBoundExceeded(s, t);
		else if (interleaving == "cutbased")
			hfc.findCutsUntilBalancedOrFlowBoundExceeded(s, t);
		else
			throw std::runtime_error("Unknown interleaving option");
		hfc.timer.stop();
		hfc.timer.report(std::cout);
		hfc.timer.clear();
		
		std::cout << "Reread" << std::endl;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		
		std::cout << "Reset" << std::endl;
		hfc.reset();
		std::cout << "Run again" << std::endl;
		
		hfc.timer.start();
		if (interleaving == "flowbased")
			hfc.runUntilBalancedOrFlowBoundExceeded(s, t);
		else if (interleaving == "cutbased")
			hfc.findCutsUntilBalancedOrFlowBoundExceeded(s, t);
		else
			throw std::runtime_error("Unknown interleaving option");
		hfc.timer.stop();
		hfc.timer.report(std::cout);
		hfc.timer.clear();
		
		
	}
}

int main(int argc, const char* argv[]) {
	if (argc < 2 || argc > 3)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile interleaving-style (flowbased or cutbased)");
	std::string hgfile = argv[1];
	std::string interleavingstyle = "flowbased";
	if (argc == 3)
		interleavingstyle = argv[2];
	whfc::runSnapshotTester(hgfile, interleavingstyle);
	return 0;
}