#include <iostream>
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "algorithm/parallel_push_relabel.h"


namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		std::cout << s << " " << t << info.upperFlowBound << std::endl;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		ParallelPushRelabel pr(hg);
		Flow f = pr.computeFlow(s, t);
		std::cout << "f = " << f << std::endl;
	}
}

int main(int argc, const char* argv[]) {
	if (argc < 2 || argc > 3)
		throw std::runtime_error("Usage: ./FlowTester hypergraphfile");
	std::string hgfile = argv[1];
	whfc::runSnapshotTester(hgfile);
	return 0;
}
