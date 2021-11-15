#include <iostream>
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "datastructure/flow_hypergraph_builder.h"
#include "algorithm/parallel_push_relabel.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/task_scheduler_init.h>

namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		int threads = 1;
		tbb::task_scheduler_init tsi(threads);
		whfc::pinning_observer thread_pinner;
		thread_pinner.observe(true);

		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		// std::cout << s << " " << t << " " << info.maxBlockWeight[0] << " " << info.maxBlockWeight[1] << " " << info.upperFlowBound<< std::endl;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		int seed = 42;
		using FlowAlgorithm = ParallelPushRelabel;
		HyperFlowCutter<FlowAlgorithm> hfc(hg, seed);
		hfc.setFlowBound(info.upperFlowBound);
		for (int i = 0; i < 2; ++i)
			hfc.cs.setMaxBlockWeight(i, info.maxBlockWeight[i]);

		WHFC_IO::readRandomGeneratorState(filename, hfc.cs.rng);

		hfc.timer.start();
		bool result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t);
		hfc.timer.stop();
		// std::cout << V(result) << std::endl;
		// hfc.timer.report(std::cout);
		hfc.timer.clear();
	}
}

int main(int argc, const char* argv[]) {
	if (argc < 2 || argc > 3)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile");
	std::string hgfile = argv[1];
	whfc::runSnapshotTester(hgfile);
	return 0;
}
