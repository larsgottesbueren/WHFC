#include <iostream>
#include "io/hmetis_io.h"
#include "io/whfc_io.h"

#include "algorithm/parallel_push_relabel.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/task_scheduler_init.h>

#include "algorithm/hyperflowcutter.h"
#include "algorithm/dinic.h"

namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		static constexpr bool log = true;
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		LOGGER << "(s,t,max f) =" << s << t << info.upperFlowBound;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		LOGGER << "(n,m,p) =" << hg.numNodes() << hg.numHyperedges() << hg.numPins();


		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		TimeReporter tr;
		tr.start("push relabel");
		ParallelPushRelabel pr(hg);
		Flow f = pr.computeFlow(s, t);
		tr.stop("push relabel");
		LOGGER << "Push-Relabel f =" << f;

		int seed = 42;
		tr.start("dinitz");
		HyperFlowCutter<Dinic> hfc(hg, seed);
		for (int i = 0; i < 2; ++i) hfc.cs.setMaxBlockWeight(i, info.maxBlockWeight[i]);
		hfc.cs.initialize(s, t);
		hfc.flow_algo.exhaustFlow(hfc.cs);
		tr.stop("dinitz");
		hfc.cs.flipViewDirection();
		hfc.flow_algo.growReachable(hfc.cs);
		hfc.cs.flipViewDirection();
		LOGGER <<"Dinic. f =" <<  hfc.cs.flowValue;
		tr.report(std::cout);
	}


}

int main(int argc, const char* argv[]) {
	if (argc != 3)
		throw std::runtime_error("Usage: ./FlowTester hypergraphfile #threads");
	std::string hgfile = argv[1];
	int threads = std::stoi(argv[2]);
	tbb::task_scheduler_init tsi(threads);
	whfc::pinning_observer thread_pinner;
	thread_pinner.observe(true);

	whfc::runSnapshotTester(hgfile);
	return 0;
}
