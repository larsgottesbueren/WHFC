#include <iostream>
#include "io/hmetis_io.h"
#include "io/whfc_io.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/task_scheduler_init.h>

#include "algorithm/parallel_push_relabel.h"
#include "algorithm/sequential_push_relabel.h"

#include "algorithm/hyperflowcutter.h"
#include "algorithm/dinic.h"
// #include "algorithm/graph_push_relabel.h"

namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		static constexpr bool log = false;
		TimeReporter tr;
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		LOGGER << "(s,t,max f) =" << s << t << info.upperFlowBound;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		LOGGER << "(n,m,p) =" << hg.numNodes() << hg.numHyperedges() << hg.numPins();
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");



		int seed = 42;
		HyperFlowCutter<Dinic> hfc(hg, seed);
		for (int i = 0; i < 2; ++i) hfc.cs.setMaxBlockWeight(i, std::numeric_limits<NodeWeight>::max());
		hfc.cs.initialize(s, t);
		tr.start("dinitz");
		hfc.flow_algo.exhaustFlow(hfc.cs);
		tr.stop("dinitz");
		hfc.cs.flipViewDirection();
		hfc.flow_algo.growReachable(hfc.cs);
		hfc.cs.flipViewDirection();
		LOGGER <<"Dinic. f =" <<  hfc.cs.flowValue;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);


		tr.start("seq push relabel");
		SequentialPushRelabel spr(hg);
		Flow f_spr = spr.computeFlow(s, t);
		tr.stop("seq push relabel");
		LOGGER << "Seq Push-Relabel f =" << f_spr;

/*
		tr.start("push relabel");
		ParallelPushRelabel pr(hg);
		// pr.dinitz_flow_value = hfc.cs.flowValue;
		Flow f = pr.computeFlow(s, t);
		tr.stop("push relabel");
		LOGGER << "Push-Relabel f =" << f;
*/
/*
		tr.start("graph push relabel");
		GraphPushRelabel gpr(hg, false);
		Flow f_gpr = gpr.computeFlow(s, t);
		tr.stop("graph push relabel");
		LOGGER << "Graph Push-Relabel f =" << f_gpr;
*/
		if (hfc.cs.flowValue != f_spr) {
			std::cout << filename << " flow sequential push relabel = " << f_spr << " flow dinitz = " << hfc.cs.flowValue << std::endl;
		}
		// tr.report(std::cout);
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
