#include <iostream>
#include "io/hmetis_io.h"
#include "io/whfc_io.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/task_scheduler_init.h>

#include "algorithm/parallel_push_relabel.h"
#include "algorithm/sequential_push_relabel.h"
#include "algorithm/graph_push_relabel.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/dinic.h"

namespace whfc {
	void runSnapshotTester(const std::string& filename) {
		static constexpr bool log = false;
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		LOGGER << "(s,t,max f) =" << s << t << info.upperFlowBound;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		LOGGER << "(n,m,p) =" << hg.numNodes() << hg.numHyperedges() << hg.numPins();
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");


		int seed = 42;
	 	TimeReporter tr;
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


		/*
				tr.start("seq push relabel");
				SequentialPushRelabel spr(hg);
				Flow f_spr = spr.computeFlow(s, t);
				tr.stop("seq push relabel");
				LOGGER << "Seq Push-Relabel f =" << f_spr;
		*/
/*
		std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);

		int max_num_threads = 128;
		for (int i = 0; i < 5; ++i) {
			for (int threads = 1; threads <= max_num_threads; threads *= 2) {
				tbb::task_scheduler_init tsi(threads);
				whfc::pinning_observer thread_pinner;
				thread_pinner.observe(true);
				ParallelPushRelabel pr(hg);
				pr.computeFlow(s, t);
				std::cout << base_filename << "," << "ParPR-RL" << "," << threads << "," << pr.timer.get("push relabel").count() << std::endl;
			}
		}
*/
		tr.start("graph push relabel");
		GraphPushRelabel gpr(hg, true);
		Flow f_gpr = gpr.computeFlow(s, t);
		tr.stop("graph push relabel");
		LOGGER << "Graph Push-Relabel f =" << f_gpr;

		if (f_gpr != hfc.cs.flowValue)
			std::cout <<  "\033[1;31mflow value mismatch\033[0m\n" << filename << " flow push relabel = " << f_gpr << " flow dinitz = " << hfc.cs.flowValue << std::endl;
		// std::cout << "time dinitz " << tr.get("dinitz").count() << " s" << std::endl;
// 		tr.report(std::cout);
	}


}

int main(int argc, const char* argv[]) {
	if (argc > 3 || argc < 2)
	 	throw std::runtime_error("Usage: ./FlowTester hypergraphfile #threads");
	std::string hgfile = argv[1];
	int threads = 1;
	if (argc == 3)
		threads = std::stoi(argv[2]);
	/*
	tbb::task_scheduler_init tsi(threads);
	whfc::pinning_observer thread_pinner;
	if (argc == 3)
		thread_pinner.observe(true);
	*/
	whfc::runSnapshotTester(hgfile);
	return 0;
}
