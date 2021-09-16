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
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		std::cout << "(s,t,max f) = " << s << " " << t << " " << info.upperFlowBound << std::endl;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		std::cout << "(n,m,p) = " << hg.numNodes() << " " << hg.numHyperedges() << " " << hg.numPins() << std::endl;

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		ParallelPushRelabel pr(hg);
		Flow f = pr.computeFlow(s, t);
		std::cout << "Push-Relabel f = " << f << std::endl;



		int seed = 42;
		HyperFlowCutter<Dinic> hfc(hg, seed);
		for (int i = 0; i < 2; ++i) hfc.cs.setMaxBlockWeight(i, info.maxBlockWeight[i]);
		hfc.cs.initialize(s, t);
		hfc.flow_algo.exhaustFlow(hfc.cs);
		std::cout <<"Dinic. f = " <<  hfc.cs.flowValue << std::endl;
	}

	void runPlainHypergraph(const std::string& filename, int is, int it) {
		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		std::cout << "(n,m,p) = " << hg.numNodes() << " " << hg.numHyperedges() << " " << hg.numPins() << std::endl;

		Node s(is); Node t(it);
		std::cout << "(s,t) = " << s << " " << t << std::endl;

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		ParallelPushRelabel pr(hg);
		Flow f = pr.computeFlow(s, t);
		std::cout << "Push-Relabel f = " << f << std::endl;


		int seed = 42;
		HyperFlowCutter<Dinic> hfc(hg, seed);
		for (int i = 0; i < 2; ++i) hfc.cs.setMaxBlockWeight(i, std::numeric_limits<NodeWeight>::max());
		hfc.cs.initialize(s, t);
		hfc.flow_algo.exhaustFlow(hfc.cs);
		std::cout <<"Dinic. f = " <<  hfc.cs.flowValue << std::endl;
	}
}

int main(int argc, const char* argv[]) {
	//if (argc != 3)
	//	throw std::runtime_error("Usage: ./FlowTester hypergraphfile #threads");
	std::string hgfile = argv[1];
	int threads = std::stoi(argv[2]);
	tbb::task_scheduler_init tsi(threads);
	whfc::pinning_observer thread_pinner;
	thread_pinner.observe(true);

	// whfc::runSnapshotTester(hgfile);
	whfc::runPlainHypergraph(hgfile, std::stoi(argv[3]), std::stoi(argv[4]));
	return 0;
}
