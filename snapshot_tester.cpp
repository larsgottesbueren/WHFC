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

		int seed = 0;
		using FlowAlgorithm = ParallelPushRelabel;
		HyperFlowCutter<FlowAlgorithm> hfc(hg, seed);
		hfc.setFlowBound(info.upperFlowBound);
		for (int i = 0; i < 2; ++i)
			hfc.cs.setMaxBlockWeight(i, info.maxBlockWeight[i]);

		WHFC_IO::readRandomGeneratorState(filename, hfc.cs.rng);

		bool time_limit_exceeded = false;
		size_t measure_step = 0;
		auto start_time = std::chrono::high_resolution_clock::now();
		int time_limit = 3600; // seconds
		size_t num_cuts = 0;
		Flow last_cut = 0;

		auto on_cut = [&] {
			if (hfc.cs.flow_algo.flow_value != last_cut) {
				last_cut = hfc.cs.flow_algo.flow_value;
				num_cuts++;
			}
			if (++measure_step == 50) {
				measure_step = 0;
				auto now = std::chrono::high_resolution_clock::now();
				if (std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count() > time_limit) {
					time_limit_exceeded = true;
					return false;
				}
			}
			return true;
		};

		hfc.timer.start();
		bool result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t, on_cut);
		hfc.timer.stop();
		std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
		/*
		 * header
		 * graph,algorithm,seed,threads,improved,flow,flowbound,time,mbc_time,time_limit_exceeded,num_cuts
		 */
		std::cout << base_filename << ",FlowCutter-ParPR," << seed << "," << threads << ",";
		std::cout << (result ? "yes" : "no") << ",";
		std::cout << hfc.cs.flow_algo.flow_value << "," << info.upperFlowBound << ",";
		std::cout << hfc.timer.get("HyperFlowCutter").count() << "," << hfc.timer.get("MBMC").count() << ",";
		std::cout << (time_limit_exceeded ? "yes" : "no") << ",";
		std::cout << num_cuts;

		std::cout << std::endl;


		std::cout << V(result) << " " << V(hfc.cs.flow_algo.flow_value) << std::endl;
		hfc.timer.report(std::cout);
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
