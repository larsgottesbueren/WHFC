#include <iostream>
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "datastructure/flow_hypergraph_builder.h"
#include "algorithm/parallel_push_relabel.h"
#include "algorithm/sequential_push_relabel.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/global_control.h>

namespace whfc {
	void pin() {
		int slot = 0;
		cpu_set_t target_mask;
		CPU_ZERO(&target_mask);
		CPU_SET(slot, &target_mask);
		size_t size = CPU_ALLOC_SIZE(std::thread::hardware_concurrency());
		int err = sched_setaffinity(0, size, &target_mask);
		if (err) { std::cout << "couldnt set thread affinity" << std::endl; std::exit(1); }
	}

	void unpin() {
		cpu_set_t target_mask;
		CPU_ZERO(&target_mask);
		for (size_t core = 0; core < std::thread::hardware_concurrency(); ++core) {
			CPU_SET(core, &target_mask);
		}
		size_t size = CPU_ALLOC_SIZE(std::thread::hardware_concurrency());
		int err = sched_setaffinity(0, size, &target_mask);
		if (err) { std::cout << "couldnt reset thread affinity" << std::endl; std::exit(1); }
	}


	void runSnapshotTester(const std::string& filename, int max_threads) {
		pin();
		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		// std::cout << s << " " << t << " " << info.maxBlockWeight[0] << " " << info.maxBlockWeight[1] << " " << info.upperFlowBound<< std::endl;
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		unpin();

		for (int threads = 32; threads <= 32; threads *= 2) {
			auto gc = tbb::global_control{tbb::global_control::max_allowed_parallelism, threads};
			whfc::pinning_observer thread_pinner;
			thread_pinner.observe(true);

			for (int rep = 0; rep < 1; ++rep) {

				int seed = 0;
				using FlowAlgorithm = ParallelPushRelabel;
				// using FlowAlgorithm = SequentialPushRelabel;
				HyperFlowCutter<FlowAlgorithm> hfc(hg, seed);
				hfc.setFlowBound(info.upperFlowBound);
				hfc.forceSequential(false);
				hfc.setBulkPiercing(true);
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

				std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);

				hfc.timer.start();
				bool result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t, on_cut);
				hfc.timer.stop();
				/*
				 * header
				 * graph,algorithm,seed,threads,improved,flow,flowbound,time,mbc_time,time_limit_exceeded,num_cuts,discharge,global relabel,update,source cut,saturate,assimilate,pierce
				 */
				std::cout << base_filename << ",FlowCutter,";
				// std::cout << seed << ",";
				std::cout << rep << ",";
				std::cout << threads << ",";
				std::cout << (result ? "yes" : "no") << ",";
				std::cout << hfc.cs.flow_algo.flow_value << "," << info.upperFlowBound << ",";
				std::cout << hfc.timer.get("HyperFlowCutter").count() << "," << hfc.timer.get("MBMC").count() << ",";
				std::cout << (time_limit_exceeded ? "yes" : "no") << ",";
				std::cout << num_cuts;

				auto& f = hfc.cs.flow_algo;
				std::cout << "," <<  f.discharge_time << "," << f.global_relabel_time << "," << f.update_time << "," << f.source_cut_time << "," << f.saturate_time;
				std::cout << "," << hfc.assimilate_time << "," << hfc.pierce_time;

				std::cout << std::endl;
			}

			/*
			std::cout << V(result) << " " << V(hfc.cs.flow_algo.flow_value) << std::endl;
			hfc.timer.report(std::cout);
			hfc.timer.clear();
			 */
		}


	}
}

int main(int argc, const char* argv[]) {
	if (argc < 2 || argc > 4)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile [#threads optional]");
	std::string hgfile = argv[1];
	int threads = 1;
	if (argc == 3) {
		threads = std::stoi(argv[2]);
	}
	whfc::runSnapshotTester(hgfile, threads);
	return 0;
}
