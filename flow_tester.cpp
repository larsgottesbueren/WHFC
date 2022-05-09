#include <iostream>
#include "io/hmetis_io.h"
#include "io/whfc_io.h"

#include "util/tbb_thread_pinning.h"
#include <tbb/task_scheduler_init.h>

#include "algorithm/parallel_push_relabel.h"
#include "algorithm/parallel_push_relabel_block.h"
#include "algorithm/sequential_push_relabel.h"

namespace whfc {
	void runSnapshotTester(const std::string& filename, int max_num_threads) {
		static constexpr bool log = false;
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s; Node t = info.t;
		LOGGER << "(s,t,max f) =" << s << t << info.upperFlowBound;

		FlowHypergraphBuilder hg;
		HMetisIO::readFlowHypergraphWithBuilder(hg, filename);
		LOGGER << "(n,m,p) =" << hg.numNodes() << hg.numHyperedges() << hg.numPins();
		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");

		std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);

		for (int threads = 32; threads <= 32; threads *= 2) {
			tbb::task_scheduler_init tsi(threads);
			whfc::pinning_observer thread_pinner;
			thread_pinner.observe(true);
			for (int i = 0; i < 1; ++i) {
				TimeReporter timer;
				ParallelPushRelabel pr(hg);
				timer.start("ParPR-RL");
				Flow f_pr = pr.computeMaxFlow(s, t);
				timer.stop("ParPR-RL");
				//std::cout << base_filename << "," << i << ",ParPR-RL," << threads << "," << timer.get("ParPR-RL").count() << std::endl;

				/*
				 * header
				 * graph,algorithm,seed,threads,time,discharge,global relabel,update,saturate
				 */
                std::cout << base_filename << ",ParPR-RL,";
                std::cout << i << ",";
                std::cout << threads << ",";
                std::cout << timer.get("ParPR-RL").count() << ",";
                std::cout << "," <<  pr.discharge_time << "," << pr.global_relabel_time << "," << pr.update_time << "," << pr.saturate_time;
                std::cout << std::endl;

                /*
                ParallelPushRelabelBlock prb(hg);
                timer.start("ParPR-Block");
                Flow f_pr_block = pr.computeMaxFlow(s, t);
                timer.stop("ParPR-Block");
                std::cout << base_filename << "," << i << ",ParPR-Block," << threads << "," << timer.get("ParPR-Block").count() << std::endl;

                if (f_pr != f_pr_block) {
                    std::cout << "flow not equal " << base_filename << " " << V(f_pr) << " " << V(f_pr_block) << std::endl;
                }
                 */
			}
		}

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
	whfc::runSnapshotTester(hgfile, threads);
	return 0;
}
