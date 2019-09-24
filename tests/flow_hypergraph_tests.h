#pragma once

#include "../io/hmetis_io.h"
#include "../logger.h"

namespace whfc {
namespace Test {
	
	class FlowHypergraphTests {
	public:
		static constexpr bool debug = true;
		void run() {
			std::string hgfile = "../test_hypergraphs/twocenters.hgr";
			FlowHypergraph hg = HMetisIO::readFlowHypergraph(hgfile);
			
			
			std::string output_hg_file = "../test_hypergraphs/twocenters_written.hgr";
			HMetisIO::writeFlowHypergraph(hg, output_hg_file);
		}
	};
}
}