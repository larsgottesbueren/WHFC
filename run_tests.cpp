#include "tests/subset_sum_tests.h"
#include "tests/flow_hypergraph_tests.h"

int main(int argc, char* argv[]) {
	Test::SubsetSumTests().run();
	Test::FlowHypergraphTests().run();

	return 0;
}