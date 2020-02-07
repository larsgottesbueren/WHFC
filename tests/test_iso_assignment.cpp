#include "../algorithm/cutter_state.h"
#include "../algorithm/dinic.h"

namespace whfc {
	NodeWeight parse(const char* s) {
		return NodeWeight::fromOtherValueType(std::stoul(s));
	}
	
	void testIsoAssignment(int argc, const char* argv[]) {
		if (argc != 7) {
			std::cout << "Usage ./test_iso_assignment a max_a b max_b sr.from sr.to" << std::endl;
			std::exit(-1);
		}
		
		std::array<NodeWeight, 6> input;
		for (size_t i = 0; i < input.size(); ++i) {
			input[i] = parse(argv[i + 1]);
		}
		
		auto [x_iso, imbalance] = CutterState<Dinic>::isolatedWeightAssignmentToFirstMinimizingImbalance(input[0], input[1], input[2], input[3], IsolatedNodes::SummableRange(input[4], input[5]));
		std::cout << "x_iso= " << x_iso << " imbalance= " << imbalance << std::endl;
	}
}

int main(int argc, const char* argv[]) {
	whfc::testIsoAssignment(argc, argv);
}