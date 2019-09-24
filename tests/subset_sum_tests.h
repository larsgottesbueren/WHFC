#pragma once
#include "create_dummy_hypergraphs.h"
#include "../datastructure/isolated_nodes.h"
#include "../datastructure/bitvector.h"

namespace whfc {
namespace Test {

	class SubsetSumTests {
	public:
		using NW = whfc::NodeWeight;
		using SR = whfc::IsolatedNodes::SummableRange;

		std::vector<NW> nodeWeights = {
			NW(2), NW(2), NW(3), NW(4), NW(2), NW(5), NW(3)
		};
		size_t mbw = 12;
		whfc::FlowHypergraph hg = CreateDummyHypergraphs::noHyperedgesGivenNodeWeights(nodeWeights);



		whfc::BitVector convertDPTableIntoBitvector(const whfc::IsolatedNodes& iso, const size_t ub) {
			whfc::BitVector res(static_cast<size_t>(ub+1));
			for (NW x(0); x <= ub; ++x) {
				if (iso.isSummable(x))
					res.set(static_cast<size_t>(x));
			}
			return res;
		}

							//Copy actual			//No need to copy expected
		bool compareRanges(std::vector<SR> actual, std::vector<SR>& expected) {
			std::sort(actual.begin(), actual.end());
			std::sort(expected.begin(), expected.end());
			return actual == expected;
		}


		void run() {
			whfc::IsolatedNodes iso(hg, NW::fromOtherValueType(mbw));
			whfc::BitVector expected(mbw+1);
			expected.set(0);

			AssertMsg(iso.isDPTableUpToDate(), "DP Table should be up to date without insertions");
			AssertMsg(convertDPTableIntoBitvector(iso, mbw) == expected, "No sums available yet");

			iso.updateDPTable();
			AssertMsg(iso.isDPTableUpToDate(), "DP Table should be up to date");
			AssertMsg(convertDPTableIntoBitvector(iso, mbw) == expected, "No sums available yet");
			AssertMsg(iso.getSumRanges().size() == 1 && iso.getSumRanges().front().from == 0 && iso.getSumRanges().front().to == 0, "Only range should be the one containing zero");

			iso.add(whfc::Node(2));
			iso.updateDPTable();
			expected.set(3);
			//items = [3] --> sums = [0,3]
			AssertMsg(iso.isSummable(NW(3)), "3 should be summable. Nothing else");
			AssertMsg(convertDPTableIntoBitvector(iso, mbw) == expected, "Check");
			std::vector<SR> expectedRanges = { SR(NW(0), NW(0)), SR(NW(3),NW(3)) };
			AssertMsg(compareRanges(iso.getSumRanges(), expectedRanges), "Two ranges, second one should be 3 only");
			
			iso.add(whfc::Node(1));
			expected.set(2);
			expected.set(5);
			//items = [3,2] --> sums = [0,2,3,5]
			//tests extend range on the left
			iso.updateDPTable();
			AssertMsg(convertDPTableIntoBitvector(iso, mbw) == expected, "Check");
			expectedRanges = { SR(NW(0), NW(0)), SR(NW(2), NW(3)), SR(NW(5), NW(5)) };
			AssertMsg(compareRanges(iso.getSumRanges(), expectedRanges), "Should be 3 ranges. [0], [2,3], [5]");
			
			iso.add(whfc::Node(5));
			iso.updateDPTable();
			//items = [3,2,5] --> sums = [0,2,3,5,7,8,10]
			//tests extend range on the right (7 will be summed before 8)
			expected.set(7);
			expected.set(8);
			expected.set(10);
			Assert(convertDPTableIntoBitvector(iso, mbw) == expected);
			expectedRanges = { SR(NW(0), NW(0)), SR(NW(2), NW(3)), SR(NW(5), NW(5)), SR(NW(7), NW(8)), SR(NW(10), NW(10)) };
			AssertMsg(compareRanges(iso.getSumRanges(), expectedRanges), "Should be 5 ranges. [0], [2,3], [5], [7,8], [10]");

			iso.add(whfc::Node(0));
			iso.updateDPTable();
			//items = [2,2,3,5] --> sums = [0,2,3,4,5,7,8,9,10,12]
			//tests range merging
			expected.set(4);
			expected.set(9);
			expected.set(12);
			AssertMsg(convertDPTableIntoBitvector(iso, mbw) == expected, "Check");
			expectedRanges = { SR(NW(0), NW(0)), SR(NW(2), NW(5)), SR(NW(7), NW(10)),  SR(NW(12), NW(12)) };
			AssertMsg(compareRanges(iso.getSumRanges(), expectedRanges), "Should be 4 ranges. [0], [2,5], [7-10], [12]");
		};

	};
}
}