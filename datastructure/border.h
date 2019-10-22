#pragma once

#include <vector>
#include "bitvector.h"
#include "../definitions.h"

#include "../util/filter.h"
#include "flow_hypergraph.h"

namespace whfc {
	template<typename T, bool trackElements>
	class Border {
	public:
		explicit Border(size_t nT) : addedToSourceSideBorder(nT), addedToTargetSideBorder(nT) { }

		void flipViewDirection() {
			std::swap(sourceSideBorder, targetSideBorder);
			std::swap(addedToSourceSideBorder, addedToTargetSideBorder);
		}

		BitVector addedToSourceSideBorder, addedToTargetSideBorder;
		std::vector<T> sourceSideBorder, targetSideBorder;

		inline bool wasAdded(const T x) const {
			return addedToSourceSideBorder[x];
		}

		inline void add(T x) {
			Assert(!wasAdded(x));
			addedToSourceSideBorder.set(x);
			if constexpr (trackElements)
				sourceSideBorder.push_back(x);
		}

		inline void remove(size_t i) {
			if constexpr (trackElements) {
				sourceSideBorder[i] = sourceSideBorder.back();
				sourceSideBorder.pop_back();
			}
		}

		template<typename Predicate>
		void cleanUp(Predicate pred) {
			if constexpr (trackElements) {
				util::remove_if_inplace(sourceSideBorder, pred);
			}
		}
		
		void reset(const size_t newN) {
			addedToSourceSideBorder.reset(0, newN);
			addedToTargetSideBorder.reset(0, newN);
			sourceSideBorder.clear();
			targetSideBorder.clear();
		}
	
	};

	using NodeBorder = Border<Node, true>;

	//track hyperedges only for assertions in debug mode
#ifndef NDEBUG
	using HyperedgeCutBase = Border<Hyperedge, true>;
#else
	using HyperedgeCutBase = Border<Hyperedge, false>;
#endif
	
	class HyperedgeCut : public HyperedgeCutBase {
	public:
		using Base = HyperedgeCutBase;
		explicit HyperedgeCut(const size_t nHyperedges) : Base(nHyperedges), hasSettledSourcePins(nHyperedges), hasSettledTargetPins(nHyperedges) { }
		BitVector hasSettledSourcePins, hasSettledTargetPins;	//set in CutterState::settleNode //TODO check if hasSettledSourcePins == ReachableHyperedges::areFlowSendingPinsSources()

		inline bool isHyperedgeMixed(const Hyperedge e) const {
			return hasSettledSourcePins[e] && hasSettledTargetPins[e];
		}

		
		Flow weight(const FlowHypergraph& hg) const {
			Flow w = 0;
			for (const Hyperedge e : sourceSideBorder) {
				Assert(hg.isSaturated(e));
				w += hg.capacity(e);
			}
			return w;
		}

		void flipViewDirection() {
			Base::flipViewDirection();
			std::swap(hasSettledSourcePins, hasSettledTargetPins);
		}
		
		void reset(const size_t newN) {
			Base::reset(newN);
			hasSettledSourcePins.reset(0, newN);
			hasSettledTargetPins.reset(0, newN);
		}
		
	};
}
