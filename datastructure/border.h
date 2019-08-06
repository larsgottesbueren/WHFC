#pragma once

#include <vector>
#include "bitvector.h"
#include "../definitions.h"

namespace whfc {
	template<typename T, bool trackElements>
	class Border {
	public:
		explicit Border(const size_t nT) : addedToSourceSideBorder(nT), addedToTargetSideBorder(nT) { }

		void flipViewDirection() {
			std::swap(sourceSideBorder, targetSideBorder);
			std::swap(addedToSourceSideBorder, addedToTargetSideBorder);
		}

		//bool trackElements = true;
		BitVector addedToSourceSideBorder, addedToTargetSideBorder;
		std::vector<T> sourceSideBorder, targetSideBorder;

		inline bool wasAdded(const T x) const {
			return addedToSourceSideBorder[x];
		}

		inline void add(T x) {
			assert(!wasAdded(x));
			addedToSourceSideBorder.set(x);
			if constexpr (trackElements)
				sourceSideBorder.push_back(x);
		}

		inline void remove(size_t i) {
			if constexpr (trackElements) {
				sourceSideBorder[i] = sourceSideBorder.back();
				//std::swap(sourceSideBorder[i], sourceSideBorder.back());
				sourceSideBorder.pop_back();
			}
		}

		template<typename Predicate>
		void filter(Predicate pred) {
			if constexpr (trackElements) {
				auto new_end = std::remove_if(sourceSideBorder.begin(), sourceSideBorder.end(), pred);
				sourceSideBorder.erase(new_end, sourceSideBorder.end());
			}
		}
	};

	using NodeBorder = Border<Node, true>;

	class HyperedgeCut : public Border<Hyperedge, false> {
	public:
		using Base = Border<Hyperedge, false>;
		explicit HyperedgeCut(const size_t nHyperedges) : Base(nHyperedges) { }
		BitVector hasSettledSourcePins, hasSettledTargetPins;	//set in CutterState::settleNode
		size_t sourceMixed = 0, targetMixed = 0;	//equal if both cut-fronts were built. but they aren't.

		inline bool isHyperedgeMixed(const Hyperedge e) const {
			return hasSettledSourcePins[e] && hasSettledTargetPins[e];
		}

		/*

		 At the moment we don't need this

		//HyperedgeSet = FlowAlgorithm::ReachableHyperedges. Alternative: template the HyperedgeCut class and store reference
		template<class HyperedgeSet>
		void deleteNonCutHyperedges(const HyperedgeSet& h) {
			filter([&](const Hyperedge &e) {
				return isHyperedgeMixed(e) || h.areAllPinsSources(e);
			});
		}

		 */

		void flipViewDirection() {
			Base::flipViewDirection();
			std::swap(sourceMixed, targetMixed);
		}
	};
}