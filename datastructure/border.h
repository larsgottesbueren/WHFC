#pragma once

#include <vector>
#include "bitvector.h"
#include "../definitions.h"

namespace whfc {
	template<typename T>
	class Border {
	public:
		explicit Border(size_t nT) : addedToSourceSide(nT), addedToTargetSide(nT) { }

		void flipViewDirection() {
			std::swap(sourceSide, targetSide);
			std::swap(addedToSourceSide, addedToTargetSide);
		}

		bool trackElements = true;
		BitVector addedToSourceSide, addedToTargetSide;
		std::vector<T> sourceSide, targetSide;

		inline bool wasAdded(const T x) const { return addedToSourceSide[x]; }

		inline void add(T x) {
			assert(!wasAdded(x));
			addedToSourceSide.set(x);
			if (trackElements)
				sourceSide.push_back(x);
		}
		inline void remove(size_t i) {
			std::swap(sourceSide[i], sourceSide.back());
			sourceSide.pop_back();
		}
	};

	using NodeBorder = Border<Node>;

	class HyperedgeCut : public Border<Hyperedge> {
	public:
		using Base = Border<Hyperedge>;
		explicit HyperedgeCut(size_t nHyperedges) : Base(nHyperedges) { trackElements = false; }
		BitVector hasSettledSourcePins, hasSettledTargetPins;	//set in CutterState::settleNode
		size_t sourceMixed = 0, targetMixed = 0;	//equal if both cut-fronts were built. but they aren't.

		inline bool isHyperedgeMixed(const Hyperedge e) const { return hasSettledSourcePins[e] && hasSettledTargetPins[e]; }

		//HyperedgeSet = FlowAlgorithm::ReachableHyperedges. Alternative: template the HyperedgeCut class and store reference
		template<class HyperedgeSet>
		void deleteNonCutHyperedges(const HyperedgeSet& h) {
			for (size_t i = 0; i < sourceSide.size(); i++) {
				Hyperedge e = sourceSide[i];
				if (isHyperedgeMixed(e) || h.areAllPinsSources(e)) {
					sourceMixed += static_cast<size_t>(isHyperedgeMixed(e));
					remove(i);
					i--;
				}
			}
		}

		void flipViewDirection() {
			Base::flipViewDirection();
			std::swap(sourceMixed, targetMixed);
		}
	};
}