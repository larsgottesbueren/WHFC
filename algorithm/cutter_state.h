#pragma once

#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/isolated_nodes.h"
#include "../datastructure/bitset_reachable_sets.h"

namespace whfc {

	template<typename FlowAlgorithm>
	class CutterState {
	public:
		using Pin = FlowHypergraph::Pin;

		int viewDirection = 0;	//TODO wrap viewDirection in class?
		FlowHypergraph& hg;
		Flow flowValue = 0;

		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = FlowAlgorithm::ReachableNodes;
		//using ReachableHyperedges = FlowAlgorithm::ReachableHyperedges;
		ReachableNodes n;
		ReachableHyperedges h;
		std::vector<Node> sourcePiercingNodes, targetPiercingNodes;
		bool augmentingPathAvailableFromPiercing = true;
		bool hasCut = false;
		HyperedgeCut cut;
		NodeBorder borderNodes;
		IsolatedNodes isolatedNodes;
		NodeWeight maxBlockWeight;


		bool isBalanced() const {
			const NodeWeight 	sw = n.sourceReachableWeight,
								tw = n.targetReachableWeight,
								total = hg.totalNodeWeight(),
								iso = isolatedNodes.weight;

			//Beware: the isolated nodes are now weighted. Therefore the subset sum problem is no longer simple :(
			//but at least the instances are tiny
			return (sw <= maxBlockWeight && total - sw - iso <= maxBlockWeight) || (tw <= maxBlockWeight && total - tw - iso <= maxBlockWeight);
		}

		//what is this good for again?
		int moreBalancedCutSide() const {
			const NodeWeight
					sw = n.sourceReachableWeight,
					tw = n.targetReachableWeight,
					total = hg.totalNodeWeight(),
					iso = isolatedNodes.weight;

			//TODO incorporate iso
			const NodeWeight s_diff = std::max(maxBlockWeight - (total - sw), maxBlockWeight - sw);
			const NodeWeight t_diff = std::max(maxBlockWeight - (total - tw), maxBlockWeight - tw);
			return s_diff <= t_diff ? currentViewDirection() : oppositeViewDirection();
		}



		//TODO write functions to get the largest possible balanced node weight of a side, or maxBlockWeight+1 if that's not possible


		inline bool isIsolated(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && isolatedNodes.isCandidate(u); }
		inline bool canBeSettled(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && !isIsolated(u); }

		inline bool hasSourcePin(const Hyperedge e) const { return cut.hasSettledSourcePins[e]; }
		inline bool hasTargetPin(const Hyperedge e) const { return cut.hasSettledTargetPins[e]; }

		inline bool shouldBeAddedToCut(const Hyperedge e) const { return !h.areAllPinsSourceReachable(e) && !cut.wasAdded(e) && hg.isSaturated(e); }	//the first condition is just an optimization, not really necessary

		inline void addToCut(const Hyperedge e) {
			assert(shouldBeAddedToCut(e));
			for (const Pin& px : hg.pinsOf(e))
				if (canBeSettled(px.pin) && !borderNodes.wasAdded(px.pin))
					borderNodes.add(px.pin);
			cut.add(e);
		}

		void settleNode(const Node u) {
			assert(canBeSettled(u));
			if (!n.isSourceReachable(u))
				n.reach(u);
			n.settle(u);

			isolatedNodes.settleNode(u);

			for (const auto& he_inc : hg.hyperedgesOf(u)) {
				const Hyperedge e = he_inc.e;
				if (!hasSourcePin(e)) {
					cut.hasSettledSourcePins.set(e);
					if (hasTargetPin(e)) {	//e just became mixed
						isolatedNodes.accommodateNewlyMixedHyperedge(e, isIsolated);
					}
				}
			}
		}

		void flipViewDirection() {
			viewDirection = 1 - viewDirection;
			n.flipViewDirection();
			h.flipViewDirection();
			sourcePiercingNodes.swap(targetPiercingNodes);
			cut.flipViewDirection();
			borderNodes.flipViewDirection();
		}

		int currentViewDirection() const {
			return viewDirection;
		}

		int oppositeViewDirection() const {
			return 1 - viewDirection;
		}

		void clearForSearch() {
			if (augmentingPathAvailableFromPiercing) {
				n.resetSourceReachableToSource();
				h.resetSourceReachableToSource();
			}
		}

		void filterBorder() {
			borderNodes.filter([&](const Node& x) { return !canBeSettled(x);} );
		}

		void filterCut() {
			//not necessary at the moment
			//cut.deleteNonCutHyperedges(h);
		}

	};


}