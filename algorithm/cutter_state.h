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

		inline bool isIsolated(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && isolatedNodes.isCandidate(u); }
		inline bool canBeSettled(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && !isIsolated(u); }

		inline NodeWeight unclaimedNodeWeight() const {
			return hg.totalNodeWeight() - n.sourceReachableWeight - n.targetReachableWeight - isolatedNodes.weight;
		}

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


		//Due to isolated nodes, this function is SUPER involved. If you're reimplementing this in a different context, I suggest you just settle isolated nodes.
		//This function contains a lot of premature optimization. I have no clue how bad (or rather not) things would get without these little optimizations.
		bool isBalanced() {

			const NodeWeight
					sw = n.sourceReachableWeight,
					tw = n.targetReachableWeight,
					uw = unclaimedNodeWeight(),
					total = hg.totalNodeWeight(),
					iso = isolatedNodes.weight,
					suw = sw + uw,
					tuw = tw + uw;

			{	//quick checks to determine early that balance is not possible
				//even though checking for balance is SubsetSum, once we are roughly in the realm of balance, it is likely that a solution exists
				//in order to save running time, we therefore want to make sure to invoke SubsetSum as rarely as possible

				if (sw > maxBlockWeight || tw > maxBlockWeight)					//this is good at late and early stages
					return false;

				if (sw + uw > maxBlockWeight && tw + uw > maxBlockWeight)		//this is good at early stages
					return false;
				//find some more!
			}

			{	//quick checks to determine whether balance is possible without invoking SubsetSum, i.e. don't split the isolated nodes.
				//this should be possible often enough.
				bool balanced = false;
				balanced |= sw + uw + iso <= maxBlockWeight;
				balanced |= tw + uw + iso <= maxBlockWeight;
				balanced |= sw + uw <= maxBlockWeight && tw + iso <= maxBlockWeight;
				balanced |= tw + uw <= maxBlockWeight && sw + iso <= maxBlockWeight;
				if (balanced)
					return true;


				//TODO figure out some ideas for when there are very homogenous node weights, in particular for the fine levels of the ML hierarchy.
			}

			{	//reuse cached values from the SubsetSum invokation.
				//requires that either no new isolated node was added
				//other possibilities are: store the smallest subset-sum since last update/query

			}

			if (tw + uw > maxBlockWeight) {
				//then sw must take uw

			}

			isolatedNodes.updateDPTable();


			//TODO incorporate iso
			const NodeWeight s_diff = std::max(maxBlockWeight - (total - sw), maxBlockWeight - sw);
			const NodeWeight t_diff = std::max(maxBlockWeight - (total - tw), maxBlockWeight - tw);
			return false;
		}

		//TODO write functions to get the largest possible balanced node weight of a side, or maxBlockWeight+1 if that's not possible

	};


}