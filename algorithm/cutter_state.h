#pragma once

#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/isolated_nodes.h"
#include "../datastructure/bitset_reachable_sets.h"

namespace whfc {

	template<typename FlowAlgorithm>
	class CutterState {
		static constexpr bool debug = true;
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
		NodeWeight maxBlockWeight;
		IsolatedNodes isolatedNodes;

		CutterState(FlowHypergraph& _hg, NodeWeight _maxBlockWeight) :
				hg(_hg),
				n(_hg),
				h(static_cast<size_t>(_hg.numHyperedges())),
				cut(static_cast<size_t>(_hg.numHyperedges())),
				borderNodes(static_cast<size_t>(_hg.numNodes())),
				maxBlockWeight(_maxBlockWeight),
				isolatedNodes(hg, _maxBlockWeight)
		{

		}

		inline bool isIsolated(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && isolatedNodes.isCandidate(u); }
		inline bool canBeSettled(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && !isIsolated(u); }

		inline bool isInfeasible() {
			if (hg.totalNodeWeight() > 2 * maxBlockWeight)
				return true;

			LOG << "Starting expensive SubsetSum based feasibility check";
			//TODO make IsolatedNodes reusable in a way that we can copy isolatedNodes and add the unclaimed nodes here.
			return false;
		}

		inline NodeWeight unclaimedNodeWeight() const { return hg.totalNodeWeight() - n.sourceReachableWeight - n.targetReachableWeight - isolatedNodes.weight; }
		inline bool hasSourcePin(const Hyperedge e) const { return cut.hasSettledSourcePins[e]; }
		inline bool hasTargetPin(const Hyperedge e) const { return cut.hasSettledTargetPins[e]; }
		inline bool shouldBeAddedToCut(const Hyperedge e) const { return !h.areAllPinsSourceReachable(e) && !cut.wasAdded(e) && hg.isSaturated(e); }	//the first condition is just an optimization, not really necessary

		inline void addToCut(const Hyperedge e) {
			Assert(shouldBeAddedToCut(e));
			for (const Pin& px : hg.pinsOf(e))
				if (canBeSettled(px.pin) && !borderNodes.wasAdded(px.pin))
					borderNodes.add(px.pin);
			cut.add(e);
		}

		void settleNode(const Node u) {
			Assert(canBeSettled(u));
			if (!n.isSourceReachable(u))
				n.reach(u);
			n.settle(u);

			isolatedNodes.settleNode(u);

			for (const auto& he_inc : hg.hyperedgesOf(u)) {
				const Hyperedge e = he_inc.e;
				if (!hasSourcePin(e)) {
					cut.hasSettledSourcePins.set(e);
					if (hasTargetPin(e)) {	//e just became mixed
						for (const auto& px : hg.pinsOf(e)) {
							const Node p = px.pin;
							isolatedNodes.mixedIncidentHyperedges[p]++;
							/*
							 * Previously, we identified candidates, via mixedIncidentHyperedges[p] == hg.degree(p), and later accepted them as isolated, only if they were not settled.
							 * This was particularly necessary for piercing hyperedges.
							 * However, it is suboptimal since the later settled candidates could just as well be moved between the blocks without modifying the cut.
							 * Immediately isolating them can yield balance faster, requires substantially less code, and should be easier to understand.
							 * TODO this does not seem quite done yet. think about it and then come back.
							 * However, we have to move the ones that were marked reachable out of that set. Also make sure that it is IMPOSSIBLE to actually ever reach an isolated node.
							 *
							 * Doing this incurs another optimization problem. In which order should nodes be settled to create as many isolated nodes as possible?
							 *
							 * We also have to stop the growAssimilated code from settling these candidates.
							 *
							 * Additionally, the smaller side might change. This shouldn't matter much but is quite weird.
							 *
							 */
							if (isIsolated(p)) {

								isolatedNodes.add(p);
							}
						}
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
		//Same goes for the SubsetSum solver in datastructure/isolated_nodes.h
		bool isBalanced() {

			const NodeWeight
					sw = n.sourceReachableWeight,		//cannot be split
					tw = n.targetReachableWeight,		//cannot be split
					uw = unclaimedNodeWeight(),			//cannot be split (in current stages. if we integrate proper PCKP heuristics for MBMC this would change)
					total = hg.totalNodeWeight(),
					iso = isolatedNodes.weight;			//can be split


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
			}

			{	//reuse cached values from the previous SubsetSum invokation.
				//requires that either no new isolated node was added
				//other possibilities are: store the smallest subset-sum since last update/query

			}

			isolatedNodes.updateDPTable();

			const NodeWeight
					sRem = maxBlockWeight - sw,
					tRem = maxBlockWeight - tw,
					suw = sw + uw,
					tuw = tw + uw,
					suwRem = suw <= maxBlockWeight ? maxBlockWeight - suw : NodeWeight::Invalid(),
					tuwRem = tuw <= maxBlockWeight ? maxBlockWeight - tuw : NodeWeight::Invalid();

			//instead of iterating over the ranges, we could do a segment query. that is hopefully overkill!

			//sides: (S + U, T) + <ISO> and (S, T + U) + <ISO>
			for (const IsolatedNodes::SummableRange& sr : isolatedNodes.getSumRanges()) {
				if (suwRem.isValid()) {
					//S+U not overloaded. Therefore, try (S + U, T) + <ISO>

					//allocate as much as possible to S+U, i.e. x = min(suwRem, sr.to), the rest, i.e. iso - x has to go to T
					if (suwRem >= sr.from && tw + (iso - std::min(suwRem, sr.to)) <= maxBlockWeight)
							return true;
					//analogously, allocate as much as possible to T
					if (tRem >= sr.from && suw + (iso - std::min(tRem, sr.to)) <= maxBlockWeight)
						return true;
				}

				if (tuwRem.isValid()) {
					//T+U not overloaded. Therefore, try (S, T + U) + <ISO>
					if (tuwRem >= sr.from && sw + (iso - std::min(tuwRem, sr.to)) <= maxBlockWeight)
						return true;
					if (sRem >= sr.from && tuw + (iso - std::min(sRem, sr.to)) <= maxBlockWeight)
						return true;
				}
			}

			return false;
		}


		//TODO should we cache result of isBalanced() to reuse it easily?
		//TODO write functions to get the largest possible balanced node weight of a side, or maxBlockWeight+1 if that's not possible
	};


}