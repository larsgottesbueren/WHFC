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

		int viewDirection = 0;
		FlowHypergraph& hg;
		Flow cutSize = 0;

		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = FlowAlgorithm::ReachableNodes;
		//using ReachableHyperedges = FlowAlgorithm::ReachableHyperedges;
		ReachableNodes n;
		ReachableHyperedges h;
		std::vector<Node> sourcePiercingNodes, targetPiercingNodes;
		bool augmentingPathAvailable = true;
		HyperedgeCut cut;
		NodeBorder borderNodes;
		IsolatedNodes isolatedNodes;

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

			for (const auto& he_inc : hg.hyperedgesOf(u)) {
				const Hyperedge e = he_inc.e;
				if (!hasSourcePin(e)) {
					cut.hasSettledSourcePins.set(e);
					if (hasTargetPin(e)) {	//e just became mixed
						for (const Pin& px : hg.pinsOf(e)) {
							const Node p = px.pin;
							isolatedNodes.mixedIncidentHyperedges[p]++;
							/*
							 * Previously, we identified candidates, via mixedIncidentHyperedges[p] == hg.degree(p), and later accepted them as isolated, only if they were not settled.
							 * This was particularly necessary for piercing hyperedges.
							 * However, it is suboptimal since the later settled candidates could just as well be moved between the blocks without modifying the cut.
							 * Immediately isolating them can yield better final cuts and is substantially less code.
							 * TODO this does not seem quite done yet.
							 * However, we have to move the ones that were marked reachable out of that set. Also make sure that it is IMPOSSIBLE to actually ever reach an isolated node.
							 */
							if (isIsolated(p))
								isolatedNodes.add(p);
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
			if (augmentingPathAvailable) {
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