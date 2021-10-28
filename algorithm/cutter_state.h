#pragma once

#include "../datastructure/queue.h"
#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/node_border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../util/math.h"
#include "../util/random.h"


// TODO factor out verification code, maybe even balance checking?

namespace whfc {

	struct SimulatedNodeAssignment {
		bool assignUnclaimedToSource = true;
		double imbalanceSourceBlock = std::numeric_limits<double>::max(), imbalanceTargetBlock = std::numeric_limits<double>::max();
		size_t numberOfTrackedMoves = 0;
		int direction = 0;

		double imbalance() const {
			return std::max(imbalanceSourceBlock, imbalanceTargetBlock);
		}

		bool isPerfectlyBalanced() const {
			return std::abs(imbalanceSourceBlock - imbalanceTargetBlock) < 1e-9;
		}
	};

	struct Move {
		enum class Type : uint8_t { SettleNode, SettleAllPins, SettleFlowSendingPins };
		Node node;
		Hyperedge hyperedge;
		int direction;
		Type t;
		Move(Node node, int dir, Type t) : node(node), hyperedge(invalidHyperedge), direction(dir), t(t) {
			assert(t == Type::SettleNode);
		}
		Move(Hyperedge hyperedge, int dir, Type t) : node(invalidNode), hyperedge(hyperedge), direction(dir), t(t) {
			assert(t == Type::SettleAllPins || t == Type::SettleFlowSendingPins);
		}
	};

	struct PiercingNode {
		Node node;
		bool isReachableFromOppositeSide;
		PiercingNode(const Node node, bool isReachableFromOppositeSide) : node(node), isReachableFromOppositeSide(isReachableFromOppositeSide) { }
	};

	struct NonDynamicCutterState {
		std::vector<PiercingNode> sourcePiercingNodes, targetPiercingNodes;
		int direction;
	};

	template<typename FlowAlgorithm>
	class CutterState {
	public:
		static constexpr bool log = false;

		using Pin = FlowHypergraph::Pin;

		FlowAlgorithm flow_algo;
		int viewDirection = 0;
		FlowHypergraph& hg;
		Flow flowValue = 0;

		using ReachableNodes = typename FlowAlgorithm::ReachableNodes;
		using ReachableHyperedges = typename FlowAlgorithm::ReachableHyperedges;
		ReachableNodes n;
		ReachableHyperedges h;
		std::vector<PiercingNode> sourcePiercingNodes, targetPiercingNodes;
		std::vector<Move> trackedMoves;

		bool augmentingPathAvailableFromPiercing = true;
		bool hasCut = false;
		bool mostBalancedCutMode = false;
		HyperedgeCuts cuts;
		NodeBorders borderNodes;
		std::array<NodeWeight, 2> maxBlockWeightPerSide;
		bool partitionWrittenToNodeSet = false;
		TimeReporter& timer;
		Randomizer rng;

		CutterState(FlowHypergraph& _hg, TimeReporter& timer) :
				flow_algo(_hg),
				hg(_hg),
				n(_hg),
				h(_hg),
				cuts(_hg.numHyperedges()),
				borderNodes(_hg.numNodes()),
				maxBlockWeightPerSide({NodeWeight(0), NodeWeight(0)}),
				timer(timer)
		{
			timer.registerCategory("Balance Check");
		}

		inline bool canBeSettled(const Node u) const {
			return !n.isSource(u) && !n.isTarget(u);
		}

		inline NodeWeight unclaimedNodeWeight() const {
			return hg.totalNodeWeight() - n.sourceReachableWeight - n.targetReachableWeight;
		}

		inline NodeWeight notSettledNodeWeight() const {
			return hg.totalNodeWeight() - n.sourceWeight - n.targetWeight;
		}

		inline bool shouldBeAddedToCut(const Hyperedge e) const {
			return !h.areAllPinsSourceReachable(e) && !cuts.sourceSide.wasAdded(e) && hg.isSaturated(e); // the first condition is just an optimization, not really necessary
		}

		inline void addToCut(const Hyperedge e) {
			//Note: the current implementation of selecting piercing nodes relies on not inserting target-reachable nodes during most balanced cut mode
			assert(shouldBeAddedToCut(e));
			for (const Pin& px : hg.pinsOf(e)) {
				if (canBeSettled(px.pin) && !borderNodes.sourceSide->wasAdded(px.pin) && (!mostBalancedCutMode || !n.isTargetReachable(px.pin))) {
					borderNodes.sourceSide->add(px.pin, n.isTargetReachable(px.pin));
				}
			}
			cuts.sourceSide.add(e);
		}

		void setMaxBlockWeight(int side, NodeWeight mw) {
			maxBlockWeightPerSide[side] = mw;
		}

		NodeWeight maxBlockWeight(int side) const {
			return maxBlockWeightPerSide[side];
		}

		NodeWeight maxBlockWeight() const {
			return maxBlockWeight(currentViewDirection());
		}


		void settleNode(const Node u, bool check = true) {
			assert(!n.isSource(u) && !n.isTarget(u));
			unused(check);

			if (!n.isSourceReachable(u))
				n.reach(u);
			n.settle(u);

			if (mostBalancedCutMode) {
				trackedMoves.emplace_back(u, currentViewDirection(), Move::Type::SettleNode);
				return;
			}

		}

		// TODO rename to report settling flow sending pins!
		void settleFlowSendingPins(const Hyperedge e) {
			if (mostBalancedCutMode) {
				trackedMoves.emplace_back(e, currentViewDirection(), Move::Type::SettleFlowSendingPins);
			}
			h.settleFlowSendingPins(e);
		}

		void settleAllPins(const Hyperedge e) {
			if (mostBalancedCutMode) {
				trackedMoves.emplace_back(e, currentViewDirection(), Move::Type::SettleAllPins);
			}
			h.settleAllPins(e);
		}

		void flipViewDirection() {
			viewDirection = 1 - viewDirection;
			hg.flipViewDirection();
			n.flipViewDirection();
			h.flipViewDirection();
			sourcePiercingNodes.swap(targetPiercingNodes);
			cuts.flipViewDirection();
			borderNodes.flipViewDirection();
		}

		int currentViewDirection() const {
			return viewDirection;
		}

		int oppositeViewDirection() const {
			return 1 - viewDirection;
		}

		void clearForSearch() {
			n.resetSourceReachableToSource(augmentingPathAvailableFromPiercing);
			h.resetSourceReachableToSource(augmentingPathAvailableFromPiercing);
		}

		void reset() {		// TODO could consolidate with initialize
			viewDirection = 0;
			flowValue = 0;
			flow_algo.reset();
			n.fullReset();
			h.fullReset();
			sourcePiercingNodes.clear(); targetPiercingNodes.clear();
			trackedMoves.clear();
			augmentingPathAvailableFromPiercing = true;
			hasCut = false;
			mostBalancedCutMode = false;
			cuts.reset(hg.numHyperedges());			//this requires that FlowHypergraph is reset before resetting the CutterState
			borderNodes.reset(hg.numNodes());
			partitionWrittenToNodeSet = false;
		}

		void initialize(const Node s, const Node t) {
			if (hg.nodeWeight(s) > maxBlockWeight(0) || hg.nodeWeight(t) > maxBlockWeight(1)) {
				throw std::runtime_error("Terminal weight already exceeds max block weight at initialization. Consider setting max block weights per side via hfc.cs.setMaxBlockWeight(  side  )");
			}
			assert(sourcePiercingNodes.empty() && targetPiercingNodes.empty());
			sourcePiercingNodes.emplace_back(s,false);
			settleNode(s, false);
			targetPiercingNodes.emplace_back(t,false);
			flipViewDirection();
			settleNode(t, false);
			flipViewDirection();
		}

		int sideToGrow() const {
			const double imb_s = static_cast<double>(n.sourceReachableWeight) / static_cast<double>(maxBlockWeight(currentViewDirection()));
			const double imb_t = static_cast<double>(n.targetReachableWeight) / static_cast<double>(maxBlockWeight(oppositeViewDirection()));
			return imb_s <= imb_t ? currentViewDirection() : oppositeViewDirection();
		}

		bool isBalanced() {
			assert(hasCut);
			assert(!partitionWrittenToNodeSet && "Cannot call isBalanced() once the partition has been written");

			const NodeWeight
					sw = n.sourceReachableWeight,		//cannot be split
					tw = n.targetReachableWeight,		//cannot be split
					uw = unclaimedNodeWeight();			//cannot be split (in current stages. if we integrate proper PCKP heuristics for MBMC this would change)

			const NodeWeight
					s_mbw = maxBlockWeight(currentViewDirection()),
					t_mbw = maxBlockWeight(oppositeViewDirection());

			if (sw > s_mbw || tw > t_mbw)					//this is good at late and early stages
				return false;
			if (sw + uw > s_mbw && tw + uw > t_mbw)			//this is good at early stages
				return false;


			bool balanced = false;
			balanced |= sw + uw <= s_mbw && tw <= t_mbw;
			balanced |= tw + uw <= t_mbw && sw <= s_mbw;
			return balanced;
		}

		NonDynamicCutterState enterMostBalancedCutMode() {
			assert(!mostBalancedCutMode);
			assert(trackedMoves.empty());
			assert(hasCut);
			mostBalancedCutMode = true;	// activates move tracking
			borderNodes.enterMostBalancedCutMode();
			cuts.enterMostBalancedCutMode();
			return { sourcePiercingNodes, targetPiercingNodes, currentViewDirection() };
		}

		void resetToFirstBalancedState(NonDynamicCutterState& nds) {
			if (currentViewDirection() != nds.direction) {
				flipViewDirection();
			}
			sourcePiercingNodes = nds.sourcePiercingNodes;
			targetPiercingNodes = nds.targetPiercingNodes;
			revertMoves(0);
			borderNodes.resetForMostBalancedCut();
			cuts.resetForMostBalancedCut();
		}

		SimulatedNodeAssignment mostBalancedAssignment() {
			auto block_imb = [&](NodeWeight a, NodeWeight max_a) {
				return (static_cast<double>(a) / static_cast<double>(max_a)) - 1.0;
			};
			SimulatedNodeAssignment suw;
			suw.imbalanceSourceBlock = block_imb(hg.totalNodeWeight() - n.targetReachableWeight, maxBlockWeight(currentViewDirection()));
			suw.imbalanceTargetBlock = block_imb(n.targetReachableWeight, maxBlockWeight(oppositeViewDirection()));

			SimulatedNodeAssignment tuw;
			tuw.imbalanceSourceBlock = block_imb(n.sourceReachableWeight, maxBlockWeight(currentViewDirection()));
			tuw.imbalanceTargetBlock = block_imb(hg.totalNodeWeight() - n.sourceReachableWeight(), maxBlockWeight(oppositeViewDirection()));

			SimulatedNodeAssignment sol = suw.imbalance() < tuw.imbalance() ? suw : tuw;

			sol.numberOfTrackedMoves = trackedMoves.size();
			sol.direction = currentViewDirection();
			return sol;
		}

		// takes the information from mostBalancedIsolatedNodesAssignment()
		// can be an old run, since the DP solution for trackedIsolatedWeight only contains nodes that were isolated during that run
		void writePartition(const SimulatedNodeAssignment& r) {
			assert(!partitionWrittenToNodeSet);
			assert(isBalanced());
			if (currentViewDirection() != r.direction)
				flipViewDirection();

			for (const Node u : hg.nodeIDs()) {
				if (n.isSourceReachable(u) && !n.isSource(u))
					n.settle(u);

				if (n.isTargetReachable(u) && !n.isTarget(u))
					n.settleTarget(u);

				if (!n.isSourceReachable(u) && !n.isTargetReachable(u)) {
					if (r.assignUnclaimedToSource) {
						n.reach(u); n.settle(u);
					}
					else {
						n.reachTarget(u); n.settleTarget(u);
					}
				}
			}

			if (currentViewDirection() != 0)
				flipViewDirection();

			assert(n.sourceWeight + n.targetWeight == hg.totalNodeWeight());
			partitionWrittenToNodeSet = true;
		}

		void writePartition() {
			writePartition(mostBalancedAssignment());
		}

		void revertMoves(const size_t numberOfTrackedMoves) {
			while (trackedMoves.size() > numberOfTrackedMoves) {
				Move& m = trackedMoves.back();
				if (m.node != invalidNode) {
					assert(m.hyperedge == invalidHyperedge);
					if (m.direction == currentViewDirection())
						n.unsettleSource(m.node);
					else
						n.unsettleTarget(m.node);
				}
				else {
					assert(m.node == invalidNode);
					//for timestamp and distance reachable sets, we would only need unsettleAllPins and unsettleFlowSendingPins, since S and T are disjoint by nature.
					if (currentViewDirection() == m.direction) {
						if (m.t == Move::Type::SettleAllPins)
							h.unsettleAllPins(m.hyperedge);
						else
							h.unsettleFlowSendingPins(m.hyperedge);
					}
					else {
						if (m.t == Move::Type::SettleAllPins)
							h.unsettleAllPinsTarget(m.hyperedge);
						else
							h.unsettleFlowSendingPinsTarget(m.hyperedge);
					}
				}
				trackedMoves.pop_back();
			}
		}

		void applyMoves(const std::vector<Move>& moves) {
			for (const Move& m : moves) {
				if (m.node != invalidNode) {
					assert(!n.isSourceReachable(m.node) && !n.isTargetReachable(m.node));
					if (m.direction == currentViewDirection()) {
						n.reach(m.node); n.settle(m.node);
					}
					else {
						n.reachTarget(m.node); n.settleTarget(m.node);
					}
				}
				else {
					assert(m.hyperedge != invalidHyperedge);
					if (currentViewDirection() == m.direction) {
						if (m.t == Move::Type::SettleAllPins)
							h.settleAllPins(m.hyperedge);
						else
							h.settleFlowSendingPins(m.hyperedge);
					}
					else {
						if (m.t == Move::Type::SettleAllPins)
							h.settleAllPinsTarget(m.hyperedge);
						else
							h.settleFlowSendingPinsTarget(m.hyperedge);
					}
				}
			}
		}

		std::string toString(bool skip_iso_and_unclaimed = false) {
			std::stringstream os;
			os << " cut= " << flowValue
			   << " s=" << n.sourceWeight << "|" << n.sourceReachableWeight
			   << " t=" << n.targetWeight << "|" << n.targetReachableWeight;
			if (!skip_iso_and_unclaimed)
			   os << " u=" << unclaimedNodeWeight();
			os << " mbw=[" << maxBlockWeight(currentViewDirection()) << " " << maxBlockWeight(oppositeViewDirection()) << "]"
			   << " total=" << hg.totalNodeWeight()
			   ;
			return os.str();
		}



		void verifyCutPostConditions() {
			assert(hasCut);

#ifndef NDEBUG
			cuts.sourceSide.cleanUp([&](const Hyperedge& e) { return h.areAllPinsSources(e); });
			Flow expected_flow = 0;
			for (const Hyperedge& e : cuts.sourceSide.entries()) {
				assert(hg.isSaturated(e));
				expected_flow += hg.capacity(e);
			}
			assert(flowValue == expected_flow);


#endif
			verifyExtractedCutHyperedgesActuallySplitHypergraph();
			verifyCutInducedByPartitionMatchesExtractedCutHyperedges();
		}

		void verifyCutInducedByPartitionMatchesExtractedCutHyperedges() {
#ifndef NDEBUG
			std::vector<Hyperedge> cut_from_partition;
			for (Hyperedge e : hg.hyperedgeIDs()) {
				bool hasSource = false;
				bool hasOther = false;
				for (Pin& p : hg.pinsOf(e)) {
					Node v = p.pin;
					hasSource |= n.isSource(v);
					hasOther |= !n.isSource(v);
				}
				if (hasSource && hasOther) {
					cut_from_partition.push_back(e);
					assert(h.areFlowSendingPinsSources(e));
				}

				if (hasSource && !hasOther)
					assert(h.areAllPinsSources(e));
			}
			std::vector<Hyperedge> sorted_cut = cuts.sourceSide.copy();
			std::sort(sorted_cut.begin(), sorted_cut.end());
			assert(sorted_cut == cut_from_partition);
#endif
		}

		void verifyCutInducedByPartitionMatchesFlowValue() {
#ifndef NDEBUG
			Flow cut_weight = 0;
			for (Hyperedge e : hg.hyperedgeIDs()) {
				bool hasSource = false;
				bool hasOther = false;
				for (Pin& p : hg.pinsOf(e)) {
					hasSource |= n.isSource(p.pin);
					hasOther |= !n.isSource(p.pin);
				}
				if (hasSource && hasOther) {
					cut_weight += hg.capacity(e);
				}
			}
			assert(flowValue == cut_weight);
#endif
		}


		void verifyExtractedCutHyperedgesActuallySplitHypergraph() {
#ifndef NDEBUG
			BitVector he_seen(hg.numHyperedges()), node_seen(hg.numNodes());
			LayeredQueue<Node> queue(hg.numNodes());
			for (Node u : hg.nodeIDs()) {
				if (n.isSource(u)) {
					queue.push(u);
					node_seen.set(u);
				}
			}

			for (Hyperedge e : cuts.sourceSide.entries())
				he_seen.set(e);

			while (!queue.empty()) {
				Node u = queue.pop();
				for (auto& he_inc : hg.hyperedgesOf(u)) {
					Hyperedge e = he_inc.e;
					if (!he_seen[e]) {
						he_seen.set(e);
						for (auto& pin : hg.pinsOf(e)) {
							Node v = pin.pin;
							assert(!n.isTargetReachable(v));
							assert(n.isSourceReachable(v));
							if (!node_seen[v]) {
								node_seen.set(v);
								queue.push(v);
							}
						}
					}
				}
			}

			for (Node u : hg.nodeIDs()) {
				if (n.isTargetReachable(u))
					assert(!node_seen[u]);
				if (n.isSourceReachable(u))
					assert(node_seen[u]);
			}

			queue.clear();
			he_seen.reset();
			node_seen.reset();
			for (Node u : hg.nodeIDs()) {
				if (n.isTarget(u)) {
					queue.push(u);
					node_seen.set(u);
				}
			}

			for (Hyperedge e : cuts.sourceSide.entries())
				he_seen.set(e);

			while (!queue.empty()) {
				Node u = queue.pop();
				for (auto& he_inc : hg.hyperedgesOf(u)) {
					Hyperedge e = he_inc.e;
					if (!he_seen[e]) {
						he_seen.set(e);
						for (auto& pin : hg.pinsOf(e)) {
							Node v = pin.pin;
							assert(!n.isSourceReachable(v));
							//no assert(n.isTargetReachable(v)) since we removed the source-side cut
							if (!node_seen[v]) {
								node_seen.set(v);
								queue.push(v);
							}
						}
					}
				}
			}

			for (Node u : hg.nodeIDs()) {
				if (n.isTargetReachable(u))
					assert(node_seen[u]);
				if (n.isSourceReachable(u))
					assert(!node_seen[u]);
			}
#endif
		}

	};

}
