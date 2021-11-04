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
		NodeWeight source_weight, target_weight, source_reachable_weight, target_reachable_weight;
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
			return !flow_algo.isSource(u) && !flow_algo.isTarget(u);
		}

		inline NodeWeight unclaimedNodeWeight() const {
			return hg.totalNodeWeight() - source_reachable_weight - target_reachable_weight;
		}

		inline NodeWeight notSettledNodeWeight() const {
			return hg.totalNodeWeight() - source_weight - target_weight;
		}

		void addToSourceSideCut(const Hyperedge e) {
			//Note: the current implementation of selecting piercing nodes relies on not inserting target-reachable nodes during most balanced cut mode
			if (!cuts.sourceSide.wasAdded(e)) {
				cuts.sourceSide.add(e);
				for (const Pin& px : hg.pinsOf(e)) {
					if (canBeSettled(px.pin) && !borderNodes.sourceSide->wasAdded(px.pin) && (!mostBalancedCutMode || !flow_algo.isTargetReachable(px.pin))) {
						borderNodes.sourceSide->add(px.pin, flow_algo.isTargetReachable(px.pin));
					}
				}
			}
		}

		void addToTargetSideCut(const Hyperedge e) {
			if (!cuts.targetSide.wasAdded(e)) {
				cuts.targetSide.add(e);
				for (const Pin& px : hg.pinsOf(e)) {
					if (canBeSettled(px.pin) && !borderNodes.targetSide->wasAdded(px.pin) && (!mostBalancedCutMode || !flow_algo.isSourceReachable(px.pin))) {
						borderNodes.targetSide->add(px.pin, flow_algo.isSourceReachable(px.pin));
					}
				}
			}
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

		void assimilate() {
			bool source_side = sideToGrow() == 0;
			if (source_side) {
				source_weight = source_reachable_weight;
				for (Node u : flow_algo.sourceReachableNodes()) {
					if (!flow_algo.isSource(u)) {
						if (mostBalancedCutMode && flow_algo.isHypernode(u)) {
							trackedMoves.emplace_back(u, 0, Move::Type::SettleNode);
						} else if (flow_algo.isInNode(u)) {
							Hyperedge e = flow_algo.inNodeToEdge(u);
							Node out_node = flow_algo.edgeToOutNode(e);
							if (!flow_algo.isSourceReachable(out_node)) {		// in node visited but not out node --> cut hyperedge
								addToSourceSideCut(e);
							}
						}
					}
					flow_algo.makeSource(u);
				}
			} else {
				target_weight = target_reachable_weight;
				for (Node u : flow_algo.sourceReachableNodes()) {
					if (!flow_algo.isTarget(u)) {
						if (mostBalancedCutMode && flow_algo.isHypernode(u)) {
							trackedMoves.emplace_back(u, 1, Move::Type::SettleNode);
						} else if (flow_algo.isOutNode(u)) {
							Hyperedge e = flow_algo.outNodeToEdge(u);
							Node in_node = flow_algo.edgeToInNode(e);
							if (!flow_algo.isTargetReachable(in_node)) {		// out node visited but not in node --> cut hyperedge
								addToTargetSideCut(e);
							}
						}
					}
					flow_algo.makeTarget(u);
				}
			}
		}

		void flipViewDirection() {
			viewDirection = 1 - viewDirection;
			hg.flipViewDirection();
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

		void reset() {		// TODO could consolidate with initialize
			viewDirection = 0;
			flowValue = 0;
			flow_algo.reset();
			h.fullReset();
			sourcePiercingNodes.clear();
			targetPiercingNodes.clear();
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
			source_weight = hg.nodeWeight(s);
			source_reachable_weight = source_weight;
			flow_algo.makeSource(s);

			targetPiercingNodes.emplace_back(t,false);
			target_weight = hg.nodeWeight(t);
			target_reachable_weight = target_weight;
			flow_algo.makeTarget(t);
		}

		int sideToGrow() const {
			const double imb_s = static_cast<double>(source_reachable_weight) / static_cast<double>(maxBlockWeight(currentViewDirection()));
			const double imb_t = static_cast<double>(target_reachable_weight) / static_cast<double>(maxBlockWeight(oppositeViewDirection()));
			return imb_s <= imb_t ? currentViewDirection() : oppositeViewDirection();
		}

		bool isBalanced() {
			assert(hasCut);
			assert(!partitionWrittenToNodeSet && "Cannot call isBalanced() once the partition has been written");

			const NodeWeight
					sw = source_reachable_weight,		//cannot be split
					tw = target_reachable_weight,		//cannot be split
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
			suw.imbalanceSourceBlock = block_imb(hg.totalNodeWeight() - target_reachable_weight, maxBlockWeight(currentViewDirection()));
			suw.imbalanceTargetBlock = block_imb(target_reachable_weight, maxBlockWeight(oppositeViewDirection()));

			SimulatedNodeAssignment tuw;
			tuw.imbalanceSourceBlock = block_imb(source_reachable_weight, maxBlockWeight(currentViewDirection()));
			tuw.imbalanceTargetBlock = block_imb(hg.totalNodeWeight() - source_reachable_weight, maxBlockWeight(oppositeViewDirection()));

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

			using result_t = std::pair<NodeWeight, NodeWeight>;
			result_t zero(0,0);
			result_t extra = tbb::parallel_reduce(
					tbb::blocked_range<Node>(Node(0), Node(hg.numNodes())), zero,
					[&](const auto& r, result_t sum) -> result_t {
						for (Node u = r.begin(); u < r.end(); ++u) {
							if (flow_algo.isSourceReachable(u) && !flow_algo.isSource(u)) {
								flow_algo.makeSource(u);
								sum.first += hg.nodeWeight(u);
							} else if (flow_algo.isTargetReachable(u) && !flow_algo.isTarget(u)) {
								flow_algo.settleTarget(u);
								sum.second += hg.nodeWeight(u);
							} else if (!flow_algo.isSourceReachable(u) && !flow_algo.isTargetReachable(u)) {
								if (r.assignUnclaimedToSource) {
									flow_algo.makeSource(u);
									sum.first += hg.nodeWeight(u);
								} else {
									flow_algo.makeTarget(u);
									sum.second += hg.nodeWeight(u);
								}
							}
						}
						return sum;
					}, [](const auto& l, const auto& r) -> result_t {
						return {l.first + r.first, l.second + r.second};
					});

			source_weight += extra.first;
			target_weight += extra.second;


			assert(source_weight + target_weight == hg.totalNodeWeight());
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
					assert(!flow_algo.isSourceReachable(m.node) && !flow_algo.isTargetReachable(m.node));
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
			   << " s=" << source_weight << "|" << source_reachable_weight
			   << " t=" << target_weight << "|" << target_reachable_weight;
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
					hasSource |= flow_algo.isSource(v);
					hasOther |= !flow_algo.isSource(v);
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
					hasSource |= flow_algo.isSource(p.pin);
					hasOther |= !flow_algo.isSource(p.pin);
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
				if (flow_algo.isSource(u)) {
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
							assert(!flow_algo.isTargetReachable(v));
							assert(flow_algo.isSourceReachable(v));
							if (!node_seen[v]) {
								node_seen.set(v);
								queue.push(v);
							}
						}
					}
				}
			}

			for (Node u : hg.nodeIDs()) {
				if (flow_algo.isTargetReachable(u))
					assert(!node_seen[u]);
				if (flow_algo.isSourceReachable(u))
					assert(node_seen[u]);
			}

			queue.clear();
			he_seen.reset();
			node_seen.reset();
			for (Node u : hg.nodeIDs()) {
				if (flow_algo.isTarget(u)) {
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
							assert(!flow_algo.isSourceReachable(v));
							//no assert(flow_algo.isTargetReachable(v)) since we removed the source-side cut
							if (!node_seen[v]) {
								node_seen.set(v);
								queue.push(v);
							}
						}
					}
				}
			}

			for (Node u : hg.nodeIDs()) {
				if (flow_algo.isTargetReachable(u))
					assert(node_seen[u]);
				if (flow_algo.isSourceReachable(u))
					assert(!node_seen[u]);
			}
#endif
		}

	};

}
