#pragma once

#include "../datastructure/queue.h"
#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/node_border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../util/math.h"
#include "../util/random.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>

namespace whfc {
	template<typename T>
	using vec = std::vector<T, tbb::scalable_allocator<T> >;

	struct SimulatedNodeAssignment {
		bool assignUnclaimedToSource = true;
		double imbalanceSourceBlock = std::numeric_limits<double>::max(), imbalanceTargetBlock = std::numeric_limits<double>::max();
		size_t numberOfTrackedMoves = 0;

		double imbalance() const {
			return std::max(imbalanceSourceBlock, imbalanceTargetBlock);
		}

		bool isPerfectlyBalanced() const {
			return std::abs(imbalanceSourceBlock - imbalanceTargetBlock) < 1e-9;
		}
	};

	struct Move {
		Node node;
		int direction;
		Move(Node node, int dir) : node(node), direction(dir) { }
	};

	struct NonDynamicCutterState {
		vec<Node> sourcePiercingNodes, targetPiercingNodes;
	};

	template<typename FlowAlgorithm>
	class CutterState {
	public:
		static constexpr bool log = false;

		using Pin = FlowHypergraph::Pin;

		FlowAlgorithm flow_algo;
		int side_to_pierce = 0;
		FlowHypergraph& hg;

		NodeWeight source_weight, target_weight, source_reachable_weight, target_reachable_weight;
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
					if (canBeSettled(px.pin) && !borderNodes.sourceSide.wasAdded(px.pin) && (!mostBalancedCutMode || !flow_algo.isTargetReachable(px.pin))) {
						borderNodes.sourceSide.add(px.pin, flow_algo.isTargetReachable(px.pin));
					}
				}
			}
		}

		void addToTargetSideCut(const Hyperedge e) {
			if (!cuts.targetSide.wasAdded(e)) {
				cuts.targetSide.add(e);
				for (const Pin& px : hg.pinsOf(e)) {
					if (canBeSettled(px.pin) && !borderNodes.targetSide.wasAdded(px.pin) && (!mostBalancedCutMode || !flow_algo.isSourceReachable(px.pin))) {
						borderNodes.targetSide.add(px.pin, flow_algo.isSourceReachable(px.pin));
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

		bool reachableFromSideNotToPierce(const Node u) const {
			return side_to_pierce == 0 ? flow_algo.isTargetReachable(u) : flow_algo.isSourceReachable(u);
		}

		void setPiercingNode(const Node piercingNode) {
			augmentingPathAvailableFromPiercing = reachableFromSideNotToPierce(piercingNode);
			if (side_to_pierce == 0) {
				flow_algo.source_piercing_nodes.clear();
				flow_algo.source_piercing_nodes.emplace_back(piercingNode);
				flow_algo.makeSource(piercingNode);
				source_weight += hg.nodeWeight(piercingNode);
			} else {
				flow_algo.target_piercing_nodes.clear();
				flow_algo.target_piercing_nodes.emplace_back(piercingNode);
				flow_algo.makeTarget(piercingNode);
				target_weight += hg.nodeWeight(piercingNode);
			}
			hasCut = false;
		}

		void computeReachableWeights() {
			if (augmentingPathAvailableFromPiercing) {
				tbb::parallel_invoke([&]{ computeSourceReachableWeight(); }, [&]{ computeTargetReachableWeight(); });
			} else {
				// no flow increased --> one side didn't change
				if (side_to_pierce == 0) {
					computeSourceReachableWeight();
				} else {
					computeTargetReachableWeight();
				}
			}
		}

		void computeSourceReachableWeight() {
			auto sr = flow_algo.sourceReachableNodes();
			source_reachable_weight = source_weight + tbb::parallel_reduce(
					tbb::blocked_range<size_t>(0, sr.size()), 0, [&](const auto& r, NodeWeight sum) -> NodeWeight {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					Node u = sr[i];	// next_active container is for source side. active for target side
					if (flow_algo.isHypernode(u) && !flow_algo.isSource(u)) {
						sum += hg.nodeWeight(u);
					}
				}
				return sum;
			}, std::plus<>());
		}

		void computeTargetReachableWeight() {
			auto tr = flow_algo.targetReachableNodes();
			target_reachable_weight = target_weight + tbb::parallel_reduce(
					tbb::blocked_range<size_t>(0, tr.size()), 0, [&](const auto& r, NodeWeight sum) -> NodeWeight {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					Node u = tr[i];
					if (flow_algo.isHypernode(u) && !flow_algo.isTarget(u)) {
						sum += hg.nodeWeight(u);
					}
				}
				return sum;
			}, std::plus<>());
		}

		void assimilateSourceSide() {
			source_weight = source_reachable_weight;
			for (Node u : flow_algo.sourceReachableNodes()) {
				if (!flow_algo.isSource(u)) {
					if (mostBalancedCutMode) {
						trackedMoves.emplace_back(u, 0);
					}
					if (flow_algo.isInNode(u)) {
						Hyperedge e = flow_algo.inNodeToEdge(u);
						Node out_node = flow_algo.edgeToOutNode(e);
						if (!flow_algo.isSourceReachable(out_node)) {		// in node visited but not out node --> cut hyperedge
							addToSourceSideCut(e);
						}
					}
					flow_algo.makeSource(u);
				}
			}
		}

		void assimilateTargetSide() {
			target_weight = target_reachable_weight;
			for (Node u : flow_algo.targetReachableNodes()) {
				if (!flow_algo.isTarget(u)) {
					if (mostBalancedCutMode) {
						trackedMoves.emplace_back(u, 1);
					}
					if (flow_algo.isOutNode(u)) {
						Hyperedge e = flow_algo.outNodeToEdge(u);
						Node in_node = flow_algo.edgeToInNode(e);
						if (!flow_algo.isTargetReachable(in_node)) {		// out node visited but not in node --> cut hyperedge
							addToTargetSideCut(e);
						}
					}
					flow_algo.makeTarget(u);
				}
			}
		}

		void assimilate() {
			computeReachableWeights();

			side_to_pierce = sideToGrow();

			if (side_to_pierce == 0 /* source side */) {
				assimilateSourceSide();
			} else {
				assimilateTargetSide();
			}
		}


		void reset() {		// TODO could consolidate with initialize
			flow_algo.reset();
			trackedMoves.clear();
			augmentingPathAvailableFromPiercing = true;
			hasCut = false;
			mostBalancedCutMode = false;
			cuts.reset(hg.numHyperedges());			//this requires that FlowHypergraph is reset before resetting the CutterState
			borderNodes.reset(hg.numNodes());
			partitionWrittenToNodeSet = false;
		}

		void initialize(Node s, Node t) {
			if (hg.nodeWeight(s) > maxBlockWeight(0) || hg.nodeWeight(t) > maxBlockWeight(1)) {
				throw std::runtime_error("Terminal weight already exceeds max block weight at initialization. Consider setting max block weights per side via hfc.cs.setMaxBlockWeight(  side  )");
			}

			flow_algo.initialize(s, t);

			source_weight = hg.nodeWeight(s);
			source_reachable_weight = source_weight;
			target_weight = hg.nodeWeight(t);
			target_reachable_weight = target_weight;
		}

		int sideToGrow() const {
			const double imb_s = static_cast<double>(source_reachable_weight) / static_cast<double>(maxBlockWeight(0));
			const double imb_t = static_cast<double>(target_reachable_weight) / static_cast<double>(maxBlockWeight(1));
			return imb_s <= imb_t ? 0 : 1;
		}

		bool isBalanced() {
			assert(hasCut);
			assert(!partitionWrittenToNodeSet && "Cannot call isBalanced() once the partition has been written");

			const NodeWeight
					sw = source_reachable_weight,		//cannot be split
					tw = target_reachable_weight,		//cannot be split
					uw = unclaimedNodeWeight();			//cannot be split (in current stages. if we integrate proper PCKP heuristics for MBMC this would change)

			const NodeWeight
					s_mbw = maxBlockWeight(0),
					t_mbw = maxBlockWeight(1);

			if (sw > s_mbw || tw > t_mbw)					//this is good at late and early stages
				return false;
			if (sw + uw > s_mbw && tw + uw > t_mbw)			//this is good at early stages
				return false;

			bool balanced = false;
			balanced |= sw + uw <= s_mbw && tw <= t_mbw;
			balanced |= tw + uw <= t_mbw && sw <= s_mbw;
			return balanced;
		}

		bool rejectPiercingIfAugmenting() const {
			return mostBalancedCutMode || flow_algo.flow_value == flow_algo.upper_flow_bound;
		}

		bool betterBalanceImpossible() const {
			if (unclaimedNodeWeight() == 0) {
				return false;
			}
			if (sideToGrow() == 0) {
				const double imb_S_U = static_cast<double>(hg.totalNodeWeight() - target_reachable_weight) / static_cast<double>(maxBlockWeight(0));
				const double imb_T = static_cast<double>(target_reachable_weight) / static_cast<double>(maxBlockWeight(1));
				return imb_S_U <= imb_T;
			} else {
				const double imb_S = static_cast<double>(source_reachable_weight) / static_cast<double>(maxBlockWeight(0));
				const double imb_T_U = static_cast<double>(hg.totalNodeWeight() - source_reachable_weight) / static_cast<double>(maxBlockWeight(1));
				return imb_T_U <= imb_S;
			}
		}

		NonDynamicCutterState enterMostBalancedCutMode() {
			assert(!mostBalancedCutMode);
			assert(trackedMoves.empty());
			assert(hasCut);
			mostBalancedCutMode = true;	// activates move tracking
			borderNodes.enterMostBalancedCutMode();
			cuts.enterMostBalancedCutMode();
			return { flow_algo.source_piercing_nodes, flow_algo.target_piercing_nodes };
		}

		void resetToFirstBalancedState(NonDynamicCutterState& nds) {
			flow_algo.source_piercing_nodes = nds.sourcePiercingNodes;
			flow_algo.target_piercing_nodes = nds.targetPiercingNodes;
			revertMoves(0);
			borderNodes.resetForMostBalancedCut();
			cuts.resetForMostBalancedCut();
		}

		SimulatedNodeAssignment mostBalancedAssignment() {
			auto block_imb = [&](NodeWeight a, NodeWeight max_a) {
				return (static_cast<double>(a) / static_cast<double>(max_a)) - 1.0;
			};
			SimulatedNodeAssignment suw;
			suw.imbalanceSourceBlock = block_imb(hg.totalNodeWeight() - target_reachable_weight, maxBlockWeight(0));
			suw.imbalanceTargetBlock = block_imb(target_reachable_weight, maxBlockWeight(1));

			SimulatedNodeAssignment tuw;
			tuw.imbalanceSourceBlock = block_imb(source_reachable_weight, maxBlockWeight(0));
			tuw.imbalanceTargetBlock = block_imb(hg.totalNodeWeight() - source_reachable_weight, maxBlockWeight(1));

			SimulatedNodeAssignment sol = suw.imbalance() < tuw.imbalance() ? suw : tuw;

			sol.numberOfTrackedMoves = trackedMoves.size();
			return sol;
		}

		// takes the information from mostBalancedIsolatedNodesAssignment()
		// can be an old run, since the DP solution for trackedIsolatedWeight only contains nodes that were isolated during that run
		void writePartition(const SimulatedNodeAssignment& assignment) {
			assert(!partitionWrittenToNodeSet);
			assert(isBalanced());

			using result_t = std::pair<NodeWeight, NodeWeight>;
			result_t zero(0,0);
			result_t extra = tbb::parallel_reduce(
					tbb::blocked_range<size_t>(0, hg.numNodes()), zero,
					[&](const auto& range, result_t sum) -> result_t {
						for (Node u(range.begin()); u < range.end(); ++u) {
							if (flow_algo.isSourceReachable(u) && !flow_algo.isSource(u)) {
								flow_algo.makeSource(u);
								sum.first += hg.nodeWeight(u);
							} else if (flow_algo.isTargetReachable(u) && !flow_algo.isTarget(u)) {
								flow_algo.makeTarget(u);
								sum.second += hg.nodeWeight(u);
							} else if (!flow_algo.isSourceReachable(u) && !flow_algo.isTargetReachable(u)) {
								if (assignment.assignUnclaimedToSource) {
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
			// only in most balanced cut mode --> no need for parallelism
			while (trackedMoves.size() > numberOfTrackedMoves) {
				Move& m = trackedMoves.back();
				flow_algo.unreach(m.node);
				if (flow_algo.isHypernode(m.node)) {
					if (m.direction == 0) source_weight -= hg.nodeWeight(m.node);
					else target_weight -= hg.nodeWeight(m.node);
				}
				trackedMoves.pop_back();
			}
			source_reachable_weight = source_weight;
			target_reachable_weight = target_weight;
		}

		void applyMoves(const std::vector<Move>& moves) {
			for (const Move& m : moves) {
				if (m.direction == 0) {
					flow_algo.makeSource(m.node);
					if (flow_algo.isHypernode(m.node)) source_weight += hg.nodeWeight(m.node);
				}
				else {
					flow_algo.makeTarget(m.node);
					if (flow_algo.isHypernode(m.node)) target_weight += hg.nodeWeight(m.node);
				}
			}
			source_reachable_weight = source_weight;
			target_reachable_weight = target_weight;
		}

		std::string toString(bool skip_iso_and_unclaimed = false) {
			std::stringstream os;
			os << " cut= " << flow_algo.flow_value
			   << " s=" << source_weight << "|" << source_reachable_weight
			   << " t=" << target_weight << "|" << target_reachable_weight;
			if (!skip_iso_and_unclaimed)
			   os << " u=" << unclaimedNodeWeight();
			os << " mbw=[" << maxBlockWeight(0) << " " << maxBlockWeight(1) << "]"
			   << " total=" << hg.totalNodeWeight()
			   ;
			return os.str();
		}



		void verifyCutPostConditions() {
			assert(hasCut);

#ifndef NDEBUG
			Flow expected_flow = 0;
			if (side_to_pierce == 0) {
				cuts.sourceSide.cleanUp([&](const Hyperedge& e) { return flow_algo.isSource(flow_algo.edgeToOutNode(e)); });
				for (const Hyperedge& e : cuts.sourceSide.entries()) {
					assert(flow_algo.flow[flow_algo.bridgeEdgeIndex(e)] == hg.capacity(e));
					expected_flow += hg.capacity(e);
				}
			} else {
				cuts.targetSide.cleanUp([&](const Hyperedge& e) { return flow_algo.isTarget(flow_algo.edgeToInNode(e)); });
				for (const Hyperedge& e : cuts.targetSide.entries()) {
					assert(flow_algo.flow[flow_algo.bridgeEdgeIndex(e)] == hg.capacity(e));
					expected_flow += hg.capacity(e);
				}
			}
			assert(flow_algo.flow_value == expected_flow);
#endif
			verifyExtractedCutHyperedgesActuallySplitHypergraph();
			verifyCutInducedByPartitionMatchesExtractedCutHyperedges();
		}

		void verifyCutInducedByPartitionMatchesExtractedCutHyperedges() {
#ifndef NDEBUG
			std::vector<Hyperedge> cut_from_partition;
			if (side_to_pierce == 0) {
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
						assert(flow_algo.isSource(flow_algo.edgeToInNode(e)));
					}

					if (hasSource && !hasOther) {
						assert(flow_algo.isSource(flow_algo.edgeToOutNode(e)));
					}
				}
				auto sorted_cut = cuts.sourceSide.copy();
				std::sort(sorted_cut.begin(), sorted_cut.end());
				assert(sorted_cut == cut_from_partition);
			} else {
				for (Hyperedge e : hg.hyperedgeIDs()) {
					bool hasTarget = false;
					bool hasOther = false;
					for (Pin& p : hg.pinsOf(e)) {
						Node v = p.pin;
						hasTarget |= flow_algo.isTarget(v);
						hasOther |= !flow_algo.isTarget(v);
					}
					if (hasTarget && hasOther) {
						cut_from_partition.push_back(e);
						assert(flow_algo.isTarget(flow_algo.edgeToOutNode(e)));
					}

					if (hasTarget && !hasOther) {
						assert(flow_algo.isTarget(flow_algo.edgeToInNode(e)));
					}
				}
				auto sorted_cut = cuts.targetSide.copy();
				std::sort(sorted_cut.begin(), sorted_cut.end());
				assert(sorted_cut == cut_from_partition);
			}
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
			assert(flow_algo.flow_value == cut_weight);
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

			if (side_to_pierce == 0) {
				for (Hyperedge e : cuts.sourceSide.entries()) he_seen.set(e);
			} else {
				for (Hyperedge e : cuts.targetSide.entries()) he_seen.set(e);
			}

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

			if (side_to_pierce == 0) {
				for (Hyperedge e : cuts.sourceSide.entries()) he_seen.set(e);
			} else {
				for (Hyperedge e : cuts.targetSide.entries()) he_seen.set(e);
			}

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
