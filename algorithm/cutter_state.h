#pragma once

#include "../datastructure/queue.h"
#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/isolated_nodes.h"
#include "../datastructure/bitset_reachable_sets.h"
#include "../util/math.h"

namespace whfc {
	
	struct SimulatedIsolatedNodesAssignment {
		bool assignUnclaimedToSource = true;
		bool assignTrackedIsolatedWeightToSource = true;
		NodeWeight trackedIsolatedWeight = NodeWeight::Invalid();
		NodeWeight blockWeightDiff = NodeWeight::Invalid();
		size_t numberOfTrackedMoves = 0;
		int direction;
	};
	
	struct Move {
		enum class Type : uint8_t { Isolate, Settle };
		Node node;
		int direction;
		Type t;
		Move(Node node, int dir, Type t) : node(node), direction(dir), t(t) { }
	};
	
	struct PiercingNode {
		Node node;
		bool isReachableFromOppositeSide;
		PiercingNode(const Node node, bool isReachableFromOppositeSide) : node(node), isReachableFromOppositeSide(isReachableFromOppositeSide) { }
	};
	
	template<typename FlowAlgorithm>
	class CutterState {
	public:
		using Pin = FlowHypergraph::Pin;
		static constexpr bool log = true;
		
		int viewDirection = 0;	//potential prettyfication: ViewDirection class
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
		bool prohibitIsolation = false;
		HyperedgeCut cut;
		NodeBorder borderNodes;
		NodeWeight maxBlockWeight;
		IsolatedNodes isolatedNodes;
		bool partitionWrittenToNodeSet = false;
		TimeReporter& timer;

		CutterState(FlowHypergraph& _hg, NodeWeight _maxBlockWeight, TimeReporter& timer) :
				hg(_hg),
				n(_hg),
				h(_hg),
				cut(_hg.numHyperedges()),
				borderNodes(_hg.numNodes()),
				maxBlockWeight(_maxBlockWeight),
				isolatedNodes(hg, _maxBlockWeight),
				timer(timer)
		{
			timer.registerCategory("Balance Check");
			timer.registerCategory("Output Balanced Partition");
		}
		
		inline bool isIsolated(const Node u) const {
			return !n.isSource(u) && !n.isTarget(u) && isolatedNodes.isCandidate(u);
		}
		
		inline bool canBeSettled(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && !isIsolated(u); }

		inline NodeWeight unclaimedNodeWeight() const { return hg.totalNodeWeight() - n.sourceReachableWeight - n.targetReachableWeight - isolatedNodes.weight; }
		inline NodeWeight notSettledNodeWeight() const { return hg.totalNodeWeight() - n.sourceWeight - n.targetWeight - isolatedNodes.weight; }
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

			if (mostBalancedCutMode)
				trackedMoves.emplace_back(u, currentViewDirection(), Move::Type::Settle);
			
			if (prohibitIsolation)
				return;
			
			for (const auto& he_inc : hg.hyperedgesOf(u)) {
				const Hyperedge e = he_inc.e;
				if (!hasSourcePin(e)) {
					cut.hasSettledSourcePins.set(e);
					if (hasTargetPin(e)) {	//e just became mixed
						for (const auto& px : hg.pinsOf(e)) {
							const Node p = px.pin;
							isolatedNodes.mixedIncidentHyperedges[p]++;
							if (isIsolated(p)) {
								isolatedNodes.add(p);
								if (n.isSourceReachable(p))
									n.unreachSource(p);
								if (n.isTargetReachable(p))
									n.unreachTarget(p);
								if (mostBalancedCutMode)
									trackedMoves.emplace_back(u, currentViewDirection(), Move::Type::Isolate);
							}
						}
					}
				}
			}
		}

		void flipViewDirection() {
			viewDirection = 1 - viewDirection;
			hg.flipViewDirection();
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
			n.resetSourceReachableToSource(augmentingPathAvailableFromPiercing);
			h.resetSourceReachableToSource(augmentingPathAvailableFromPiercing);
		}
		
		void reset() {
			viewDirection = 0;
			flowValue = 0;
			n.fullReset();
			h.fullReset();
			sourcePiercingNodes.clear(); targetPiercingNodes.clear();
			trackedMoves.clear();
			augmentingPathAvailableFromPiercing = true;
			hasCut = false;
			mostBalancedCutMode = false;
			prohibitIsolation = false;
			cut.reset(hg.numHyperedges());			//this requires that FlowHypergraph is reset before resetting the CutterState
			borderNodes.reset(hg.numNodes());
			isolatedNodes.reset();
			partitionWrittenToNodeSet = false;
		}
		
		void initialize(const Node s, const Node t) {
			Assert(sourcePiercingNodes.empty() && targetPiercingNodes.empty());
			sourcePiercingNodes.emplace_back(s,false);
			settleNode(s);
			targetPiercingNodes.emplace_back(t,false);
			flipViewDirection();
			settleNode(t);
			flipViewDirection();
		}
		
		void cleanUpBorder() {
			borderNodes.cleanUp([&](const Node& x) { return !canBeSettled(x); });
		}

		void cleanUpCut() {
			cut.cleanUp([&](const Hyperedge& e) { return h.areAllPinsSources(e); });
		}
		
		bool isBalanced() {
			Assert(hasCut);
			AssertMsg(!partitionWrittenToNodeSet, "Cannot call isBalanced() once the partition has been written");
			
			const NodeWeight
					sw = n.sourceReachableWeight,		//cannot be split
					tw = n.targetReachableWeight,		//cannot be split
					uw = unclaimedNodeWeight(),			//cannot be split (in current stages. if we integrate proper PCKP heuristics for MBMC this would change)
					iso = isolatedNodes.weight;			//can be split

			if (sw > maxBlockWeight || tw > maxBlockWeight)					//this is good at late and early stages
				return false;
			if (sw + uw > maxBlockWeight && tw + uw > maxBlockWeight)		//this is good at early stages
				return false;

			{	//quick checks to determine whether balance is possible without invoking SubsetSum, i.e. don't split the isolated nodes.
				bool balanced = false;
				balanced |= sw + uw + iso <= maxBlockWeight;
				balanced |= tw + uw + iso <= maxBlockWeight;
				balanced |= sw + uw <= maxBlockWeight && tw + iso <= maxBlockWeight;
				balanced |= tw + uw <= maxBlockWeight && sw + iso <= maxBlockWeight;
				if (balanced)
					return true;
			}
			
			timer.start("Balance Check");
			
			isolatedNodes.updateDPTable();

			const NodeWeight
					sRem = maxBlockWeight - sw,
					tRem = maxBlockWeight - tw,
					suw = sw + uw,
					tuw = tw + uw,
					suwRem = suw <= maxBlockWeight ? maxBlockWeight - suw : NodeWeight::Invalid(),
					tuwRem = tuw <= maxBlockWeight ? maxBlockWeight - tuw : NodeWeight::Invalid();
			
			bool balanced = false;
			
			//sides: (S + U, T) + <ISO> and (S, T + U) + <ISO>
			for (const IsolatedNodes::SummableRange& sr : isolatedNodes.getSumRanges()) {
				if (suwRem.isValid()) {
					//S+U not overloaded. Therefore, try (S + U, T) + <ISO>

					//allocate as much as possible to S+U, i.e. x = min(suwRem, sr.to), the rest, i.e. iso - x has to go to T
					balanced |= suwRem >= sr.from && tw + (iso - std::min(suwRem, sr.to)) <= maxBlockWeight;
					//analogously, allocate as much as possible to T
					balanced |= tRem >= sr.from && suw + (iso - std::min(tRem, sr.to)) <= maxBlockWeight;
				}

				if (tuwRem.isValid()) {
					//T+U not overloaded. Therefore, try (S, T + U) + <ISO>
					balanced |= tuwRem >= sr.from && sw + (iso - std::min(tuwRem, sr.to)) <= maxBlockWeight;
					balanced |= sRem >= sr.from && tuw + (iso - std::min(sRem, sr.to)) <= maxBlockWeight;
				}

				if (balanced)
					break;
			}

			timer.stop("Balance Check");
			return balanced;
		}

		// TODO consolidate isBalanced() and mostBalancedIsolatedNodesAssignment(). they do the same thing with slightly different purposes
		// still provide two functions that use the same core
		// e.g. mostBalancedNodesAssignment with a bound on the node weight diff as parameter
		// which is set to 0 when finding the best assignment and set to 2 * maxBlockWeight - hg.totalNodeWeight() when checking balance
		// stops searching when a solution lies below that bound

		
		/*
		 * Simulates settling all isolated and unclaimed nodes to achieve the most balanced partition possible with the current sides
		 * Write that assignment to ReachableNodes if write_partition = true
		 * returns the smallest possible block weight difference and the necessary assignment information
		 */
		SimulatedIsolatedNodesAssignment mostBalancedIsolatedNodesAssignment() {
			timer.start("Assign Isolated Nodes");
			
			isolatedNodes.updateDPTable();

			const NodeWeight
					sw = n.sourceReachableWeight,
					tw = n.targetReachableWeight,
					uw = unclaimedNodeWeight(),
					suw = sw + uw,
					tuw = tw + uw;

			
			bool assignUnclaimedToSource = true;
			bool assignTrackedIsolatedWeightToSource = true;
			NodeWeight trackedIsolatedWeight = NodeWeight::Invalid();
			NodeWeight blockWeightDiff = NodeWeight::Invalid();

			for (const IsolatedNodes::SummableRange& sr : isolatedNodes.getSumRanges()) {

				auto a = isolatedWeightAssignmentToFirst(suw, tw, sr);
				if (a.second < blockWeightDiff) {
					blockWeightDiff = a.second;
					assignUnclaimedToSource = true;
					assignTrackedIsolatedWeightToSource = true;
					trackedIsolatedWeight = a.first;
				}

				auto b = isolatedWeightAssignmentToFirst(tw, suw, sr);
				if (b.second < blockWeightDiff) {
					blockWeightDiff = b.second;
					assignUnclaimedToSource = true;
					assignTrackedIsolatedWeightToSource = false;
					trackedIsolatedWeight = b.first;
				}

				auto c = isolatedWeightAssignmentToFirst(sw, tuw, sr);
				if (c.second < blockWeightDiff) {
					blockWeightDiff = c.second;
					assignUnclaimedToSource = false;
					assignTrackedIsolatedWeightToSource = true;
					trackedIsolatedWeight = c.first;
				}

				auto d = isolatedWeightAssignmentToFirst(tuw, sw, sr);
				if (d.second < blockWeightDiff) {
					blockWeightDiff = d.second;
					assignUnclaimedToSource = false;
					assignTrackedIsolatedWeightToSource = false;
					trackedIsolatedWeight = d.first;
				}
				
			}
			
#ifndef NDEBUG
			const NodeWeight iso = isolatedNodes.weight;
			NodeWeight s = (assignUnclaimedToSource ? suw : sw) + (assignTrackedIsolatedWeightToSource ? trackedIsolatedWeight : iso - trackedIsolatedWeight);
			NodeWeight t = (assignUnclaimedToSource ? tw : tuw) + (assignTrackedIsolatedWeightToSource ? iso - trackedIsolatedWeight : trackedIsolatedWeight);
			const bool is_balanced = isBalanced();
			LOGGER << V(is_balanced);
			AssertMsg(s <= maxBlockWeight, "computed assignment violates max block weight on source side");
			AssertMsg(t <= maxBlockWeight, "computed assignment violates max block weight on target side");
			AssertMsg(isolatedNodes.isSummable(trackedIsolatedWeight), "isolated weight is not summable");
#endif
			
			timer.stop("Assign Isolated Nodes");
			return {assignUnclaimedToSource, assignTrackedIsolatedWeightToSource, trackedIsolatedWeight, blockWeightDiff, trackedMoves.size(), currentViewDirection()};
		}
		
		// takes the information from mostBalancedIsolatedNodesAssignment()
		// can be an old run, since the DP solution for trackedIsolatedWeight only contains nodes that were isolated during that run
		void writePartition(const SimulatedIsolatedNodesAssignment& r) {
			AssertMsg(isBalanced(), "Not balanced yet");
			AssertMsg(!partitionWrittenToNodeSet, "Partition was already written");
			
			if (currentViewDirection() != r.direction)
				flipViewDirection();
			
			auto isoSubset = isolatedNodes.extractSubset(r.trackedIsolatedWeight);
			for (const Node u : isoSubset) {
				Assert(!n.isSourceReachable(u) && !n.isTargetReachable(u) && isIsolated(u));
				if (r.assignTrackedIsolatedWeightToSource) {
					n.reach(u); n.settle(u);
				}
				else {
					n.reachTarget(u); n.settleTarget(u);
				}
			}
			
			for (const Node u : hg.nodeIDs()) {
				if (n.isSourceReachable(u) && !n.isSource(u))
					n.settle(u);
				
				if (n.isTargetReachable(u) && !n.isTarget(u))
					n.settleTarget(u);
				
				if (!n.isSourceReachable(u) && !n.isTargetReachable(u) && !isIsolated(u)) {
					if (r.assignUnclaimedToSource) {
						n.reach(u); n.settle(u);
					}
					else {
						n.reachTarget(u); n.settleTarget(u);
					}
				}
				
				if (isIsolated(u)) {
					if (r.assignTrackedIsolatedWeightToSource) {
						n.reach(u); n.settle(u);
					}
					else {
						n.reachTarget(u); n.settleTarget(u);
					}
				}
			}
			
			if (currentViewDirection() != 0)
				flipViewDirection();
			
			Assert(n.sourceWeight + n.targetWeight == hg.totalNodeWeight());
			partitionWrittenToNodeSet = true;
		}
		
		void writePartition() {
			writePartition(mostBalancedIsolatedNodesAssignment());
		}
		
		void revertMoves(SimulatedIsolatedNodesAssignment r) {
			timer.start("Revert Moves");
			while (trackedMoves.size() > r.numberOfTrackedMoves) {
				Move& m = trackedMoves.back();
				if (m.t == Move::Type::Isolate) {
					prohibitIsolation = true;	//in the next trials, we cannot isolate nodes, since we destroyed the datastructure
					isolatedNodes.mixedIncidentHyperedges[m.node] -= 1;
				}
				else {
					if (m.direction == currentViewDirection())
						n.unsettleSource(m.node);
					else
						n.unsettleTarget(m.node);
				}
				trackedMoves.pop_back();
			}
			timer.stop("Revert Moves");
		}
		
		std::string toString() {
			std::stringstream os;
			bool flipIt = currentViewDirection() != 0;
			if (flipIt)
				flipViewDirection();
			os << " cut= " << flowValue
			   << " s=" << n.sourceWeight << "|" << n.sourceReachableWeight
			   << " t=" << n.targetWeight << "|" << n.targetReachableWeight
			   << " iso=" << isolatedNodes.weight
			   << " u=" << unclaimedNodeWeight()
			   << " total=" << hg.totalNodeWeight();
			if (flipIt)
				flipViewDirection();
			return os.str();
		}


	private:

		//side of a gets x, side of b gets isolatedNodes.weight - x
		//result.first = x, result.second = block weight difference
		inline std::pair<NodeWeight, NodeWeight> isolatedWeightAssignmentToFirst(const NodeWeight a, NodeWeight b, const IsolatedNodes::SummableRange& sr) const {
			b += isolatedNodes.weight;
			const NodeWeight x = (a < b) ? std::max(std::min(NodeWeight((b-a)/2), sr.to), sr.from) : sr.from;
			return std::make_pair(x, Math::absdiff(a + x, b - x));
		}
		
	public:
		
		void verifyFlowConstraints() {
#ifndef NDEBUG
			Flow sourceExcess = 0, targetExcess = 0;
			for (Node u : hg.nodeIDs()) {
				Flow excess = 0;
				for (auto& he_inc : hg.hyperedgesOf(u))
					excess += hg.flowSent(he_inc);
				if (n.isSource(u))
					sourceExcess += excess;
				else if (n.isTarget(u))
					targetExcess += excess;
				else
					Assert(excess == 0);
			}
			Assert(sourceExcess >= 0 && targetExcess <= 0);
			Assert(hg.flowSent(sourceExcess) == hg.flowReceived(targetExcess));
			Assert(sourceExcess == flowValue);
			
			for (Hyperedge e : hg.hyperedgeIDs()) {
				Flow flow_in = 0, flow_out = 0;
				for (Pin& p : hg.pinsOf(e)) {
					Assert(std::abs(hg.flowSent(p)) <= hg.capacity(e));
					flow_in += hg.absoluteFlowSent(p);
					flow_out += hg.absoluteFlowReceived(hg.getInHe(p));
				}
				Assert(flow_in >= 0);
				Assert(flow_in == flow_out);
				Assert(flow_in == std::abs(hg.flow(e)));
				Assert(flow_in <= hg.capacity(e));
				
				for (Pin& p : hg.pinsSendingFlowInto(e))
					Assert(hg.flowSent(p) > 0);
				for (Pin& p : hg.pinsReceivingFlowFrom(e))
					Assert(hg.flowReceived(p) > 0);
				for (Pin& p : hg.pinsWithoutFlow(e))
					Assert(hg.flowSent(p) == 0);
			}
#endif
		}
		
		void verifyCutPostConditions() {
			Assert(hasCut);
			verifySetInvariants();
			verifyCutInducedByPartitionMatchesExtractedCutHyperedges();
			verifyExtractedCutHyperedgesActuallySplitHypergraph();
			Assert(flowValue == cut.weight(hg));
		}
		
		void verifySetInvariants() {
#ifndef NDEBUG
			n.verifyDisjoint();
			n.verifySettledIsSubsetOfReachable();
			h.verifyDisjoint();
			h.verifySettledIsSubsetOfReachable();
			for (Hyperedge e : hg.hyperedgeIDs()) {
				Assert(!h.areAllPinsSources(e) ||
					   std::all_of(hg.pinsOf(e).begin(), hg.pinsOf(e).end(),
								   [&](const Pin& p) { return n.isSource(p.pin) || isIsolated(p.pin); }));
				Assert(!h.areAllPinsSourceReachable(e) ||
					   std::all_of(hg.pinsOf(e).begin(), hg.pinsOf(e).end(),
								   [&](const Pin& p) { return n.isSourceReachable(p.pin); }));
				Assert(!h.areFlowSendingPinsSources(e) ||
					   std::all_of(hg.pinsSendingFlowInto(e).begin(), hg.pinsSendingFlowInto(e).end(),
								   [&](const Pin& p) { return n.isSource(p.pin) || isIsolated(p.pin); }));
				Assert(!h.areFlowSendingPinsSourceReachable(e) ||
					   std::all_of(hg.pinsSendingFlowInto(e).begin(), hg.pinsSendingFlowInto(e).end(),
								   [&](const Pin& p) { return n.isSourceReachable(p.pin); }));
			}
#endif
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
					Assert(h.areFlowSendingPinsSources(e));
				}
				
				if (hasSource && !hasOther)
					Assert(h.areAllPinsSources(e));
			}
			
			std::vector<Hyperedge> sorted_cut = cut.sourceSideBorder;
			std::sort(sorted_cut.begin(), sorted_cut.end());
			Assert(sorted_cut == cut_from_partition);
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
			
			for (Hyperedge e : cut.sourceSideBorder)
				he_seen.set(e);
			
			while (!queue.empty()) {
				Node u = queue.pop();
				for (auto& he_inc : hg.hyperedgesOf(u)) {
					Hyperedge e = he_inc.e;
					if (!he_seen[e]) {
						he_seen.set(e);
						for (auto& pin : hg.pinsOf(e)) {
							Node v = pin.pin;
							Assert(!n.isTargetReachable(v));
							Assert(n.isSourceReachable(v));
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
					Assert(!node_seen[u]);
				if (n.isSourceReachable(u))
					Assert(node_seen[u]);
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
			
			for (Hyperedge e : cut.sourceSideBorder)
				he_seen.set(e);
			
			while (!queue.empty()) {
				Node u = queue.pop();
				for (auto& he_inc : hg.hyperedgesOf(u)) {
					Hyperedge e = he_inc.e;
					if (!he_seen[e]) {
						he_seen.set(e);
						for (auto& pin : hg.pinsOf(e)) {
							Node v = pin.pin;
							Assert(!n.isSourceReachable(v));
							//no Assert(n.isTargetReachable(v)) since we removed the source-side cut
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
					Assert(node_seen[u]);
				if (n.isSourceReachable(u))
					Assert(!node_seen[u]);
			}
#endif
		}
		
	};



}
