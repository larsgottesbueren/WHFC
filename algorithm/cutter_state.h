#pragma once

#include "../definitions.h"
#include "../datastructure/border.h"
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/isolated_nodes.h"
#include "../datastructure/bitset_reachable_sets.h"
#include "../util/math.h"

namespace whfc {

	template<typename FlowAlgorithm>
	class CutterState {
		static constexpr bool debug = true;
	public:
		using Pin = FlowHypergraph::Pin;

		int viewDirection = 0;	//potential prettyfication: ViewDirection class
		FlowHypergraph& hg;
		Flow flowValue = 0;

		//using ReachableNodes = BitsetReachableNodes;
		//using ReachableHyperedges = BitsetReachableHyperedges;
		using ReachableNodes = typename FlowAlgorithm::ReachableNodes;
		using ReachableHyperedges = typename FlowAlgorithm::ReachableHyperedges;
		ReachableNodes n;
		ReachableHyperedges h;
		
		struct PiercingNode {
			Node node;
			bool isReachableFromOppositeSide;
			PiercingNode(const Node node, bool isReachableFromOppositeSide) : node(node), isReachableFromOppositeSide(isReachableFromOppositeSide) { }
		};
		std::vector<PiercingNode> sourcePiercingNodes, targetPiercingNodes;
		
		bool augmentingPathAvailableFromPiercing = true;
		bool hasCut = false;
		HyperedgeCut cut;
		NodeBorder borderNodes;
		NodeWeight maxBlockWeight;
		IsolatedNodes isolatedNodes;
		bool partitionWrittenToNodeSet = false;

		CutterState(FlowHypergraph& _hg, NodeWeight _maxBlockWeight) :
				hg(_hg),
				n(_hg),
				h(_hg),
				cut(_hg.numHyperedges()),
				borderNodes(_hg.numNodes()),
				maxBlockWeight(_maxBlockWeight),
				isolatedNodes(hg, _maxBlockWeight)
		{

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

		inline bool isIsolated(const Node u) const { return !n.isSource(u) && !n.isTarget(u) && isolatedNodes.isCandidate(u); }
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

			//isolatedNodes.settleNode(u);		only do this, if we ever decide to eject entries in IsolatedNodes, which would require resolving from scratch

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
			if (augmentingPathAvailableFromPiercing) {
				n.resetSourceReachableToSource();
				h.resetSourceReachableToSource();
			}
		}

		void cleanUpBorder() {
			borderNodes.cleanUp([&](const Node& x) { return !canBeSettled(x); });
		}

		void cleanUpCut() {
			cut.deleteNonCutHyperedges(h);
		}


		bool isBalanced() {
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

			isolatedNodes.updateDPTable();

			const NodeWeight
					sRem = maxBlockWeight - sw,
					tRem = maxBlockWeight - tw,
					suw = sw + uw,
					tuw = tw + uw,
					suwRem = suw <= maxBlockWeight ? maxBlockWeight - suw : NodeWeight::Invalid(),
					tuwRem = tuw <= maxBlockWeight ? maxBlockWeight - tuw : NodeWeight::Invalid();

			//sides: (S + U, T) + <ISO> and (S, T + U) + <ISO>
			for (const IsolatedNodes::SummableRange& sr : isolatedNodes.getSumRanges()) {
				bool balanced = false;
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
					return true;
			}

			return false;
		}


		/*
		 * Settles all nodes to their respective sides in the output partition
		 * maybe a different interface is better?
		 */
		void outputMostBalancedPartition() {
			AssertMsg(isolatedNodes.isDPTableUpToDate(), "DP Table not up to date");
			AssertMsg(isBalanced(), "Not balanced yet");
			AssertMsg(!partitionWrittenToNodeSet, "Partition was already written");

			const NodeWeight
					sw = n.sourceReachableWeight,
					tw = n.targetReachableWeight,
					uw = unclaimedNodeWeight(),
					iso = isolatedNodes.weight,
					suw = sw + uw,
					tuw = tw + uw;

			bool assignUnclaimedToSource = true;
			bool assignTrackedIsolatedWeightToSource = true;
			NodeWeight trackedIsolatedWeight = NodeWeight::Invalid();

			NodeWeight blockWeightDiff = NodeWeight::Invalid();

			for (const IsolatedNodes::SummableRange& sr : isolatedNodes.getSumRanges()) {

				{
					auto a = isolatedWeightAssignmentToFirst(suw, tw, sr);
					if (a.second < blockWeightDiff) {
						assignUnclaimedToSource = true;
						assignTrackedIsolatedWeightToSource = true;
						trackedIsolatedWeight = a.first;
					}
				}

				{
					auto b = isolatedWeightAssignmentToFirst(tw, suw, sr);
					if (b.second < blockWeightDiff) {
						assignUnclaimedToSource = true;
						assignTrackedIsolatedWeightToSource = false;
						trackedIsolatedWeight = b.first;
					}
				}

				{
					auto c = isolatedWeightAssignmentToFirst(sw, tuw, sr);
					if (c.second < blockWeightDiff) {
						assignUnclaimedToSource = false;
						assignTrackedIsolatedWeightToSource = true;
						trackedIsolatedWeight = c.first;
					}
				}

				{
					auto d = isolatedWeightAssignmentToFirst(tuw, sw, sr);
					if (d.second < blockWeightDiff) {
						assignUnclaimedToSource = false;
						assignTrackedIsolatedWeightToSource = false;
						trackedIsolatedWeight = d.first;
					}
				}

			}

			NodeWeight s = (assignUnclaimedToSource ? sw : suw) + (assignTrackedIsolatedWeightToSource ? trackedIsolatedWeight : iso - trackedIsolatedWeight);
			NodeWeight t = (assignUnclaimedToSource ? tuw : tw) + (assignTrackedIsolatedWeightToSource ? iso - trackedIsolatedWeight : trackedIsolatedWeight);
			AssertMsg(s <= maxBlockWeight, "computed assignment violates max block weight on source side");
			AssertMsg(t <= maxBlockWeight, "computed assignment violates max block weight on target side");
			AssertMsg(isolatedNodes.isSummable(trackedIsolatedWeight), "isolated weight is not summable");

			auto [sIso, tIso] = isolatedNodes.extractBipartition(trackedIsolatedWeight);		//commodity. could be handled otherwise.
			if (!assignTrackedIsolatedWeightToSource)
				std::swap(sIso, tIso);


			//TODO figure out desired output


			partitionWrittenToNodeSet = true;
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



		NodeWeight maxBlockWeightDiff() {
			if (!partitionWrittenToNodeSet)
				throw std::runtime_error("Partition wasn't written yet. Call outputMostBalancedPartition() first");
			return maxBlockWeight - std::min(n.sourceWeight,n.targetWeight);
		}

/*
		bool isInfeasible() {
			//call rarely. only intended for debugging purposes

			if (hg.totalNodeWeight() > 2 * maxBlockWeight)
				return true;
			LOG << "Starting expensive SubsetSum based feasibility check";
			IsolatedNodes subsetSumCopy = isolatedNodes;
			for (const Node u : hg.nodeIDs())
				if (canBeSettled(u))
					subsetSumCopy.add(u);
			subsetSumCopy.updateDPTable();

			//this is unfortunately as good as it gets. removing one element from a given subset sum instance is just as hard as re-solving from scratch
			std::swap(isolatedNodes, subsetSumCopy);
			const bool res = !isBalanced();
			std::swap(isolatedNodes, subsetSumCopy);
			return res;
	}
*/
	};



}