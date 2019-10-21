#pragma once

#include "reachable_sets_base.h"

namespace whfc {
	class DistanceRange {
	public:
		using DistanceT = uint32_t;
		DistanceT base, upper_bound;
		DistanceRange(const DistanceT d) : base(d), upper_bound(d) { }
		DistanceRange(const DistanceT base, const DistanceT upper_bound) : base(base), upper_bound(upper_bound) { }
		bool contains(const DistanceT d) const { return d >= base && d < upper_bound; }
		bool operator==(const DistanceRange& o) const { return o.base == base && o.upper_bound == upper_bound; }
	};
	
	class DistanceReachableHyperedges;
	
	class DistanceReachableNodes : public ReachableNodesBase {
	public:
		friend class DistanceReachableHyperedges;
		
		using Base = ReachableNodesBase;
		using DistanceT = DistanceRange::DistanceT;
		
		static constexpr bool log = true;
		
		DistanceReachableNodes(const FlowHypergraph& hg) : Base(hg), distance(hg.numNodes(), unreachableDistance), s(sourceSettledDistance), t(targetSettledDistance) {
			Assert(4 + hg.numNodes() * 2 < std::numeric_limits<DistanceT>::max());
		}
		
		inline size_t capacity() const { return distance.size(); }
		
		inline bool isSource(const Node u) const { return distance[u] == sourceSettledDistance; }
		inline bool isTarget(const Node u) const { return distance[u] == targetSettledDistance; }
		inline bool isSourceReachable__unsafe__(const Node u) const { return isSource(u) || distance[u] >= s.base; }	//saves a comparison in Dinic
		inline bool isSourceReachable(const Node u) const { return isSource(u) || s.contains(distance[u]); }
		inline bool isTargetReachable(const Node u) const { return isTarget(u) || t.contains(distance[u]); }
		inline void reach(const Node u) { Assert(u < hg.numNodes()); Assert(!isSourceReachable(u)); distance[u] = runningDistance; Base::reach(u); }
		inline void settle(const Node u) { Assert(!isSource(u)); distance[u] = sourceSettledDistance; Base::settle(u); }
		inline void reachTarget(const Node u) { Assert(!isSourceReachable(u) && !isTargetReachable(u)); distance[u] = t.base; Base::reachTarget(u); }
		inline void settleTarget(const Node u) { Assert(!isSourceReachable(u) && isTargetReachable(u)); distance[u] = targetSettledDistance; Base::settleTarget(u); }
		inline void unreachSource(const Node u) { Assert(isSourceReachable(u) && !isTargetReachable(u)); distance[u] = unreachableDistance; Base::unreachSource(u); }
		inline void unreachTarget(const Node u) { Assert(isTargetReachable(u) && !isSourceReachable(u)); distance[u] = unreachableDistance; Base::unreachTarget(u); }
		
		
		void fullReset() {
			distance.assign(hg.numNodes(), unreachableDistance);
			sourceSettledDistance = 1;
			targetSettledDistance = 2;
			runningDistance = resetBaseDistance;
			s = DistanceRange(sourceSettledDistance);
			t = DistanceRange(targetSettledDistance);
			Base::fullReset();
		}
		
		inline DistanceT hop() {
			return ++runningDistance;
		}

		void flipViewDirection() {
			std::swap(sourceSettledDistance, targetSettledDistance);
			std::swap(s,t);
			Base::flipViewDirection();
		}

		void resetSourceReachableToSource(bool augmenting_path_available) {
			if (!augmenting_path_available) {
#ifndef NDEBUG
				for (Node u : hg.nodeIDs())
					Assert(!s.contains(distance[u]));
#endif
			}
			if (!isBaseDistanceSafe()) {
				for (Node u(0); u < capacity(); ++u) {
					if (isSource(u) || isTarget(u))
						continue;
					if (isTargetReachable(u))
						distance[u] = resetBaseDistance + (distance[u] - t.base);
					else
						distance[u] = unreachableDistance;
				}
				t.upper_bound = resetBaseDistance + (t.upper_bound - t.base);
				t.base = resetBaseDistance;
				runningDistance = t.upper_bound;
			}
			Assert(isBaseDistanceSafe());
			s.base = runningDistance;
			s.upper_bound = std::numeric_limits<DistanceT>::max();
			Base::resetSourceReachableToSource();
		}

		void lockInSourceDistance() {
			verifyDistancesAreStale();
			s.upper_bound = runningDistance;
		}
		
		void setPiercingNodeDistance(const Node piercing_node, bool reset) {
			distance[piercing_node] = reset ? sourceSettledDistance : s.base;
		}
		
		void verifyDistancesAreStale() const {
			Assert(std::none_of(distance.begin(), distance.end(), [&](const DistanceT& dist) { return dist >= runningDistance; }));
		}
		
		bool isDistanceStale(const Node u) const {
			return distance[u] < s.base;
		}
		
		void verifyDisjoint() const {
			// disjoint by default
			// but check whether distance ranges are disjoint
			Assert(s.base <= s.upper_bound);
			Assert(t.base <= t.upper_bound);
			Assert(s.upper_bound <= t.base || t.upper_bound <= s.base);
		}
		
		void verifySettledIsSubsetOfReachable() const {
			// subset by default
		}
		
		std::vector<DistanceT> distance;
		
		DistanceT sourceBaseDistance() const {
			return s.base;
		}
		
		inline bool isBaseDistanceSafe() const { return hg.numNodes() + static_cast<size_t>(runningDistance) < static_cast<size_t>(std::numeric_limits<DistanceT>::max()) ; }
		static constexpr DistanceT unreachableDistance = 0;
		DistanceT sourceSettledDistance = 1, targetSettledDistance = 2;
		static constexpr DistanceT resetBaseDistance = 3;
		DistanceT runningDistance = resetBaseDistance;
		DistanceRange s, t;
		
	};

	class DistanceReachableHyperedges {
	public:
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		DistanceReachableHyperedges(const FlowHypergraph& hg) : inDistance(hg.numHyperedges(), unreachableDistance), outDistance(hg.numHyperedges(), unreachableDistance),
																hg(hg), s(sourceSettledDistance), t(targetSettledDistance) { }
		
		inline size_t capacity() const { return outDistance.size(); }
		inline bool areAllPinsSources(const Hyperedge e) const { return outDistance[e] == sourceSettledDistance; }
		inline bool areAllPinsSourceReachable__unsafe__(const Hyperedge e) const { return areAllPinsSources(e) || outDistance[e] >= s.base; }
		inline bool areAllPinsSourceReachable(const Hyperedge e) const { return areAllPinsSources(e) || s.contains(outDistance[e]); }
		inline void settleAllPins(const Hyperedge e) { Assert(!areAllPinsSources(e)); outDistance[e] = sourceSettledDistance; }
		inline void reachAllPins(const Hyperedge e) { Assert(!areAllPinsSourceReachable(e)); outDistance[e] = runningDistance; }

		inline bool areFlowSendingPinsSources(const Hyperedge e) const { return inDistance[e] == sourceSettledDistance; }
		inline bool areFlowSendingPinsSourceReachable__unsafe__(const Hyperedge e) const { return areFlowSendingPinsSources(e) || inDistance[e] >= s.base; }
		inline bool areFlowSendingPinsSourceReachable(const Hyperedge e) const { return areFlowSendingPinsSources(e) || s.contains(inDistance[e]); }
		inline void settleFlowSendingPins(const Hyperedge e) { Assert(!areFlowSendingPinsSources(e)); inDistance[e] = sourceSettledDistance; }
		inline void reachFlowSendingPins(const Hyperedge e) { Assert(!areFlowSendingPinsSourceReachable(e)); inDistance[e] = runningDistance; }

		inline DistanceT hop() { return ++runningDistance; }

		void flipViewDirection() {
			std::swap(sourceSettledDistance, targetSettledDistance);
			std::swap(s, t);
			std::swap(inDistance, outDistance);
		}

		void fullReset() {
			inDistance.assign(hg.numHyperedges(), unreachableDistance);
			outDistance.assign(hg.numHyperedges(), unreachableDistance);
			sourceSettledDistance = 1; targetSettledDistance = 2;
			runningDistance = resetBaseDistance;
			s = DistanceRange(sourceSettledDistance);
			t = DistanceRange(targetSettledDistance);
		}
		
		void resetSourceReachableToSource(bool augmenting_path_available) {
			if (!augmenting_path_available) {
#ifndef NDEBUG
				for (Hyperedge e : hg.hyperedgeIDs()) {
					Assert(!s.contains(inDistance[e]));
					Assert(!s.contains(outDistance[e]));
				}
#endif
			}
			if (!isBaseDistanceSafe()) {
				auto reset = [&](std::vector<DistanceT>& d) {
					for (Hyperedge e : hg.hyperedgeIDs()) {
						if (d[e] == sourceSettledDistance || d[e] == targetSettledDistance)
							continue;
						if (t.contains(d[e]))
							d[e] = resetBaseDistance + (d[e] - t.base);
						else
							d[e] = unreachableDistance;
					}
				};
				reset(inDistance);
				reset(outDistance);

				t.upper_bound = resetBaseDistance + (t.upper_bound - t.base);
				t.base = resetBaseDistance;
				runningDistance = t.upper_bound;
			}
			assert(isBaseDistanceSafe());
			s.base = runningDistance;
			s.upper_bound = std::numeric_limits<DistanceT>::max();
		}

		void lockInSourceDistance() {
			verifyDistancesAreStale();
			s.upper_bound = runningDistance;
		}
		
		void verifyDisjoint() const { /*disjoint by default*/ }
		void verifySettledIsSubsetOfReachable() const { /*is subset by default*/ }
		void verifyDistancesAreStale() const {
			Assert(std::none_of(outDistance.begin(), outDistance.end(), [&](const DistanceT& dist) { return dist >= runningDistance; }));
			Assert(std::none_of(inDistance.begin(), inDistance.end(), [&](const DistanceT& dist) { return dist >= runningDistance; }));
		}
		
		void compareDistances(DistanceReachableNodes& n) {
			unused(n);
			Assert(n.sourceSettledDistance == sourceSettledDistance);//same direction?
			Assert(n.runningDistance == runningDistance);
			Assert(n.s == s);
			Assert(n.t == t);
			Assert(n.isBaseDistanceSafe() == isBaseDistanceSafe());
		}
		
		std::vector<DistanceT> inDistance, outDistance;
	//protected:
		const FlowHypergraph& hg;
		static constexpr DistanceT unreachableDistance = 0;
		DistanceT sourceSettledDistance = 1, targetSettledDistance = 2;
		static constexpr DistanceT resetBaseDistance = 3;
		DistanceT runningDistance = resetBaseDistance;
		DistanceRange s, t;
		
		inline bool isBaseDistanceSafe() const { return hg.numNodes() + static_cast<size_t>(runningDistance) < static_cast<size_t>(std::numeric_limits<DistanceT>::max()) ; }
		
	};

}