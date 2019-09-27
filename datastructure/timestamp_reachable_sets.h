#pragma once

#include <limits>
#include <vector>
#include <cassert>
#include <stdexcept>
#include "reachable_sets_base.h"

namespace whfc {
	template<typename Timestamp>
	class TimestampReachableNodes : public ReachableNodesBase {
	public:
		using Base = ReachableNodesBase;
		using Type = TimestampReachableNodes<Timestamp>;
		
		TimestampReachableNodes(const FlowHypergraph& hg) : Base(hg), timestamps(hg.numNodes(), unreachableTS) { }

		inline size_t capacity() const { return timestamps.size(); }
		inline bool isSource(const Node u) const { return timestamps[u] == sourceSettledTS; }
		inline bool isSourceReachable(const Node u) const { return isSource(u) || timestamps[u] == sourceReachableTS; }
		inline bool isTarget(const Node u) const { return timestamps[u] == targetSettledTS; }
		inline bool isTargetReachable(const Node u) const { return isTarget(u) || timestamps[u] == targetReachableTS; }
		inline void reach(const Node u) { assert(u < capacity()); assert(!isSourceReachable(u)); timestamps[u] = sourceReachableTS; Base::reach(u); }
		inline void settle(const Node u) { assert(!isSourceReachable(u)); timestamps[u] = sourceSettledTS; Base::settle(u); }
		
		
		inline void unreachSource(const Node u) { Assert(isSourceReachable(u) && !isTargetReachable(u)); timestamps[u] = unreachableTS; Base::unreachSource(u); }
		inline void unreachTarget(const Node u) { Assert(isTargetReachable(u) && !isSourceReachable(u)); timestamps[u] = unreachableTS; Base::unreachTarget(u); }
		
		
		void flipViewDirection() {
			std::swap(sourceSettledTS, targetSettledTS);
			std::swap(sourceReachableTS, targetReachableTS);
			Base::flipViewDirection();
		}

		void resetSourceReachableToSource() {
			if (generation == std::numeric_limits<Timestamp>::max()) {
				for (auto& ts : timestamps) {
					if (ts != sourceSettledTS && ts != targetSettledTS && ts != targetReachableTS)
						ts = unreachableTS;
				}
				generation = initialTS;
			}
			else
				generation++;
			sourceReachableTS = generation;
			Base::resetSourceReachableToSource();
		}

		void verifyDisjoint() const {
			//disjoint by default
		}

		void verifySettledIsSubsetOfReachable() const {
			//is subset by default
		}

	protected:
		static constexpr Timestamp initialTS = 3;
		static constexpr Timestamp unreachableTS = 0;
		Timestamp generation = initialTS;
		Timestamp sourceSettledTS = 1, targetSettledTS = 2;
		Timestamp sourceReachableTS = sourceSettledTS;
		Timestamp targetReachableTS = targetSettledTS;
		std::vector<Timestamp> timestamps;

	};

	template<typename Timestamp>
	class TimestampReachableHyperedges {
	public:
		using Type = TimestampReachableHyperedges<Timestamp>;

		TimestampReachableHyperedges(const FlowHypergraph& hg) : in(hg.numHyperedges(), unreachableTS), out(hg.numHyperedges(), unreachableTS) { }
		
		inline size_t capacity() const { return out.size(); }
		inline bool areAllPinsSources(const Hyperedge e) const { return out[e] == sourceSettledTS; }
		inline bool areAllPinsSourceReachable(const Hyperedge e) const { return areAllPinsSources(e) || out[e] == sourceReachableTS; }
		inline void settleAllPins(const Hyperedge e) { assert(!areAllPinsSources(e)); assert(in[e] != targetSettledTS); out[e] = sourceSettledTS; }
		inline void reachAllPins(const Hyperedge e) { assert(!areAllPinsSourceReachable(e)); assert(in[e] != targetSettledTS); out[e] = sourceReachableTS; }

		inline bool areFlowSendingPinsSources(const Hyperedge e) const { return in[e] == sourceSettledTS; }
		inline bool areFlowSendingPinsSourceReachable(const Hyperedge e) const { return areFlowSendingPinsSources(e) || in[e] == sourceReachableTS; }
		inline void settleFlowSendingPins(const Hyperedge e) { assert(!areFlowSendingPinsSources(e)); in[e] = sourceSettledTS; }
		inline void reachFlowSendingPins(const Hyperedge e) { assert(!areFlowSendingPinsSourceReachable(e)); in[e] = sourceReachableTS; }
		
		void resetSourceReachableToSource() {
			if (generation++ == std::numeric_limits<Timestamp>::max()) {
				for (auto& ts : out) if (ts != sourceSettledTS && ts != targetSettledTS && ts != targetReachableTS) ts = unreachableTS;
				for (auto& ts : in) if (ts != sourceSettledTS && ts != targetSettledTS && ts != targetReachableTS) ts = unreachableTS;
				generation = initialTS;
			}
			sourceReachableTS = generation;
		}

		void flipViewDirection() {
			std::swap(in, out);
			std::swap(sourceSettledTS, targetSettledTS);
			std::swap(sourceReachableTS, targetReachableTS);
		}

		void verifyDisjoint() const {
			//disjoint by default
		}

		void verifySettledIsSubsetOfReachable() const {
			//is subset by default
		}

	protected:
		static constexpr Timestamp initialTS = 3;
		static constexpr Timestamp unreachableTS = 0;
		Timestamp generation = initialTS;
		Timestamp sourceSettledTS = 1, targetSettledTS = 2;
		Timestamp sourceReachableTS = sourceSettledTS;
		Timestamp targetReachableTS = targetSettledTS;
		std::vector<Timestamp> in, out;
	};


}