#pragma once

#include "bitset_reachable_sets.h"
#include "timestamp_reachable_sets.h"

namespace whfc {

	class ReachableNodesChecker : public ReachableNodesBase {
	public:
		using Base = ReachableNodesBase;
		explicit ReachableNodesChecker(const FlowHypergraph& hg) : Base(hg), bits(hg), timestamps(hg) { }
		
		inline size_t capacity() const {
			Assert(bits.capacity() == timestamps.capacity());
			return bits.capacity();
		}
		
		inline bool isSource(const Node u) const {
			Assert(bits.isSource(u) == timestamps.isSource(u));
			return bits.isSource(u);
		}
		
		inline bool isSourceReachable(const Node u) const {
			Assert(bits.isSourceReachable(u) == timestamps.isSourceReachable(u));
			return bits.isSourceReachable(u);
		}
		inline bool isTarget(const Node u) const {
			Assert(bits.isTarget(u) == timestamps.isTarget(u));
			return bits.isTarget(u);
		}
		inline bool isTargetReachable(const Node u) const {
			Assert(bits.isTargetReachable(u) == timestamps.isTargetReachable(u));
			return bits.isTargetReachable(u);
		}
		inline void reach(const Node u) {
			Base::reach(u);
			bits.reach(u);
			timestamps.reach(u);
		}
		inline void settle(const Node u) {
			Base::settle(u);
			bits.settle(u);
			timestamps.settle(u);
		}
		
		inline void unreachSource(const Node u) {
			Base::unreachSource(u);
			bits.unreachSource(u);
			timestamps.unreachSource(u);
		}
		
		inline void unreachTarget(const Node u) {
			Base::unreachTarget(u);
			bits.unreachTarget(u);
			timestamps.unreachTarget(u);
		}
		
		
		void flipViewDirection() {
			Base::flipViewDirection();
			bits.flipViewDirection();
			timestamps.flipViewDirection();
		}
		
		void resetSourceReachableToSource() {
			Base::resetSourceReachableToSource();
			bits.resetSourceReachableToSource();
			timestamps.resetSourceReachableToSource();
		}
		
		void verifyDisjoint() const {
			bits.verifyDisjoint();
			timestamps.verifyDisjoint();
		}
		
		void verifySettledIsSubsetOfReachable() const {
			bits.verifySettledIsSubsetOfReachable();
			timestamps.verifySettledIsSubsetOfReachable();
		}
		
	private:
		BitsetReachableNodes bits;
		TimestampReachableNodes<uint8_t> timestamps;
		
	};
	
}