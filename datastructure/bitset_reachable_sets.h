#pragma once

#include "reachable_sets_base.h"
#include "bitvector.h"

namespace whfc {

		//TODO Find way to wrap other ReachableSets with a BitsetReachableSet and compare


		class BitsetReachableNodes : public ReachableNodesBase {
		public:
			using Base = ReachableNodesBase;
			using Type = BitsetReachableNodes;

			explicit BitsetReachableNodes(const FlowHypergraph& hg) : Base(hg),
															 S(hg.numNodes()),
															 SR(hg.numNodes()),
															 T(hg.numNodes()),
															 TR(hg.numNodes())
			{

			}

			inline size_t capacity() const { return S.size(); }
			inline bool isSource(const Node u) const { return S[u]; }
			inline bool isSourceReachable(const Node u) const { return SR[u]; }
			inline bool isTarget(const Node u) const { return T[u]; }
			inline bool isTargetReachable(const Node u) const { return TR[u]; }
			inline void reach(const Node u) { SR.set(u); Base::reach(u); }
			inline void settle(const Node u) { S.set(u); Base::settle(u); }

			inline void unreachSource(const Node u) { Assert(!isSourceReachable(u)); SR.reset(u); Base::unreachSource(u); }
			inline void unreachTarget(const Node u) { Assert(!isTargetReachable(u)); SR.reset(u); Base::unreachTarget(u); }


			void flipViewDirection() {
				std::swap(S, T);
				std::swap(SR, TR);
				Base::flipViewDirection();
			}

			void resetSourceReachableToSource() {
				SR = S;
				Base::resetSourceReachableToSource();
			}

			void verifyDisjoint() const {
				Assert((SR & TR).none());
				Assert((S & T).none());
			}

			void verifySettledIsSubsetOfReachable() const {
				Assert(S.is_subset_of(SR));
				Assert(T.is_subset_of(TR));
			}

		protected:
			BitVector S, SR, T, TR;
		};

		class BitsetReachableHyperedges {
		public:
			using Type = BitsetReachableHyperedges;

			explicit BitsetReachableHyperedges(const size_t nHE) :
					IN_SETTLED_S(nHE),
					OUT_SETTLED_S(nHE),
					IN_REACHED_S(nHE),
					OUT_REACHED_S(nHE),
					IN_SETTLED_T(nHE),
					OUT_SETTLED_T(nHE),
					IN_REACHED_T(nHE),
					OUT_REACHED_T(nHE)
			{

			}

			inline size_t capacity() const { return IN_SETTLED_S.size(); }
			inline bool areAllPinsSources(const Hyperedge e) const { return OUT_SETTLED_S[e]; }
			inline bool areAllPinsSourceReachable(const Hyperedge e) const { return OUT_REACHED_S[e]; }
			inline void settleAllPins(const Hyperedge e) { Assert(!areAllPinsSources(e)); OUT_SETTLED_S.set(e); IN_SETTLED_S.set(e);}
			inline void reachAllPins(const Hyperedge e) { Assert(!areAllPinsSourceReachable(e)); OUT_REACHED_S.set(e); IN_REACHED_S.set(e);}

			inline bool areFlowSendingPinsSources(const Hyperedge e) const { return IN_SETTLED_S[e]; }
			inline bool areFlowSendingPinsSourceReachable(const Hyperedge e) const { return IN_REACHED_S[e]; }
			inline void settleFlowSendingPins(const Hyperedge e) { Assert(!areFlowSendingPinsSources(e)); IN_SETTLED_S.set(e); }
			inline void reachFlowSendingPins(const Hyperedge e) { Assert(!areFlowSendingPinsSourceReachable(e)); IN_REACHED_S.set(e); }

			void resetSourceReachableToSource() {
				IN_REACHED_S = IN_SETTLED_S;
				OUT_REACHED_S = OUT_SETTLED_S;
			}

			void flipViewDirection() {
				std::swap(IN_SETTLED_S, OUT_SETTLED_T);
				std::swap(OUT_SETTLED_S, IN_SETTLED_T);
				std::swap(IN_REACHED_S, OUT_REACHED_T);
				std::swap(OUT_REACHED_S, IN_REACHED_T);
			}

			void verifyDisjoint() const {
				Assert((OUT_REACHED_S & OUT_REACHED_T).none());
				Assert((OUT_SETTLED_S & OUT_SETTLED_T).none());
				Assert((IN_REACHED_S & IN_REACHED_T).none());
				Assert((IN_SETTLED_S & IN_SETTLED_T).none());
			}

			void verifySettledIsSubsetOfReachable() const {
				Assert(OUT_SETTLED_S.is_subset_of(OUT_REACHED_S));
				Assert(IN_SETTLED_S.is_subset_of(IN_REACHED_S));
				Assert(OUT_SETTLED_T.is_subset_of(OUT_REACHED_T));
				Assert(IN_SETTLED_T.is_subset_of(IN_REACHED_T));
			}

			protected:
			BitVector IN_SETTLED_S, OUT_SETTLED_S, IN_REACHED_S, OUT_REACHED_S;
			BitVector IN_SETTLED_T, OUT_SETTLED_T, IN_REACHED_T, OUT_REACHED_T;

		};

}