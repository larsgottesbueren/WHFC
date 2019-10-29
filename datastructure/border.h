#pragma once

#include <vector>
#include "bitvector.h"
#include "../definitions.h"
#include "../util/random.h"

#include "../util/filter.h"
#include "flow_hypergraph.h"

#include "../util/range.h"
#include "../util/concatenated_range.h"
#include "../util/sub_range.h"

namespace whfc {
	
	template<typename T, bool trackElements>
	class PersistentSet {
	private:
		bool persistentMode = true;
		size_t persistent_front = 0, persistent_end = 0, non_persistent_front = 0;
		BitVector was_added;
		std::vector<T> c;
		
		sub_range<std::vector<T>> persistent_entries() const {
			return sub_range(c, persistent_front, persistent_end);
		}
		
		sub_range<std::vector<T>> non_persistent_entries() const {
			return sub_range(c, non_persistent_front, c.size());
		}
		
	public:
		explicit PersistentSet(const size_t nT) : was_added(nT) { }
		
		bool wasAdded(const T x) const {
			return was_added[x];
		}
		
		void add(const T& x) {
			Assert(!wasAdded(x));
			was_added.set(x);
			if (trackElements || !persistentMode)
				c.push_back(x);
		}
		
		// delete non persistent entries, even those removed from the list
		// and recover the persistent entries, even those removed from the list
		void recover() {
			while (c.size() > persistent_end) {
				was_added.reset(c.back());
				c.pop_back();
			}
			non_persistent_front = persistent_end;
			persistent_front = 0;
		}
		
		void lockInPersistentEntries() {
			persistentMode = false;
			persistent_end = c.size();
			non_persistent_front = persistent_end;
		}
		
		template<typename Predicate>
		void cleanUp(Predicate p) {
			if (persistentMode) {
				util::remove_if_inplace(c, p);
			}
			else {
				util::move_to_front_if(c, persistent_front, persistent_end, p);
				util::move_to_front_if(c, non_persistent_front, c.size(), p);
			}
		}
		
		void reset(size_t newN) {
			was_added.reset(0, newN);
			c.clear();
			persistent_front = 0;
			persistent_end = 0;
			non_persistent_front = 0;
			persistentMode = true;
		}
		
		bool empty() const {
			return persistent_front == persistent_end && non_persistent_front == c.size();
		}
		
		T popRandomEntryPreferringPersistent() {
			Assert(!persistentMode);
			Assert(!empty());
			if (persistent_front < persistent_end) {
				size_t ind = Random::randomIndex(persistent_front, persistent_end - 1);
				std::swap(c[ind], c[persistent_front]);
				return c[persistent_front++];
			}
			else {
				size_t ind = Random::randomIndex(non_persistent_front, c.size() - 1);
				std::swap(c[ind], c[non_persistent_front]);
				return c[non_persistent_front++];
			}
		}
		
		auto entries() const {
			return concatenated_range< sub_range<std::vector<T>>, T  >(persistent_entries(), non_persistent_entries());
		}
		
		const std::vector<T>& entries_in_persistent_mode() const {
			Assert(persistentMode);
			return c;
		}
		
		std::vector<T> copy() const {
			std::vector<T> c;
			for (const T& x : entries())
				c.push_back(x);
			return c;
		}
	};
	
	template<typename T, bool trackElements>
	class Borders {
	public:
		explicit Borders(size_t nT) : sourceSide(nT), targetSide(nT) { }

		void flipViewDirection() {
			std::swap(sourceSide, targetSide);
		}

		PersistentSet<T, trackElements> sourceSide, targetSide;

		void reset(const size_t newN) {
			sourceSide.reset(newN);
			targetSide.reset(newN);
		}
		
		void enterMostBalancedCutMode() {
			sourceSide.lockInPersistentEntries();
			targetSide.lockInPersistentEntries();
		}
		
		void resetForMostBalancedCut() {
			sourceSide.recover();
			targetSide.recover();
		}
	};

	using NodeBorders = Borders<Node, true>;

	//track hyperedges only for assertions in debug mode
#ifndef NDEBUG
	using HyperedgeCuts = Borders<Hyperedge, true>;
#else
	using HyperedgeCuts = Borders<Hyperedge, false>;
#endif

}
