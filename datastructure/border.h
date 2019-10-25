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
	class Border {
	public:
		explicit Border(size_t nT) : was_added(nT) { }
		static constexpr bool log = true;
		using BorderRange = concatenated_range< sub_range<std::vector<T>>, T >;
		
		BitVector was_added;
		std::vector<T> persistent_entries, non_persistent_entries;
		size_t currentNumberOfPersistentEntries = 0;
		size_t currentNumberOfNonPersistentEntries = 0;
		bool nonPersistentMode = false;
		
		bool wasAdded(const T x) const {
			return was_added[x];
		}
		
		void add(T x) {
			Assert(!wasAdded(x));
			was_added.set(x);
			if (nonPersistentMode) {
				non_persistent_entries.push_back(x);
				//set number of non persistent entries.
				
				//TODO push back and moving stale elements to the back does not work! think of something different.
			}
			else {
				if constexpr (trackElements) {
					persistent_entries.push_back(x);
				}
			}
		}
		
		template<typename Predicate>
		void cleanUp(Predicate pred) {
			if (nonPersistentMode) {
				LOGGER << "non persistent mode clean up";
				LOGGER << "clean persistent entries" << V(persistent_entries.size()) << V(currentNumberOfPersistentEntries);
				util::move_to_end_if(persistent_entries, currentNumberOfPersistentEntries, pred);
				LOGGER << "clean non persistent entries" << V(non_persistent_entries.size()) << V(currentNumberOfNonPersistentEntries);
				util::move_to_end_if(non_persistent_entries, currentNumberOfNonPersistentEntries, pred);
			}
			else {
				if constexpr (trackElements) {
					LOGGER << "clean up in non persistent mode. old size" << persistent_entries.size();
					util::remove_if_inplace(persistent_entries, pred);
					LOGGER << "new size" << persistent_entries.size();
				}
			}
		}
		
		void reset(const size_t newN) {	// not a proper clearlist, we're deleting some entries without resetting the was_added bit
			was_added.reset(0, newN);
			persistent_entries.clear();
			non_persistent_entries.clear();
			currentNumberOfPersistentEntries = 0;
			currentNumberOfNonPersistentEntries = 0;
			nonPersistentMode = false;
		}
		
		bool empty() const {
			if (nonPersistentMode)
				return currentNumberOfNonPersistentEntries == 0 && currentNumberOfPersistentEntries == 0;
			else
				return persistent_entries.empty();
		}
		
		T popRandomEntryPreferringPersistent() {
			Assert(nonPersistentMode);
			Assert(!empty());
			if (currentNumberOfPersistentEntries > 0) {
				size_t ind = Random::randomIndex(0, --currentNumberOfPersistentEntries );
				std::swap(persistent_entries[ind], persistent_entries[currentNumberOfPersistentEntries]);
				return persistent_entries[currentNumberOfPersistentEntries];
			}
			else {
				size_t ind = Random::randomIndex(0, --currentNumberOfNonPersistentEntries);
				std::swap(non_persistent_entries[ind], non_persistent_entries[currentNumberOfNonPersistentEntries]);
				return non_persistent_entries[currentNumberOfNonPersistentEntries];
			}
		}
		
		void enterNonPersistentMode() {
			reinsertPersistentEntries();
			nonPersistentMode = true;
		}
		
		void removeNonPersistentEntries() {
			for (const T x : non_persistent_entries)
				was_added.reset(x);
			non_persistent_entries.clear();
			currentNumberOfNonPersistentEntries = 0;
		}
		
		void reinsertPersistentEntries() {
			currentNumberOfPersistentEntries = persistent_entries.size();
		}
		
		void resetToPersistentState() {
			removeNonPersistentEntries();
			reinsertPersistentEntries();
		}
		
		BorderRange entries() const {
			if (nonPersistentMode)
				return BorderRange( sub_range(persistent_entries, 0, currentNumberOfPersistentEntries), sub_range(non_persistent_entries, 0, currentNumberOfNonPersistentEntries) );
			else
				return BorderRange( sub_range(persistent_entries, 0, persistent_entries.size()), sub_range(non_persistent_entries, 0 , 0) );
		}
		
		std::vector<T> copy() {
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

		Border<T, trackElements> sourceSide, targetSide;

		void reset(const size_t newN) {
			sourceSide.reset(newN);
			targetSide.reset(newN);
		}
		
		void enterMostBalancedCutMode() {
			sourceSide.enterNonPersistentMode();
			targetSide.enterNonPersistentMode();
		}
		
		void resetForMostBalancedCut() {
			sourceSide.resetToPersistentState();
			targetSide.resetToPersistentState();
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
