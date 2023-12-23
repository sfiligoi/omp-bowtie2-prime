/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALIGNER_SW_COMMON_H_
#define ALIGNER_SW_COMMON_H_

#include "aligner_result.h"

/**
 * Encapsulates the result of a dynamic programming alignment.  In our
 * case, the result is a combination of:
 *
 * 1. All the nucleotide edits
 * 2. All the "edits" where an ambiguous reference char is resolved to
 *    an unambiguous char.
 * 3. The score of the best alginment
 * 4. The score of the second-best alignment
 *
 * Having scores for the best and second-best alignments gives us an
 * idea of where gaps may make reassembly beneficial.
 */
struct SwResult {

	SwResult() :
		alres(),
		sws(0),
		swcups(0),
		swrows(0),
		swskiprows(0),
		swskip(0),
		swsucc(0),
		swfail(0),
		swbts(0)
	{ }

	void set_alloc(BTAllocator *alloc, bool propagate_alloc=true) {
		alres.set_alloc(alloc,propagate_alloc);
	}

	void set_alloc(std::pair<BTAllocator *, bool> arg) {
		alres.set_alloc(arg);
	}

	/**
	 * Clear all contents.
	 */
	void reset() {
		sws = swcups = swrows = swskiprows = swskip = swsucc =
		swfail = swbts = 0;
		alres.reset();
	}
	
	/**
	 * Reverse all edit lists.
	 */
	void reverse() {
		alres.reverseEdits();
	}
	
	/**
	 * Return true iff no result has been installed.
	 */
	bool empty() const {
		return alres.empty();
	}
	
#ifndef NDEBUG
	/**
	 * Check that result is internally consistent.
	 */
	bool repOk() const {
		assert(alres.repOk());
		return true;
	}

	/**
	 * Check that result is internally consistent w/r/t read.
	 */
	bool repOk(const Read& rd) const {
		assert(alres.repOk(rd));
		return true;
	}
#endif

	AlnRes alres;
	uint64_t sws;    // # DP problems solved
	uint64_t swcups; // # DP cell updates
	uint64_t swrows; // # DP row updates
	uint64_t swskiprows; // # skipped DP row updates (b/c no valid alignments can go thru row)
	uint64_t swskip; // # DP problems skipped by sse filter
	uint64_t swsucc; // # DP problems resulting in alignment
	uint64_t swfail; // # DP problems not resulting in alignment
	uint64_t swbts;  // # DP backtrace steps
};

// The various ways that one might backtrack from a later cell (either oall,
// rdgap or rfgap) to an earlier cell
enum {
	SW_BT_OALL_DIAG,         // from oall cell to oall cell
	SW_BT_OALL_REF_OPEN,     // from oall cell to oall cell
	SW_BT_OALL_READ_OPEN,    // from oall cell to oall cell
	SW_BT_RDGAP_EXTEND,      // from rdgap cell to rdgap cell
	SW_BT_RFGAP_EXTEND       // from rfgap cell to rfgap cell
};

#endif /*def ALIGNER_SW_COMMON_H_*/
