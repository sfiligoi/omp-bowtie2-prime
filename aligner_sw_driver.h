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

/*
 *  aligner_sw_driver.h
 *
 * REDUNDANT SEED HITS
 *
 * We say that two seed hits are redundant if they trigger identical
 * seed-extend dynamic programming problems.  Put another way, they both lie on
 * the same diagonal of the overall read/reference dynamic programming matrix.
 * Detecting redundant seed hits is simple when the seed hits are ungapped.  We
 * do this after offset resolution but before the offset is converted to genome
 * coordinates (see uses of the seenDiags1_/seenDiags2_ fields for examples).
 *
 * REDUNDANT ALIGNMENTS
 *
 * In an unpaired context, we say that two alignments are redundant if they
 * share any cells in the global DP table.  Roughly speaking, this is like
 * saying that two alignments are redundant if any read character aligns to the
 * same reference character (same reference sequence, same strand, same offset)
 * in both alignments.
 *
 * In a paired-end context, we say that two paired-end alignments are redundant
 * if the mate #1s are redundant and the mate #2s are redundant.
 *
 * How do we enforce this?  In the unpaired context, this is relatively simple:
 * the cells from each alignment are checked against a set containing all cells
 * from all previous alignments.  Given a new alignment, for each cell in the
 * new alignment we check whether it is in the set.  If there is any overlap,
 * the new alignment is rejected as redundant.  Otherwise, the new alignment is
 * accepted and its cells are added to the set.
 *
 * Enforcement in a paired context is a little trickier.  Consider the
 * following approaches:
 *
 * 1. Skip anchors that are redundant with any previous anchor or opposite
 *    alignment.  This is sufficient to ensure no two concordant alignments
 *    found are redundant.
 *
 * 2. Same as scheme 1, but with a "transitive closure" scheme for finding all
 *    concordant pairs in the vicinity of an anchor.  Consider the AB/AC
 *    scenario from the previous paragraph.  If B is the anchor alignment, we
 *    will find AB but not AC.  But under this scheme, once we find AB we then
 *    let B be a new anchor and immediately look for its opposites.  Likewise,
 *    if we find any opposite, we make them anchors and continue searching.  We
 *    don't stop searching until every opposite is used as an anchor.
 *
 * 3. Skip anchors that are redundant with any previous anchor alignment (but
 *    allow anchors that are redundant with previous opposite alignments).
 *    This isn't sufficient to avoid redundant concordant alignments.  To avoid
 *    redundant concordants, we need an additional procedure that checks each
 *    new concordant alignment one-by-one against a list of previous concordant
 *    alignments to see if it is redundant.
 *
 * We take approach 1.
 */

#ifndef ALIGNER_SW_DRIVER_H_
#define ALIGNER_SW_DRIVER_H_

#include <utility>
#include "ds.h"
#include "aligner_seed.h"
#include "aligner_sw.h"
#include "aligner_cache.h"
#include "reference.h"
#include "group_walk.h"
#include "bt2_idx.h"
#include "mem_ids.h"
#include "aln_sink.h"
#include "pe.h"
#include "ival_list.h"
#include "simple_func.h"
#include "random_util.h"

// Hardcode for now. May want to pass in Makefile
// The chosen values are appropriate for the default SHOTGUN paramers
#define ALN_MAX_ITER 1600


struct SeedPos {

	SeedPos() : fw(false), offidx(0), rdoff(0), seedlen(0) { }

	SeedPos(
		bool fw_,
		uint32_t offidx_,
		uint32_t rdoff_,
		uint32_t seedlen_)
	{
		init(fw_, offidx_, rdoff_, seedlen_);
	}
	
	void init(
		bool fw_,
		uint32_t offidx_,
		uint32_t rdoff_,
		uint32_t seedlen_)
	{
		fw      = fw_;
		offidx  = offidx_;
		rdoff   = rdoff_;
		seedlen = seedlen_;
	}
	
	bool operator<(const SeedPos& o) const {
		if(offidx < o.offidx)   return true;
		if(offidx > o.offidx)   return false;
		if(rdoff < o.rdoff)     return true;
		if(rdoff > o.rdoff)     return false;
		if(seedlen < o.seedlen) return true;
		if(seedlen > o.seedlen) return false;
		if(fw && !o.fw)         return true;
		if(!fw && o.fw)         return false;
		return false;
	}
	
	bool operator>(const SeedPos& o) const {
		if(offidx < o.offidx)   return false;
		if(offidx > o.offidx)   return true;
		if(rdoff < o.rdoff)     return false;
		if(rdoff > o.rdoff)     return true;
		if(seedlen < o.seedlen) return false;
		if(seedlen > o.seedlen) return true;
		if(fw && !o.fw)         return false;
		if(!fw && o.fw)         return true;
		return false;
	}
	
	bool operator==(const SeedPos& o) const {
		return fw == o.fw && offidx == o.offidx &&
		       rdoff == o.rdoff && seedlen == o.seedlen;
	}

	bool fw;
	uint32_t offidx;
	uint32_t rdoff;
	uint32_t seedlen;
};

/**
 * An SATuple along with the associated seed position.
 */
struct SATupleAndPos {
	
	SATuple sat;    // result for this seed hit
	SeedPos pos;    // seed position that yielded the range this was taken from
	size_t  origSz; // size of range this was taken from
	
	bool operator<(const SATupleAndPos& o) const {
		if(sat < o.sat) return true;
		if(sat > o.sat) return false;
		return pos < o.pos;
	}

	bool operator==(const SATupleAndPos& o) const {
		return sat == o.sat && pos == o.pos;
	}
};

/**
 * Return values from extendSeeds and extendSeedsPaired.
 */
enum {
	// All end-to-end and seed hits were examined
	// The policy does not need us to look any further
	EXTEND_EXHAUSTED_CANDIDATES = 1,
	EXTEND_POLICY_FULFILLED,
	// We stopped because we reached a point where the only remaining
	// alignments of interest have perfect scores, but we already investigated
	// perfect alignments
	EXTEND_PERFECT_SCORE,
	// We stopped because we ran up against a limit on how much work we should
	// do for one set of seed ranges, e.g. the limit on number of consecutive
	// unproductive DP extensions
	EXTEND_EXCEEDED_SOFT_LIMIT,
	// We stopped because we ran up against a limit on how much work we should
	// do for overall before giving up on a mate
	EXTEND_EXCEEDED_HARD_LIMIT
};

/**
 * Data structure encapsulating a range that's been extended out in two
 * directions.
 */
struct ExtendRange {

	void init(size_t off_, size_t len_, size_t sz_) {
		off = off_; len = len_; sz = sz_;
	}

	size_t off; // offset of extended region
	size_t len; // length between extremes of extended region
	size_t sz;  // # of elements in SA range
};

class SwDriver {

public:

	SwDriver() :
		seenDiags1_(DP_CAT),
		seenDiags2_(DP_CAT),
		redAnchor_(DP_CAT),
		gwstate_(GW_CAT),
		reportOverhangs(gReportOverhangs) { }

	// Allow move operator
	SwDriver(SwDriver&& other) = default;
	SwDriver& operator=(SwDriver&& other) noexcept = delete;

	// but not copy
	SwDriver(const SwDriver& other) = delete;
	SwDriver& operator=(const SwDriver& other) = delete;

	void set_alloc(BTAllocator *alloc, bool propagate_alloc=true) {
		seenDiags1_.set_alloc(alloc,propagate_alloc);
		seenDiags2_.set_alloc(alloc,propagate_alloc);
		redAnchor_.set_alloc(alloc,propagate_alloc);
		resGap_.set_alloc(alloc,propagate_alloc);
		gwstate_.set_alloc(alloc,propagate_alloc);
	}

	void set_alloc(std::pair<BTAllocator *, bool> arg) {
		set_alloc(arg.first, arg.second);
	}

	/**
	 * Given seed results, set up all of our state for resolving and keeping
	 * track of reference offsets for hits.
	 */
	void prioritizeSATups(
		const SeedResults& sh,       // seed hits to extend into full alignments
		const Ebwt& ebwtFw,          // BWT
		const BitPairReference& ref, // Reference strings
		int seedmms,                 // # seed mismatches allowed
		bool lensq,                  // square extended length
		bool szsq,                   // square SA range size
		AlignmentCacheInterface ca,  // alignment cache for seed hits
		bool all);                   // report all hits?

	/**
	 * Given a collection of SeedHits for a single read, extend seed alignments
	 * into full alignments.  Where possible, try to avoid redundant offset
	 * lookups and dynamic programming problems.  Optionally report alignments
	 * to a AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeeds(
		const SeedResults& sh,       // seed hits to extend into full alignments
		const Ebwt& ebwtFw,          // BWT
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dynamic programming aligner
		const Scoring& sc,           // scoring scheme
		const int seedmms,           // # mismatches allowed in seed
		const int seedlen,           // length of seed
		const int seedival,          // interval between seeds
		const TAlScore minsc,        // minimum score for anchor
		const int nceil,             // maximum # Ns permitted in ref portion
		const size_t maxhalf,        // maximum width on one side of DP table
		const bool doExtend,         // do seed extension
		const bool enable8,          // use 8-bit SSE where possible
		PerReadMetrics& prm,         // per-read metrics
		AlnSinkWrap* mhs,            // HitSink for multiseed-style aligner
		bool reportImmediately,      // whether to report hits immediately to mhs
		bool& exhaustive);

#ifdef SUPPORT_PAIRED
// NOTE: Unsupported, likely does not work
	/**
	 * Given a collection of SeedHits for a read pair, extend seed
	 * alignments into full alignments and then look for the opposite
	 * mate using dynamic programming.  Where possible, try to avoid
	 * redundant offset lookups.  Optionally report alignments to a
	 * AlnSinkWrap object as they are discovered.
	 *
	 * If 'reportImmediately' is true, returns true iff a call to
	 * mhs->report() returned true (indicating that the reporting
	 * policy is satisfied and we can stop).  Otherwise, returns false.
	 */
	int extendSeedsPaired(
		Read& rd,                    // mate to align as anchor
		Read& ord,                   // mate to align as opposite
		bool anchor1,                // true iff anchor mate is mate1
		bool oppFilt,                // true iff opposite mate was filtered out
		SeedResults& sh,             // seed hits for anchor
		const Ebwt& ebwtFw,          // BWT
		const Ebwt* ebwtBw,          // BWT'
		const BitPairReference& ref, // Reference strings
		SwAligner& swa,              // dyn programming aligner for anchor
		SwAligner& swao,             // dyn programming aligner for opposite
		const Scoring& sc,           // scoring scheme
		const PairedEndPolicy& pepol,// paired-end policy
		int seedmms,                 // # mismatches allowed in seed
		int seedlen,                 // length of seed
		int seedival,                // interval between seeds
		TAlScore& minsc,             // minimum score for anchor
		TAlScore& ominsc,            // minimum score for opposite
		int nceil,                   // max # Ns permitted in ref for anchor
		int onceil,                  // max # Ns permitted in ref for opposite
		bool nofw,                   // don't align forward read
		bool norc,                   // don't align revcomp read
		size_t maxhalf,              // maximum width on one side of DP table
		bool doUngapped,             // do ungapped alignment
		size_t maxIters,             // stop after this many seed-extend loop iters
		size_t maxUg,                // max # ungapped extends
		size_t maxDp,                // max # DPs
		size_t maxEeStreak,          // stop after streak of this many end-to-end fails
		size_t maxUgStreak,          // stop after streak of this many ungap fails
		size_t maxDpStreak,          // stop after streak of this many dp fails
		size_t maxMateStreak,        // stop seed range after N mate-find fails
		bool doExtend,               // do seed extension
		bool enable8,                // use 8-bit SSE where possible
		int tighten,                 // -M score tightening mode
		AlignmentCacheInterface cs,  // alignment cache for seed hits
		RandomSource& rnd,           // pseudo-random source
		PerReadMetrics& prm,         // per-read metrics for anchor
		AlnSinkWrap* msink,          // AlnSink wrapper for multiseed-style aligner
		bool swMateImmediately,      // whether to look for mate immediately
		bool reportImmediately,      // whether to report hits immediately to msink
		bool discord,                // look for discordant alignments?
		bool mixed,                  // look for unpaired as well as paired alns?
		bool& exhaustive);
#endif

	/**
	 * Prepare for a new read.
	 */
	void nextRead(size_t mate1len) {
		redAnchor_.reset();
		seenDiags1_.reset();
		seenDiags2_.reset();
		seedExRangeFw_.clear();
		seedExRangeRc_.clear();
		size_t maxlen = mate1len;
		redAnchor_.init(maxlen);
	}

protected:

	// if range as <= nsm elts, it's "small"
	constexpr static size_t nsm = 5;

	size_t nelt_;  // set by prioritizeSATups

	DList<SATupleAndPos, ALN_MAX_ITER> satpos_;  // holds SATuple, SeedPos pairs
	DList<GroupWalk2S<TSlice>, 1 > gws_;   // list of GroupWalks; only one used at a time
	
	// Ranges that we've extended through when extending seed hits
	EList<ExtendRange> seedExRangeFw_;
	EList<ExtendRange> seedExRangeRc_;

	// Data structures encapsulating the diagonals that have already been used
	// to seed alignment for mate 1 and mate 2.
	EIvalMergeListBinned seenDiags1_;
	EIvalMergeListBinned seenDiags2_;

	// For weeding out redundant alignments
	RedundantAlns  redAnchor_;  // database of cells used for anchor alignments

	// For holding results for anchor (res_) and opposite (ores_) mates
	SwResult       resGap_;    // temp holder for alignment result
	
	GroupWalkState gwstate_;   // some per-thread state shared by all GroupWalks
	const bool reportOverhangs;

};

#endif /*ALIGNER_SW_DRIVER_H_*/
