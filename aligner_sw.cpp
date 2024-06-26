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

#include <limits>
// -- BTL remove --
//#include <stdlib.h>
//#include <sys/time.h>
// -- --
#include "aligner_sw.h"
#include "aligner_result.h"
#include "search_globals.h"
#include "scoring.h"
#include "mask.h"

/**
 * Initialize with a new read.
 */
void SwAligner::initRead(
	const BTDnaString& rdfw, // forward read sequence
	const BTDnaString& rdrc, // revcomp read sequence
	const BTString& qufw,    // forward read qualities
	const BTString& qurc,    // reverse read qualities
	size_t rdi,              // offset of first read char to align
	size_t rdf,              // offset of last read char to align
	const Scoring& sc)       // scoring scheme
{
	assert_gt(rdf, rdi);
	int nceil = sc.nCeil.f<int>((double)rdfw.length());
	rdfw_    = &rdfw;      // read sequence
	rdrc_    = &rdrc;      // read sequence
	qufw_    = &qufw;      // read qualities
	qurc_    = &qurc;      // read qualities
	rdi_     = rdi;        // offset of first read char to align
	rdf_     = rdf;        // offset of last read char to align
	sc_      = &sc;        // scoring scheme
	nceil_   = nceil;      // max # Ns allowed in ref portion of aln
	readSse16_ = false;    // true -> sse16 from now on for this read
	initedRead_ = true;
#ifndef NO_SSE
	sseU8fwBuilt_  = false;  // built fw query profile, 8-bit score
	sseU8rcBuilt_  = false;  // built rc query profile, 8-bit score
	sseI16fwBuilt_ = false;  // built fw query profile, 16-bit score
	sseI16rcBuilt_ = false;  // built rc query profile, 16-bit score
#endif
}

/**
 * Initialize with a new alignment problem.
 */
void SwAligner::initRef(
	bool fw,               // whether to forward or revcomp read is aligning
	TRefId refidx,         // id of reference aligned against
	const DPRect& rect,    // DP rectangle
	char *rf,              // reference sequence
	size_t rfi,            // offset of first reference char to align to
	size_t rff,            // offset of last reference char to align to
	TRefOff reflen,        // length of reference sequence
	const Scoring& sc,     // scoring scheme
	TAlScore minsc,        // minimum score
	bool enable8,          // use 8-bit SSE if possible?
	size_t cminlen,        // minimum length for using checkpointing scheme
	size_t cpow2,          // interval b/t checkpointed diags; 1 << this
	bool doTri,            // triangular mini-fills?
	bool extend)           // is this a seed extension?
{
	size_t readGaps = sc.maxReadGaps(minsc, rdfw_->length());
	size_t refGaps  = sc.maxRefGaps(minsc, rdfw_->length());
	assert_geq(readGaps, 0);
	assert_geq(refGaps, 0);
	assert_gt(rff, rfi);
	rdgap_       = readGaps;  // max # gaps in read
	rfgap_       = refGaps;   // max # gaps in reference
	state_       = STATE_INITED;
	fw_          = fw;       // orientation
	rd_          = fw ? rdfw_ : rdrc_; // read sequence
	qu_          = fw ? qufw_ : qurc_; // quality sequence
	refidx_      = refidx;   // id of reference aligned against
	rf_          = rf;       // reference sequence
	rfi_         = rfi;      // offset of first reference char to align to
	rff_         = rff;      // offset of last reference char to align to
	reflen_      = reflen;   // length of entire reference sequence
	rect_        = &rect;    // DP rectangle
	minsc_       = minsc;    // minimum score
	cural_       = 0;        // idx of next alignment to give out
	initedRef_   = true;     // indicate we've initialized the ref portion
	enable8_     = enable8;  // use 8-bit SSE if possible?
	extend_      = extend;   // true iff this is a seed extension
}
	
/**
 * Given a read, an alignment orientation, a range of characters in a referece
 * sequence, and a bit-encoded version of the reference, set up and execute the
 * corresponding dynamic programming problem.
 *
 * The caller has already narrowed down the relevant portion of the reference
 * using, e.g., the location of a seed hit, or the range of possible fragment
 * lengths if we're searching for the opposite mate in a pair.
 */
void SwAligner::initRef(
	bool fw,               // whether to forward or revcomp read is aligning
	TRefId refidx,         // reference aligned against
	const DPRect& rect,    // DP rectangle
	const BitPairReference& refs, // Reference strings
	TRefOff reflen,        // length of reference sequence
	const Scoring& sc,     // scoring scheme
	TAlScore minsc,        // minimum score
	bool enable8,          // use 8-bit SSE if possible?
	size_t cminlen,        // minimum length for using checkpointing scheme
	size_t cpow2,          // interval b/t checkpointed diags; 1 << this
	bool doTri,            // triangular mini-fills?
	bool extend,           // true iff this is a seed extension
	size_t  upto,          // count the number of Ns up to this offset
	size_t& nsUpto)        // output: the number of Ns up to 'upto'
{
	TRefOff rfi = rect.refl;
	TRefOff rff = rect.refr + 1;
	assert_gt(rff, rfi);
	// Capture an extra reference character outside the rectangle so that we
	// can check matches in the next column over to the right
	rff++;
	// rflen = full length of the reference substring to consider, including
	// overhang off the boundaries of the reference sequence
	const size_t rflen = (size_t)(rff - rfi);
	// Figure the number of Ns we're going to add to either side
	size_t leftNs  =
		(rfi >= 0               ? 0 : (size_t)std::abs(static_cast<long>(rfi)));
	leftNs = min(leftNs, rflen);
	size_t rightNs =
		(rff <= (TRefOff)reflen ? 0 : (size_t)std::abs(static_cast<long>(rff - reflen)));
	rightNs = min(rightNs, rflen);
	// rflenInner = length of just the portion that doesn't overhang ref ends
	assert_geq(rflen, leftNs + rightNs);
	const size_t rflenInner = rflen - (leftNs + rightNs);
#ifndef NDEBUG
	bool haveRfbuf2 = false;
	EList<char> rfbuf2(rflen);
	// This is really slow, so only do it some of the time
	if((rand() % 10) == 0) {
		TRefOff rfii = rfi;
		for(size_t i = 0; i < rflen; i++) {
			if(rfii < 0 || (TRefOff)rfii >= reflen) {
				rfbuf2.push_back(4);
			} else {
				rfbuf2.push_back(refs.getBase(refidx, (size_t)rfii));
			}
			rfii++;
		}
		haveRfbuf2 = true;
	}
#endif
	// rfbuf_ = uint32_t list large enough to accommodate both the reference
	// sequence and any Ns we might add to either side.
	rfwbuf_.resize((rflen + 16) / 4);
	int offset = refs.getStretch(
		rfwbuf_.ptr(),               // buffer to store words in
		refidx,                      // which reference
		(rfi < 0) ? 0 : (size_t)rfi, // starting offset (can't be < 0)
		rflenInner                   // length to grab (exclude overhang)
		ASSERT_ONLY(, tmp_destU32_));// for BitPairReference::getStretch()
	assert_leq(offset, 16);
	rf_ = (char*)rfwbuf_.ptr() + offset;
	// Shift ref chars away from 0 so we can stick Ns at the beginning
	if(leftNs > 0) {
		// Slide everyone down
		for(size_t i = rflenInner; i > 0; i--) {
			rf_[i+leftNs-1] = rf_[i-1];
		}
		// Add Ns
		for(size_t i = 0; i < leftNs; i++) {
			rf_[i] = 4;
		}
	}
	if(rightNs > 0) {
		// Add Ns to the end
		for(size_t i = 0; i < rightNs; i++) {
			rf_[i + leftNs + rflenInner] = 4;
		}
	}
#ifndef NDEBUG
	// Sanity check reference characters
	for(size_t i = 0; i < rflen; i++) {
		assert(!haveRfbuf2 || rf_[i] == rfbuf2[i]);
		assert_range(0, 4, (int)rf_[i]);
	}
#endif
	// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
	// nucleotides (IUPAC codes) have more than one mask bit set.  If a
	// reference scanner was provided, use it to opportunistically resolve seed
	// hits.
	nsUpto = 0;
	for(size_t i = 0; i < rflen; i++) {
		// rf_[i] gets mask version of refence char, with N=16
		if(i < upto && rf_[i] > 3) {
			nsUpto++;
		}
		rf_[i] = (1 << rf_[i]);
	}
	// Correct for having captured an extra reference character
	rff--;
	initRef(
		fw,          // whether to forward or revcomp read is aligning
		refidx,      // id of reference aligned against
		rect,        // DP rectangle
		rf_,         // reference sequence, wrapped up in BTString object
		0,           // use the whole thing
		(size_t)(rff - rfi), // ditto
		reflen,      // reference length
		sc,          // scoring scheme
		minsc,       // minimum score
		enable8,     // use 8-bit SSE if possible?
		cminlen,     // minimum length for using checkpointing scheme
		cpow2,       // interval b/t checkpointed diags; 1 << this
		doTri,       // triangular mini-fills?
		extend);     // true iff this is a seed extension
}

/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.
 */
bool SwAligner::align(
	TAlScore& best)    // best alignment score observed in DP matrix
{
	assert(initedRef() && initedRead());
	assert_eq(STATE_INITED, state_);
	state_ = STATE_ALIGNED;
	// Reset solutions lists
	btncand_.clear();
	best = std::numeric_limits<TAlScore>::min();
	sse8succ_ = sse16succ_ = false;
	int flag = 0;
	size_t rdlen = rdf_ - rdi_;
	// assert(!checkpointed)
	assert(sc_->monotone);
	{
		// End-to-end
		if(enable8_ && !readSse16_ && minsc_ >= -254) {
			// 8-bit end-to-end
			best = alignNucleotidesEnd2EndSseU8(flag, false);
			sse8succ_ = (flag == 0);
		} else {
			// 16-bit end-to-end
			best = alignNucleotidesEnd2EndSseI16(flag, false);
			sse16succ_ = (flag == 0);
		}
	}
	assert(repOk());
	cural_ = 0;
	if(best == MIN_I64 || best < minsc_) {
		return false;
	}
	// assert(!gathered)
	{
		// Look for solutions using SSE matrix
		assert(sse8succ_ || sse16succ_);
		{
			if(sse8succ_) {
				gatherCellsNucleotidesEnd2EndSseU8(best);
			} else {
				gatherCellsNucleotidesEnd2EndSseI16(best);
			}
		}
	}
	if(!btncand_.empty()) {
		btncand_.sort();
	}
	return !btncand_.empty();
}

/**
 * Populate the given SwResult with information about the "next best"
 * alignment if there is one.  If there isn't one, false is returned.  Note
 * that false might be returned even though a call to done() would have
 * returned false.
 */
bool SwAligner::nextAlignment(
	SwResult& res,
	TAlScore minsc,
	RandomSource& rnd)
{
	assert(initedRead() && initedRef());
	assert_eq(STATE_ALIGNED, state_);
	assert(repOk());
	if(done()) {
		res.reset();
		return false;
	}
	assert(!done());
	size_t off = 0, nbts = 0;
	assert_lt(cural_, btncand_.size());
	assert(res.repOk());
	// For each candidate cell that we should try to backtrack from...
	const size_t candsz = btncand_.size();
	size_t SQ = dpRows() >> 4;
	if(SQ == 0) SQ = 1;
	size_t rdlen = rdf_ - rdi_;
	// assert(!checkpointed)
	while(cural_ < candsz) {
		// Doing 'continue' anywhere in here simply causes us to move on to the
		// next candidate
		if(btncand_[cural_].score < minsc) {
			btncand_[cural_].fate = BT_CAND_FATE_FILT_SCORE;
			nbtfiltsc_++; cural_++; continue;
		}
		nbts = 0;
		assert(sse8succ_ || sse16succ_);
		size_t row = btncand_[cural_].row;
		size_t col = btncand_[cural_].col;
		assert_lt(row, dpRows());
		assert_lt((TRefOff)col, rff_-rfi_);
		if(sse16succ_) {
			SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
			if (d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
				// Skipping this candidate because a previous candidate already
				// moved through this cell
				btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
				//cerr << "  skipped becuase starting cell was covered" << endl;
				nbtfiltst_++; cural_++; continue;
			}
		} else if(sse8succ_) {
			SSEData& d = fw_ ? sseU8fw_ : sseU8rc_;
			if (d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
				// Skipping this candidate because a previous candidate already
				// moved through this cell
				btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
				//cerr << "  skipped becuase starting cell was covered" << endl;
				nbtfiltst_++; cural_++; continue;
			}
		}
		assert(sc_->monotone);
		{
			bool ret = false;
			if(sse8succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				rnd.init(reseed);
				res.reset();
				ret = backtraceNucleotidesEnd2EndSseU8(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				rnd.init(reseed+1); // debug/release pseudo-randoms in lock step
			} else if(sse16succ_) {
				uint32_t reseed = rnd.nextU32() + 1;
				res.reset();
				ret = backtraceNucleotidesEnd2EndSseI16(
						btncand_[cural_].score, // in: expected score
						res,    // out: store results (edits and scores) here
						off,    // out: store diagonal projection of origin
						nbts,   // out: # backtracks
						row,    // start in this rectangle row
						col,    // start in this rectangle column
						rnd);   // random gen, to choose among equal paths
				rnd.init(reseed); // debug/release pseudo-randoms in lock step
			}
			if(ret) {
				btncand_[cural_].fate = BT_CAND_FATE_SUCCEEDED;
				break;
			} else {
				btncand_[cural_].fate = BT_CAND_FATE_FAILED;
			}
		}
		cural_++;
	} // while(cural_ < btncand_.size())
	if(cural_ == btncand_.size()) {
		assert(res.repOk());
		return false;
	}
	assert(!res.alres.empty());
	assert(res.repOk());
	if(!fw_) {
		// All edits are currently w/r/t upstream end; if read aligned
		// to Crick strand, we need to invert them so that they're
		// w/r/t the read's 5' end instead.
		res.alres.invertEdits();
	}
	cural_++;
	assert(res.repOk());
	return true;
}

