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

/**
 * aligner_sw_sse.cpp
 *
 * Versions of key alignment functions that use vector instructions to
 * accelerate dynamic programming.  Based chiefly on the striped Smith-Waterman
 * paper and implementation by Michael Farrar.  See:
 *
 * Farrar M. Striped Smith-Waterman speeds database searches six times over
 * other SIMD implementations. Bioinformatics. 2007 Jan 15;23(2):156-61.
 * http://sites.google.com/site/farrarmichael/smith-waterman
 *
 * While the paper describes an implementation of Smith-Waterman, we extend it
 * do end-to-end read alignment as well as local alignment.  The change
 * required for this is minor: we simply let vmax be the maximum element in the
 * score domain rather than the minimum.
 *
 * The vectorized dynamic programming implementation lacks some features that
 * make it hard to adapt to solving the entire dynamic-programming alignment
 * problem.  For instance:
 *
 * - It doesn't respect gap barriers on either end of the read
 * - It just gives a maximum; not enough information to backtrace without
 *   redoing some alignment
 * - It's a little difficult to handle st_ and en_, especially st_.
 * - The query profile mechanism makes handling of ambiguous reference bases a
 *   little tricky (16 cols in query profile lookup table instead of 5)
 *
 * Given the drawbacks, it is tempting to use SSE dynamic programming as a
 * filter rather than as an aligner per se.  Here are a few ideas for how it
 * can be extended to handle more of the alignment problem:
 *
 * - Save calculated scores to a big array as we go.  We return to this array
 *   to find and backtrace from good solutions.
 */

#include <limits>
#include "aligner_sw.h"

static constexpr size_t EEI16_NWORDS_PER_REG  = NBYTES_PER_REG/2;
static constexpr size_t EEI16_NBITS_PER_WORD  = 16;
static constexpr size_t EEI16_NBYTES_PER_WORD = 2;

// In 16-bit end-to-end mode, we have the option of using signed saturated
// arithmetic.  Because we have signed arithmetic, there's no need to add/subtract
// bias when building an applying the query profile.  The lowest value we can
// use is 0x8000, and the greatest is 0x7fff.

typedef int16_t EEI16_TCScore;

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
void SwAligner::buildQueryProfileEnd2EndSseI16(bool fw) {
	bool& done = fw ? sseI16fwBuilt_ : sseI16rcBuilt_;
	if(done) {
		return;
	}
	done = true;
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	const size_t len = rd->length();
	const size_t seglen = (len + (EEI16_NWORDS_PER_REG-1)) / EEI16_NWORDS_PER_REG;
	// How many SSERegI's are needed
	size_t nsses =
		64 +                    // slack bytes, for alignment?
		(seglen * ALPHA_SIZE)   // query profile data
		* 2;                    // & gap barrier data
	assert_gt(nsses, 0);
	SSEData& d = fw ? sseI16fw_ : sseI16rc_;
	d.profbuf_.resizeNoCopy(nsses);
	assert(!d.profbuf_.empty());
	d.maxPen_      = d.maxBonus_ = 0;
	d.lastIter_    = d.lastWord_ = 0;
	d.qprofStride_ = d.gbarStride_ = 2;
	d.bias_ = 0; // no bias when words are signed
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			int16_t *qprofWords =
				reinterpret_cast<int16_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2));
			int16_t *gbarWords =
				reinterpret_cast<int16_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2) + 1);
			// For each sub-word (byte) ...
			for(size_t k = 0; k < EEI16_NWORDS_PER_REG; k++) {
				int sc = 0;
				*gbarWords = 0;
				if(j < len) {
					int readc = (*rd)[j];
					int readq = (*qu)[j];
					sc = sc_->score(readc, (int)(1 << refc), readq - 33);
					size_t j_from_end = len - j - 1;
					if(j < (size_t)sc_->gapbar ||
					   j_from_end < (size_t)sc_->gapbar)
					{
						// Inside the gap barrier
						*gbarWords = 0x8000; // add this twice
					}
				}
				if(refc == 0 && j == len-1) {
					// Remember which 128-bit word and which smaller word has
					// the final row
					d.lastIter_ = i;
					d.lastWord_ = k;
				}
				if(sc < 0) {
					if((size_t)(-sc) > d.maxPen_) {
						d.maxPen_ = (size_t)(-sc);
					}
				} else {
					if((size_t)sc > d.maxBonus_) {
						d.maxBonus_ = (size_t)sc;
					}
				}
				*qprofWords = (int16_t)sc;
				gbarWords++;
				qprofWords++;
				j += seglen; // update offset into query
			}
		}
	}
}

#ifndef NDEBUG
/**
 * Return true iff the cell has sane E/F/H values w/r/t its predecessors.
 */
static bool cellOkEnd2EndI16(
	SSEData& d,
	size_t row,
	size_t col,
	int refc,
	int readc,
	int readq,
	const Scoring& sc)     // scoring scheme
{
	EEI16_TCScore floorsc = 0x8000;
	EEI16_TCScore ceilsc = MAX_I64;
	TAlScore offsetsc = -0x7fff;
	TAlScore sc_h_cur = (TAlScore)d.mat_.helt(row, col);
	TAlScore sc_e_cur = (TAlScore)d.mat_.eelt(row, col);
	TAlScore sc_f_cur = (TAlScore)d.mat_.felt(row, col);
	if(sc_h_cur > floorsc) {
		sc_h_cur += offsetsc;
	}
	if(sc_e_cur > floorsc) {
		sc_e_cur += offsetsc;
	}
	if(sc_f_cur > floorsc) {
		sc_f_cur += offsetsc;
	}
	bool gapsAllowed = true;
	size_t rowFromEnd = d.mat_.nrow() - row - 1;
	if(row < (size_t)sc.gapbar || rowFromEnd < (size_t)sc.gapbar) {
		gapsAllowed = false;
	}
	bool e_left_trans = false, h_left_trans = false;
	bool f_up_trans   = false, h_up_trans = false;
	bool h_diag_trans = false;
	if(gapsAllowed) {
		TAlScore sc_h_left = floorsc;
		TAlScore sc_e_left = floorsc;
		TAlScore sc_h_up   = floorsc;
		TAlScore sc_f_up   = floorsc;
		if(col > 0 && sc_e_cur > floorsc && sc_e_cur <= ceilsc) {
			sc_h_left = d.mat_.helt(row, col-1) + offsetsc;
			sc_e_left = d.mat_.eelt(row, col-1) + offsetsc;
			e_left_trans = (sc_e_left > floorsc && sc_e_cur == sc_e_left - sc.readGapExtend());
			h_left_trans = (sc_h_left > floorsc && sc_e_cur == sc_h_left - sc.readGapOpen());
			assert(e_left_trans || h_left_trans);
		}
		if(row > 0 && sc_f_cur > floorsc && sc_f_cur <= ceilsc) {
			sc_h_up = d.mat_.helt(row-1, col) + offsetsc;
			sc_f_up = d.mat_.felt(row-1, col) + offsetsc;
			f_up_trans = (sc_f_up > floorsc && sc_f_cur == sc_f_up - sc.refGapExtend());
			h_up_trans = (sc_h_up > floorsc && sc_f_cur == sc_h_up - sc.refGapOpen());
			assert(f_up_trans || h_up_trans);
		}
	} else {
		assert_geq(floorsc, sc_e_cur);
		assert_geq(floorsc, sc_f_cur);
	}
	if(col > 0 && row > 0 && sc_h_cur > floorsc && sc_h_cur <= ceilsc) {
		TAlScore sc_h_upleft = d.mat_.helt(row-1, col-1) + offsetsc;
		TAlScore sc_diag = sc.score(readc, (int)refc, readq - 33);
		h_diag_trans = sc_h_cur == sc_h_upleft + sc_diag;
	}
	assert(
		sc_h_cur <= floorsc ||
		e_left_trans ||
		h_left_trans ||
		f_up_trans   ||
		h_up_trans   ||
		h_diag_trans ||
		sc_h_cur > ceilsc ||
		row == 0 ||
		col == 0);
	return true;
}
#endif /*ndef NDEBUG*/

#ifdef NDEBUG

#define assert_all_lt(x, y)

#else

#define assert_all_lt(x, y) { \
	SSERegI tmp = sse_cmplt_epi16(x, y); \
	assert_eq(SSE_MASK_ALL, sse_movemask_epi8(tmp)); \
}

#endif

/**
 * Given a filled-in DP table, populate the btncand_ list with candidate cells
 * that might be at the ends of valid alignments.  No need to do this unless
 * the maximum score returned by the align*() func is >= the minimum.
 *
 * Only cells that are exhaustively scored are candidates.  Those are the
 * cells inside the shape made of o's in this:
 *
 *  |-maxgaps-|
 *  *********************************    -
 *   ********************************    |
 *    *******************************    |
 *     ******************************    |
 *      *****************************    |
 *       **************************** read len
 *        ***************************    |
 *         **************************    |
 *          *************************    |
 *           ************************    |
 *            ***********oooooooooooo    -
 *            |-maxgaps-|
 *  |-readlen-|
 *  |-------skip--------|
 *
 * And it's possible for the shape to be truncated on the left and right sides.
 *
 * 
 */

#define sse_anygt_epi16(val1,val2,outval) { \
	SSERegI s = sse_cmpgt_epi16(val1, val2); \
        outval = (sse_movemask_epi8(s) != 0); }

/**
 * Solve the current alignment problem using SSE instructions that operate on 8
 * signed 16-bit values packed into a single 128-bit register.
 */
template<typename TIdxSize>
inline EEI16_TCScore EEI16_alignNucleotides(const SSERegI profbuf[],
					const char   rf[], const TIdxSize rfd,
					SSERegI pmat[],
					const TIdxSize iter, const size_t colstride, const size_t lastWordIdx,
					const TAlScore minsc, const size_t nrow,
					DpBtCandidate btncand[], TIdxSize& btnfilled_,
					const int16_t refGapOpen, const int16_t refGapExtend, const int16_t readGapOpen, const int16_t readGapExtend) {

	// Many thanks to Michael Farrar for releasing his striped Smith-Waterman
	// implementation:
	//
	//  http://sites.google.com/site/farrarmichael/smith-waterman
	//
	// Much of the implmentation below is adapted from Michael's code.

	// Set all elts to reference gap open penalty
	SSERegI rfgapo   = sse_setzero_siall();
	SSERegI rfgape   = sse_setzero_siall();
	SSERegI rdgapo   = sse_setzero_siall();
	SSERegI rdgape   = sse_setzero_siall();
	SSERegI vlo      = sse_setzero_siall();
	SSERegI vhilsw   = sse_setzero_siall();
	SSERegI vlolsw   = sse_setzero_siall();
	SSERegI ve       = sse_setzero_siall();
	SSERegI vf       = sse_setzero_siall();
	SSERegI vh       = sse_setzero_siall();
	SSERegI vtmp     = sse_setzero_siall();

	SSERegI vs0     = sse_setzero_siall();
	SSERegI vs1     = sse_setzero_siall();

	assert_gt(refGapOpen, 0);
	assert_leq(refGapOpen, MAX_I16);
	sse_fill_i16(refGapOpen, rfgapo);
	
	// Set all elts to reference gap extension penalty
	assert_gt(refGapExtend, 0);
	assert_leq(refGapExtend, MAX_I16);
	assert_leq(refGapExtend, refGapOpen);
	sse_fill_i16(refGapExtend, rfgape);

	// Set all elts to read gap open penalty
	assert_gt(readGapOpen, 0);
	assert_leq(readGapOpen, MAX_I16);
	sse_fill_i16(readGapOpen, rdgapo);
	
	// Set all elts to read gap extension penalty
	assert_gt(readGapExtend, 0);
	assert_leq(readGapExtend, MAX_I16);
	assert_leq(readGapExtend, readGapOpen);
	sse_fill_i16(readGapExtend, rdgape);

	// Set all elts to 0x8000 (min value for signed 16-bit)
	sse_fill_i16(0x8000, vlo);

	// vlolsw: topmost (least sig) word set to 0x8000, all other words=0
	sse_set_low_i16(0x8000, vlolsw);
	
	// vhilsw: topmost (least sig) word set to 0x7fff, all other words=0
	sse_set_low_i16(0x7fff, vhilsw);
	
	// Points to a long vector of SSERegI where each element is a block of
	// contiguous cells in the E, F or H matrix.  If the index % 3 == 0, then
	// the block of cells is from the E matrix.  If index % 3 == 1, they're
	// from the F matrix.  If index % 3 == 2, then they're from the H matrix.
	// Blocks of cells are organized in the same interleaved manner as they are
	// calculated by the Farrar algorithm.
	// const SSERegI *pvScore; // points into the query profile

	assert_eq(ROWSTRIDE, colstride / iter);

	// Initialize the H and E vectors in the first matrix column
	{
	  SSERegI *pvHTmp = pmat + SSEMatrix::TMP;
	  SSERegI *pvETmp = pmat + SSEMatrix::E;
	
	  for(size_t i = 0; i < iter; i++) {
		sse_store_siall(pvETmp, vlo);
		// Could initialize Hs to high or low.  If high, cells in the lower
		// triangle will have somewhat more legitiate scores, but still won't
		// be exhaustively scored.
		sse_store_siall(pvHTmp, vlo);
		pvETmp += ROWSTRIDE;
		pvHTmp += ROWSTRIDE;
	  }
	}

	// These are swapped just before the innermost loop
	SSERegI *pvHStore = pmat + SSEMatrix::H;
	SSERegI *pvHLoad  = pmat + SSEMatrix::TMP;
	SSERegI *pvELoad  = pmat + SSEMatrix::E;
	SSERegI *pvEStore = pmat + colstride + SSEMatrix::E;
	SSERegI *pvFStore = pmat + SSEMatrix::F;
	SSERegI *pvFTmp   = NULL;
	
	// Maximum score in final row
	EEI16_TCScore lrmax = MIN_I16;

	// keep a local copy
	TIdxSize btnfilled = 0;

	// Fill in the table as usual but instead of using the same gap-penalty
	// vector for each iteration of the inner loop, load words out of a
	// pre-calculated gap vector parallel to the query profile.  The pre-
	// calculated gap vectors enforce the gap barrier constraint by making it
	// infinitely costly to introduce a gap in barrier rows.
	//
	// AND use a separate loop to fill in the first row of the table, enforcing
	// the st_ constraints in the process.  This is awkward because it
	// separates the processing of the first row from the others and might make
	// it difficult to use the first-row results in the next row, but it might
	// be the simplest and least disruptive way to deal with the st_ constraint.

	for(TIdxSize i = 0; i < rfd; i++) {
		assert(pvFStore == (pmat + i*colstride + SSEMatrix::F));
		assert(pvHStore == (pmat + i*colstride + SSEMatrix::H));
		
		// Fetch the appropriate query profile.  Note that elements of rf must
		// be numbers, not masks.
		size_t off = (size_t)firsts5[ rf[i] ] * iter * 2;
		// points into the query profile
		const SSERegI *pvScore = profbuf + off; // even elts = query profile, odd = gap barrier
		
		// Set all cells to low value
		vf = sse_cmpeq_epi16(vf, vf);
		vf = sse_slli_epi16(vf, EEI16_NBITS_PER_WORD-1);
		vf = sse_or_siall(vf, vlolsw);
		
		// Load H vector from the final row of the previous column
		vh = sse_load_siall(pvHLoad + colstride - ROWSTRIDE);
		// Shift 2 bytes down so that topmost (least sig) cell gets 0
		vh = sse_slli_siall(vh, EEI16_NBYTES_PER_WORD);
		// Fill topmost (least sig) cell with high value
		vh = sse_or_siall(vh, vhilsw);
		
		// For each character in the reference text:
		for(TIdxSize j = 0; j < iter; j++) {
			vs0 = sse_load_siall(pvScore);
                        pvScore++;
			vs1 = sse_load_siall(pvScore);
                        pvScore++;
			// Load cells from E, calculated previously
			ve = sse_load_siall(pvELoad);
			assert_all_lt(ve, vhi);
			pvELoad += ROWSTRIDE;
			
			// Store cells in F, calculated previously
			vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
			vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Factor in query profile (matches and mismatches)
			vh = sse_adds_epi16(vh, vs0);
			
			// Update H, factoring in E and F
			vh = sse_max_epi16(vh, ve);
			vh = sse_max_epi16(vh, vf);
			
			// Save the new vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update vE value
			vtmp = vh;
			vh = sse_subs_epi16(vh, rdgapo);
			vh = sse_adds_epi16(vh, vs1); // veto some read gap opens
			vh = sse_adds_epi16(vh, vs1); // veto some read gap opens
			ve = sse_subs_epi16(ve, rdgape);
			ve = sse_max_epi16(ve, vh);

			assert_all_lt(ve, vhi);
			
			// Load the next h value
			vh = sse_load_siall(pvHLoad);
			pvHLoad += ROWSTRIDE;
			
			// Save E values
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			
			// Update vf value
			vtmp = sse_subs_epi16(vtmp, rfgapo);
			vf = sse_subs_epi16(vf, rfgape);
			assert_all_lt(vf, vhi);
			vf = sse_max_epi16(vf, vtmp);
		}
		pvScore -= (iter*2) - 1; // reset veto vector
	        vs1 = sse_load_siall(pvScore);

		// pvHStore, pvELoad, pvEStore have all rolled over to the next column
		pvFTmp = pvFStore;
		pvFStore -= colstride; // reset to start of column
		vtmp = sse_load_siall(pvFStore);
		
		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		vf = sse_slli_siall(vf, EEI16_NBYTES_PER_WORD);
		vf = sse_or_siall(vf, vlolsw);
		
		vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
		vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
		vf = sse_max_epi16(vtmp, vf);
		bool anygt;
		sse_anygt_epi16(vf,vtmp,anygt);
		// Load after computing cmp, so the result is ready by the time it is tested in while
		pvHStore -= colstride; // reset to start of column
		vh = sse_load_siall(pvHStore);
		pvHLoad = pvHStore;    // new pvHLoad = pvHStore
		
		pvEStore -= colstride; // reset to start of column
		ve = sse_load_siall(pvEStore);
		
		// If any element of vtmp is greater than H - gap-open...
		TIdxSize j = 0;
		while(anygt) {
			// Store this vf
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Update vh w/r/t new vf
			vh = sse_max_epi16(vh, vf);
			
			// Save vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update E in case it can be improved using our new vh
#if 0
#else
			vh = sse_subs_epi16(vh, rdgapo);
			vh = sse_adds_epi16(vh, vs1); // veto some read gap opens
			vh = sse_adds_epi16(vh, vs1); // veto some read gap opens
			ve = sse_max_epi16(ve, vh);
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
#endif
			pvScore += 2;
			
			assert_lt(j, iter);
			if(++j == iter) {
				pvScore -= iter*2;
				j = 0;
				pvFStore -= colstride;
				pvHStore -= colstride;
#if 0
#else
				pvEStore -= colstride;
#endif
				vf = sse_slli_siall(vf, EEI16_NBYTES_PER_WORD);
				vf = sse_or_siall(vf, vlolsw);
			}
			vs1 = sse_load_siall(pvScore);
			vtmp = sse_load_siall(pvFStore);   // load next vf ASAP
			
			// Update F with another gap extension
			vf = sse_subs_epi16(vf, rfgape);
			vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
			vf = sse_adds_epi16(vf, vs1); // veto some ref gap extensions
			vf = sse_max_epi16(vtmp, vf);
			sse_anygt_epi16(vf,vtmp,anygt);

			// Load after computing cmp, so the result is afailable by the time it is tested
			vh = sse_load_siall(pvHStore);     // load next vh ASAP
			ve = sse_load_siall(pvEStore);     // load next vh ASAP
		}
		
		// Note: we may not want to extract from the final row
		EEI16_TCScore lr = ((EEI16_TCScore*)(pvHLoad))[lastWordIdx];
		TAlScore sc = (TAlScore)(lr - 0x7fff);
		if(lr > lrmax) {
			lrmax = lr;
		}
		if(sc >= minsc) {
			// Yes, this is legit
			btncand[btnfilled].init(nrow-1, i, lr);
			btnfilled++;
		}

		// pvELoad and pvHLoad are already where they need to be
		
		// Adjust the load and store vectors here.  
		pvHStore = pvHLoad + colstride;
		pvEStore = pvELoad + colstride;
		pvFStore = pvFTmp;
	}


	btnfilled_ = btnfilled;  // pass it out
	return lrmax;
}

/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.
 */
bool SwAligner::alignEnd2EndSseI16(
	TAlScore& best)    // best alignment score observed in DP matrix
{
	constexpr bool debug = false;

	assert(initedRef() && initedRead());
	assert_eq(STATE_INITED, state_);

	assert_leq(rdf_, rd_->length());
	assert_leq(rdf_, qu_->length());
	assert_lt(rfi_, rff_);
	assert_lt(rdi_, rdf_);
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = (size_t)rfi_; i < (size_t)rff_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
#ifdef ENABLE_SSE_METRICS
	SSEMetrics& met = sseMet_;
	if(!debug) met.dp++;
#endif
	buildQueryProfileEnd2EndSseI16(fw_);
	assert(!d.profbuf_.empty());

	assert_eq(0, d.maxBonus_);
	const size_t iter =
		(dpRows() + (EEI16_NWORDS_PER_REG-1)) / EEI16_NWORDS_PER_REG; // iter = segLen
        const size_t lastWordIdx = EEI16_NWORDS_PER_REG*(d.lastIter_*ROWSTRIDE)+d.lastWord_;


	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_I16);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_I16);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_I16);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_I16);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());

	assert_gt(sc_->gapbar, 0);

	colstop_ = rff_ - 1;
	lastsolcol_ = 0;

	// should never get in here, but just in case
	if (rfi_>=rff_) {
		btncand_.clear();
#ifdef ENABLE_SSE_METRICS
		if(!debug) met.dpfail++;
#endif
		best = MIN_I64;
		return false;
	}

	d.mat_.init(dpRows(), rff_-rfi_, EEI16_NWORDS_PER_REG);

	assert_leq(iter,      (size_t)MAX_U16);
	assert_leq(rff_-rfi_, (size_t)MAX_U16);
	uint16_t btnfilled = 0;
	btncand_.resizeNoCopy(rff_-rfi_); // cannot be bigger that this

	const EEI16_TCScore lrmax = EEI16_alignNucleotides<uint16_t>(d.profbuf_.ptr(), rf_+rfi_, rff_-rfi_,
					d.mat_.ptr(),
                                        iter, d.mat_.colstride(), lastWordIdx,
					minsc_, dpRows(),
					btncand_.ptr(), btnfilled,
					sc_->refGapOpen(), sc_->refGapExtend(), sc_->readGapOpen(), sc_->readGapExtend());
	btncand_.trim(btnfilled);
	state_ = STATE_ALIGNED;
	cural_ = 0;

#ifdef ENABLE_SSE_METRICS
	// Update metrics
	if(!debug) {
		size_t ninner = (rff_ - rfi_) * iter;
		met.col   += (rff_ - rfi_);             // DP columns
		met.cell  += (ninner * EEI16_NWORDS_PER_REG); // DP cells
		met.inner += ninner;                    // DP inner loop iters
		met.fixup += 0; // deprecated nfixup;                    // DP fixup loop iters
	}
#endif
	
	// Did we find a solution?
	const TAlScore score = (TAlScore)(lrmax - 0xff);
	if(score < minsc_) {
		// no
#ifdef ENABLE_SSE_METRICS
		if(!debug) met.dpfail++;
#endif
		best = score;
		return false;
	}
	
	// Could we have saturated?
	if(lrmax == MIN_I16) {
		// yes
#ifdef ENABLE_SSE_METRICS
		if(!debug) met.dpsat++;
#endif
		best = MIN_I64;
		return false;
	}
	
	// Return largest score
#ifdef ENABLE_SSE_METRICS
	if(!debug) met.dpsucc++;
#endif

	if (btnfilled>0) {
#ifdef ENABLE_SSE_METRICS
		met.gathsol += btnfilled;
#endif
		d.mat_.initMasks();
		btncand_.sort();
	}

	best = score;
	return !btncand_.empty();
}

#define MOVE_VEC_PTR_UP(vec, rowvec, rowelt) { \
	if(rowvec == 0) { \
		rowvec += d.mat_.nvecrow_; \
		vec += d.mat_.colstride_; \
		rowelt--; \
	} \
	rowvec--; \
	vec -= ROWSTRIDE; \
}

#define MOVE_VEC_PTR_LEFT(vec, rowvec, rowelt) { vec -= d.mat_.colstride_; }

#define MOVE_VEC_PTR_UPLEFT(vec, rowvec, rowelt) { \
 	MOVE_VEC_PTR_UP(vec, rowvec, rowelt); \
 	MOVE_VEC_PTR_LEFT(vec, rowvec, rowelt); \
}

#define MOVE_ALL_LEFT() { \
	MOVE_VEC_PTR_LEFT(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_LEFT(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_LEFT(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_LEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define MOVE_ALL_UP() { \
	MOVE_VEC_PTR_UP(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_UP(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_UP(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_UP(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define MOVE_ALL_UPLEFT() { \
	MOVE_VEC_PTR_UPLEFT(cur_vec, rowvec, rowelt); \
	MOVE_VEC_PTR_UPLEFT(left_vec, left_rowvec, left_rowelt); \
	MOVE_VEC_PTR_UPLEFT(up_vec, up_rowvec, up_rowelt); \
	MOVE_VEC_PTR_UPLEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

#define NEW_ROW_COL(row, col) { \
	rowelt = row / d.mat_.nvecrow_; \
	rowvec = row % d.mat_.nvecrow_; \
	eltvec = (col * d.mat_.colstride_) + (rowvec * ROWSTRIDE); \
	cur_vec = d.mat_.matbuf_.ptr() + eltvec; \
	left_vec = cur_vec; \
	left_rowelt = rowelt; \
	left_rowvec = rowvec; \
	MOVE_VEC_PTR_LEFT(left_vec, left_rowvec, left_rowelt); \
	up_vec = cur_vec; \
	up_rowelt = rowelt; \
	up_rowvec = rowvec; \
	MOVE_VEC_PTR_UP(up_vec, up_rowvec, up_rowelt); \
	upleft_vec = up_vec; \
	upleft_rowelt = up_rowelt; \
	upleft_rowvec = up_rowvec; \
	MOVE_VEC_PTR_LEFT(upleft_vec, upleft_rowvec, upleft_rowelt); \
}

/**
 * Given the dynamic programming table and a cell, trace backwards from the
 * cell and install the edits and score/penalty in the appropriate fields
 * of res.  The RandomSource is used to break ties among equally good ways
 * of tracing back.
 *
 * Whenever we enter a cell, we check whether the read/ref coordinates of
 * that cell correspond to a cell we traversed constructing a previous
 * alignment.  If so, we backtrack to the last decision point, mask out the
 * path that led to the previously observed cell, and continue along a
 * different path; or, if there are no more paths to try, we give up.
 *
 * If an alignment is found, 'off' is set to the alignment's upstream-most
 * reference character's offset into the chromosome and true is returned.
 * Otherwise, false is returned.
 */
bool SwAligner::backtraceNucleotidesEnd2EndSseI16(
	TAlScore       escore, // in: expected score
	SwResult&      res,    // out: store results (edits and scores) here
	size_t&        off,    // out: store diagonal projection of origin
	size_t&        nbts,   // out: # backtracks
	size_t         row,    // start in this row
	size_t         col,    // start in this column
	RandomSource&  rnd)    // random gen, to choose among equal paths
{
	assert_lt(row, dpRows());
	assert_lt(col, (size_t)(rff_ - rfi_));
	SSEData& d = fw_ ? sseI16fw_ : sseI16rc_;
#ifdef ENABLE_SSE_METRICS
	SSEMetrics& met = sseMet_;
	met.bt++;
#endif
	assert(!d.profbuf_.empty());
	assert_lt(row, rd_->length());
	btnstack_.clear(); // empty the backtrack stack
	btcells_.clear();  // empty the cells-so-far list
	AlnScore score; score.score_ = 0;
	score.gaps_ = score.ns_ = 0;
	size_t origCol = col;
	size_t gaps = 0, readGaps = 0, refGaps = 0;
	res.alres.reset();
	EList<Edit>& ned = res.alres.ned();
	assert(ned.empty());
	assert_gt(dpRows(), row);
	size_t trimEnd = dpRows() - row - 1; 
	size_t trimBeg = 0;
	size_t ct = SSEMatrix::H; // cell type
	// Row and col in terms of where they fall in the SSE vector matrix
	size_t rowelt, rowvec, eltvec;
	size_t left_rowelt, up_rowelt, upleft_rowelt;
	size_t left_rowvec, up_rowvec, upleft_rowvec;
	SSERegI *cur_vec, *left_vec, *up_vec, *upleft_vec;
	NEW_ROW_COL(row, col);
	while((int)row >= 0) {
#ifdef ENABLE_SSE_METRICS
		met.btcell++;
#endif
		nbts++;
		int readc = (*rd_)[rdi_ + row];
		int refm  = (int)rf_[rfi_ + col];
		int readq = (*qu_)[row];
		assert_leq(col, origCol);
		// Get score in this cell
		bool empty = false, reportedThru, canMoveThru, branch = false;
		int cur = SSEMatrix::H;
		if(!d.mat_.reset_[row]) {
			d.mat_.resetRow(row);
		}
		reportedThru = d.mat_.reportedThrough(row, col);
		canMoveThru = true;
		if(reportedThru) {
			canMoveThru = false;
		} else {
			empty = false;
			if(row > 0) {
				assert_gt(row, 0);
				size_t rowFromEnd = d.mat_.nrow() - row - 1;
				bool gapsAllowed = true;
				if(row < (size_t)sc_->gapbar ||
				   rowFromEnd < (size_t)sc_->gapbar)
				{
					gapsAllowed = false;
				}
				const TAlScore floorsc = MIN_I64;
				const int offsetsc = -0x7fff;
				// Move to beginning of column/row
				if(ct == SSEMatrix::E) { // AKA rdgap
					assert_gt(col, 0);
					TAlScore sc_cur = ((EEI16_TCScore*)(cur_vec + SSEMatrix::E))[rowelt] + offsetsc;
					assert(gapsAllowed);
					// Currently in the E matrix; incoming transition must come from the
					// left.  It's either a gap open from the H matrix or a gap extend from
					// the E matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell to the left
					TAlScore sc_h_left = ((EEI16_TCScore*)(left_vec + SSEMatrix::H))[left_rowelt] + offsetsc;
					if(sc_h_left > floorsc && sc_h_left - sc_->readGapOpen() == sc_cur) {
						mask |= (1 << 0);
					}
					// Get E score of cell to the left
					TAlScore sc_e_left = ((EEI16_TCScore*)(left_vec + SSEMatrix::E))[left_rowelt] + offsetsc;
					if(sc_e_left > floorsc && sc_e_left - sc_->readGapExtend() == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isEMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 8) & 3;
					}
					if(mask == 3) {
#if 1
						// Pick H -> E cell
						cur = SW_BT_OALL_READ_OPEN;
						d.mat_.eMaskSet(row, col, 2); // might choose E later
#else
						if(rnd.nextU2()) {
							// Pick H -> E cell
							cur = SW_BT_OALL_READ_OPEN;
							d.mat_.eMaskSet(row, col, 2); // might choose E later
						} else {
							// Pick E -> E cell
							cur = SW_BT_RDGAP_EXTEND;
							d.mat_.eMaskSet(row, col, 1); // might choose H later
						}
#endif
						branch = true;
					} else if(mask == 2) {
						// I chose the E cell
						cur = SW_BT_RDGAP_EXTEND;
						d.mat_.eMaskSet(row, col, 0); // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_READ_OPEN;
						d.mat_.eMaskSet(row, col, 0); // done
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
					assert(!empty || !canMoveThru);
				} else if(ct == SSEMatrix::F) { // AKA rfgap
					assert_gt(row, 0);
					assert(gapsAllowed);
					TAlScore sc_h_up = ((EEI16_TCScore*)(up_vec  + SSEMatrix::H))[up_rowelt] + offsetsc;
					TAlScore sc_f_up = ((EEI16_TCScore*)(up_vec  + SSEMatrix::F))[up_rowelt] + offsetsc;
					TAlScore sc_cur  = ((EEI16_TCScore*)(cur_vec + SSEMatrix::F))[rowelt] + offsetsc;
					// Currently in the F matrix; incoming transition must come from above.
					// It's either a gap open from the H matrix or a gap extend from the F
					// matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell above
					if(sc_h_up > floorsc && sc_h_up - sc_->refGapOpen() == sc_cur) {
						mask |= (1 << 0);
					}
					// Get F score of cell above
					if(sc_f_up > floorsc && sc_f_up - sc_->refGapExtend() == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isFMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 11) & 3;
					}
					if(mask == 3) {
#if 1
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						d.mat_.fMaskSet(row, col, 2); // might choose E later
#else
						if(rnd.nextU2()) {
							// I chose the H cell
							cur = SW_BT_OALL_REF_OPEN;
							d.mat_.fMaskSet(row, col, 2); // might choose E later
						} else {
							// I chose the F cell
							cur = SW_BT_RFGAP_EXTEND;
							d.mat_.fMaskSet(row, col, 1); // might choose E later
						}
#endif
						branch = true;
					} else if(mask == 2) {
						// I chose the F cell
						cur = SW_BT_RFGAP_EXTEND;
						d.mat_.fMaskSet(row, col, 0); // done
					} else if(mask == 1) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						d.mat_.fMaskSet(row, col, 0); // done
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
					assert(!empty || !canMoveThru);
				} else {
					assert_eq(SSEMatrix::H, ct);
					TAlScore sc_cur      = ((EEI16_TCScore*)(cur_vec + SSEMatrix::H))[rowelt]    + offsetsc;
					TAlScore sc_f_up     = ((EEI16_TCScore*)(up_vec  + SSEMatrix::F))[up_rowelt] + offsetsc;
					TAlScore sc_h_up     = ((EEI16_TCScore*)(up_vec  + SSEMatrix::H))[up_rowelt] + offsetsc;
					TAlScore sc_h_left   = col > 0 ? (((EEI16_TCScore*)(left_vec   + SSEMatrix::H))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_e_left   = col > 0 ? (((EEI16_TCScore*)(left_vec   + SSEMatrix::E))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_h_upleft = col > 0 ? (((EEI16_TCScore*)(upleft_vec + SSEMatrix::H))[upleft_rowelt] + offsetsc) : floorsc;
					TAlScore sc_diag     = sc_->score(readc, refm, readq - 33);
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					if(gapsAllowed) {
						if(sc_h_up     > floorsc && sc_cur == sc_h_up   - sc_->refGapOpen()) {
							mask |= (1 << 0);
						}
						if(sc_h_left   > floorsc && sc_cur == sc_h_left - sc_->readGapOpen()) {
							mask |= (1 << 1);
						}
						if(sc_f_up     > floorsc && sc_cur == sc_f_up   - sc_->refGapExtend()) {
							mask |= (1 << 2);
						}
						if(sc_e_left   > floorsc && sc_cur == sc_e_left - sc_->readGapExtend()) {
							mask |= (1 << 3);
						}
					}
					if(sc_h_upleft > floorsc && sc_cur == sc_h_upleft + sc_diag) {
						mask |= (1 << 4);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isHMaskSet(row, col)) {
						mask = (d.mat_.masks_[row][col] >> 2) & 31;
					}
					assert(gapsAllowed || mask == (1 << 4) || mask == 0);
					int opts = alts5[mask];
					int select = -1;
					if(opts == 1) {
						select = firsts5[mask];
						assert_geq(mask, 0);
						d.mat_.hMaskSet(row, col, 0);
					} else if(opts > 1) {
#if 1
						if(       (mask & 16) != 0) {
							select = 4; // H diag
						} else if((mask & 1) != 0) {
							select = 0; // H up
						} else if((mask & 4) != 0) {
							select = 2; // F up
						} else if((mask & 2) != 0) {
							select = 1; // H left
						} else if((mask & 8) != 0) {
							select = 3; // E left
						}
#else
						select = randFromMask(rnd, mask);
#endif
						assert_geq(mask, 0);
						mask &= ~(1 << select);
						assert(gapsAllowed || mask == (1 << 4) || mask == 0);
						d.mat_.hMaskSet(row, col, mask);
						branch = true;
					} else { /* No way to backtrack! */ }
					if(select != -1) {
						if(select == 4) {
							cur = SW_BT_OALL_DIAG;
						} else if(select == 0) {
							cur = SW_BT_OALL_REF_OPEN;
						} else if(select == 1) {
							cur = SW_BT_OALL_READ_OPEN;
						} else if(select == 2) {
							cur = SW_BT_RFGAP_EXTEND;
						} else {
							assert_eq(3, select)
							cur = SW_BT_RDGAP_EXTEND;
						}
					} else {
						empty = true;
						// It's empty, so the only question left is whether we should be
						// allowed in terimnate in this cell.  If it's got a valid score
						// then we *shouldn't* be allowed to terminate here because that
						// means it's part of a larger alignment that was already reported.
						canMoveThru = (origMask == 0);
					}
				}
				assert(!empty || !canMoveThru || ct == SSEMatrix::H);
			}
		}
		d.mat_.setReportedThrough(row, col);
		assert_eq(gaps, Edit::numGaps(ned));
		assert_leq(gaps, rdgap_ + rfgap_);
		// Cell was involved in a previously-reported alignment?
		if(!canMoveThru) {
			if(!btnstack_.empty()) {
				// Remove all the cells from list back to and including the
				// cell where the branch occurred
				btcells_.resize(btnstack_.back().celsz);
				// Pop record off the top of the stack
				ned.resize(btnstack_.back().nedsz);
				//aed.resize(btnstack_.back().aedsz);
				row      = btnstack_.back().row;
				col      = btnstack_.back().col;
				gaps     = btnstack_.back().gaps;
				readGaps = btnstack_.back().readGaps;
				refGaps  = btnstack_.back().refGaps;
				score    = btnstack_.back().score;
				ct       = btnstack_.back().ct;
				btnstack_.pop_back();
				assert(!sc_->monotone || score.score() >= escore);
				NEW_ROW_COL(row, col);
				continue;
			} else {
				// No branch points to revisit; just give up
				res.reset();
#ifdef ENABLE_SSE_METRICS
				met.btfail++; // DP backtraces failed
#endif
				return false;
			}
		}
		assert(!reportedThru);
		assert(!sc_->monotone || score.score() >= minsc_);
		if(empty || row == 0) {
			assert_eq(SSEMatrix::H, ct);
			btcells_.expand();
			btcells_.back().first = row;
			btcells_.back().second = col;
			// This cell is at the end of a legitimate alignment
			trimBeg = row;
			assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
			break;
		}
		if(branch) {
			// Add a frame to the backtrack stack
			btnstack_.expand();
			btnstack_.back().init(
				ned.size(),
				0,               // aed.size()
				btcells_.size(),
				row,
				col,
				gaps,
				readGaps,
				refGaps,
				score,
				(int)ct);
		}
		btcells_.expand();
		btcells_.back().first = row;
		btcells_.back().second = col;
		switch(cur) {
			// Move up and to the left.  If the reference nucleotide in the
			// source row mismatches the read nucleotide, penalize
			// it and add a nucleotide mismatch.
			case SW_BT_OALL_DIAG: {
				assert_gt(row, 0); assert_gt(col, 0);
				int readC = (*rd_)[row];
				int refNmask = (int)rf_[rfi_+col];
				assert_gt(refNmask, 0);
				int m = matchesEx(readC, refNmask);
				ct = SSEMatrix::H;
				if(m != 1) {
					Edit e(
						(int)row,
						mask2dna[refNmask],
						"ACGTN"[readC],
						EDIT_TYPE_MM);
					assert(e.repOk());
					assert(ned.empty() || ned.back().pos >= row);
					ned.push_back(e);
					int pen = QUAL2(row, col);
					score.score_ -= pen;
					assert(!sc_->monotone || score.score() >= escore);
				} else {
					// Reward a match
					int64_t bonus = sc_->match(30);
					score.score_ += bonus;
					assert(!sc_->monotone || score.score() >= escore);
				}
				if(m == -1) {
					score.ns_++;
				}
				row--; col--;
				MOVE_ALL_UPLEFT();
				assert(VALID_AL_SCORE(score));
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_OALL_REF_OPEN:
			{
				assert_gt(row, 0);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)(*rd_)[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrix::H;
				int pen = sc_->refGapOpen();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_UP();
				break;
			}
			// Move up.  Add an edit encoding the ref gap.
			case SW_BT_RFGAP_EXTEND:
			{
				assert_gt(row, 1);
				Edit e(
					(int)row,
					'-',
					"ACGTN"[(int)(*rd_)[row]],
					EDIT_TYPE_REF_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrix::F;
				int pen = sc_->refGapExtend();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; refGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_UP();
				break;
			}
			case SW_BT_OALL_READ_OPEN:
			{
				assert_gt(col, 0);
				Edit e(
					(int)row+1,
					mask2dna[(int)rf_[rfi_+col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrix::H;
				int pen = sc_->readGapOpen();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_LEFT();
				break;
			}
			case SW_BT_RDGAP_EXTEND:
			{
				assert_gt(col, 1);
				Edit e(
					(int)row+1,
					mask2dna[(int)rf_[rfi_+col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdf_-rdi_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrix::E;
				int pen = sc_->readGapExtend();
				score.score_ -= pen;
				assert(!sc_->monotone || score.score() >= minsc_);
				gaps++; readGaps++;
				assert_eq(gaps, Edit::numGaps(ned));
				assert_leq(gaps, rdgap_ + rfgap_);
				MOVE_ALL_LEFT();
				break;
			}
			default: throw 1;
		}
	} // while((int)row > 0)
	assert_eq(0, trimBeg);
	assert_eq(0, trimEnd);
	assert_geq(col, 0);
	assert_eq(SSEMatrix::H, ct);
	// The number of cells in the backtracs should equal the number of read
	// bases after trimming plus the number of gaps
	assert_eq(btcells_.size(), dpRows() - trimBeg - trimEnd + readGaps);
	// Check whether we went through a core diagonal and set 'reported' flag on
	// each cell
	bool overlappedCoreDiag = false;
	for(size_t i = 0; i < btcells_.size(); i++) {
		size_t rw = btcells_[i].first;
		size_t cl = btcells_[i].second;
		// Calculate the diagonal within the *trimmed* rectangle, i.e. the
		// rectangle we dealt with in align, gather and backtrack.
		int64_t diagi = cl - rw;
		// Now adjust to the diagonal within the *untrimmed* rectangle by
		// adding on the amount trimmed from the left.
		diagi += rect_->triml;
		if(diagi >= 0) {
			size_t diag = (size_t)diagi;
			if(diag >= rect_->corel && diag <= rect_->corer) {
				overlappedCoreDiag = true;
				break;
			}
		}
		assert(d.mat_.reportedThrough(rw, cl));
	}
	if(!overlappedCoreDiag) {
		// Must overlap a core diagonal.  Otherwise, we run the risk of
		// reporting an alignment that overlaps (and trumps) a higher-scoring
		// alignment that lies partially outside the dynamic programming
		// rectangle.
		res.reset();
#ifdef ENABLE_SSE_METRICS
		met.corerej++;
#endif
		return false;
	}
	int readC = (*rd_)[rdi_+row];      // get last char in read
	int refNmask = (int)rf_[rfi_+col]; // get last ref char ref involved in aln
	assert_gt(refNmask, 0);
	int m = matchesEx(readC, refNmask);
	if(m != 1) {
		Edit e((int)row, mask2dna[refNmask], "ACGTN"[readC], EDIT_TYPE_MM);
		assert(e.repOk());
		assert(ned.empty() || ned.back().pos >= row);
		ned.push_back(e);
		score.score_ -= QUAL2(row, col);
		assert_geq(score.score(), minsc_);
	} else {
		score.score_ += sc_->match(30);
	}
	if(m == -1) {
		score.ns_++;
	}
	if(score.ns_ > nceil_) {
		// Alignment has too many Ns in it!
		res.reset();
#ifdef ENABLE_SSE_METRICS
		met.nrej++;
#endif
		return false;
	}
	res.reverse();
	assert(Edit::repOk(ned, (*rd_)));
	assert_eq(score.score(), escore);
	assert_leq(gaps, rdgap_ + rfgap_);
	off = col;
	assert_lt(col + (size_t)rfi_, (size_t)rff_);
	score.gaps_ = gaps;
	score.edits_ = (int)ned.size();
	score.basesAligned_ = (int)(rdf_ - rdi_ - trimBeg - trimEnd - score.edits_);
	res.alres.setScore(score);
	res.alres.setShape(
		refidx_,                  // ref id
		off + rfi_ + rect_->refl, // 0-based ref offset
		reflen_,                  // reference length
		fw_,                      // aligned to Watson?
		rdf_ - rdi_,              // read length
		true,                     // pretrim soft?
		0,                        // pretrim 5' end
		0,                        // pretrim 3' end
		true,                     // alignment trim soft?
		fw_ ? trimBeg : trimEnd,  // alignment trim 5' end
		fw_ ? trimEnd : trimBeg); // alignment trim 3' end
	size_t refns = 0;
	for(size_t i = col; i <= origCol; i++) {
		if((int)rf_[rfi_+i] > 15) {
			refns++;
		}
	}
	res.alres.setRefNs(refns);
	assert(Edit::repOk(ned, (*rd_), true, trimBeg, trimEnd));
	assert(res.repOk());
#ifndef NDEBUG
	size_t gapsCheck = 0;
	for(size_t i = 0; i < ned.size(); i++) {
		if(ned[i].isGap()) gapsCheck++;
	}
	assert_eq(gaps, gapsCheck);
	BTDnaString refstr;
	for(size_t i = col; i <= origCol; i++) {
		refstr.append(firsts5[(int)rf_[rfi_+i]]);
	}
	BTDnaString editstr;
	Edit::toRef((*rd_), ned, editstr, true, trimBeg, trimEnd);
	if(refstr != editstr) {
		cerr << "Decoded nucleotides and edits don't match reference:" << endl;
		cerr << "           score: " << score.score()
		     << " (" << gaps << " gaps)" << endl;
		cerr << "           edits: ";
		Edit::print(cerr, ned);
		cerr << endl;
		cerr << "    decoded nucs: " << (*rd_) << endl;
		cerr << "     edited nucs: " << editstr << endl;
		cerr << "  reference nucs: " << refstr << endl;
		assert(0);
	}
#endif
#ifdef ENABLE_SSE_METRICS
	met.btsucc++; // DP backtraces succeeded
#endif
	return true;
}
