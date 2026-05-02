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
#include <bit>

#include "mask.h"

#ifndef SWSSE_INLINE_ONLY

#include "aligner_sw.h"

#else

#include "aligner_types.h"
#include "sse_wrap.h"
#include "aligner_swsse.h"
#include "aligner_sw_nuc.h"

#endif

// In end-to-end mode, we start high (255) and go low (0).  Factoring in
// a query profile involves unsigned saturating subtraction, so all the
// query profile elements should be expressed as a positive penalty rather
// than a negative score.

typedef uint8_t EEU8_TCScore;


// Helper for saturating subtraction (unsigned 8-bit)
// since there is no such thing in standard C++.
constexpr uint8_t subs_u8(uint8_t a, uint8_t b) {
    return a-std::min(a,b);
}

/**
 * Build query profile look up tables for the read.  The query profile look
 * up table is organized as a 1D array indexed by [i][j] where i is the
 * reference character in the current DP column (0=A, 1=C, etc), and j is
 * the segment of the query we're currently working on.
 */
static inline void EEU8_buildQueryProfile(
		const BTDnaString* rd,
		const BTString*    qu,
		const Scoring *    psc,
		SSEData<false,ALN_MAX_ROWS,ALN_MAX_COLS>& d) {
	const size_t len = rd->length();
	const size_t seglen = d.get_sse_seglen(len);
	//const size_t nsses = get_nsses(len, false);
	//d.profbuf_.resizeNoCopy(nsses);
	assert(!d.profbuf_.empty());
	d.maxPen_      = d.maxBonus_ = 0;
	d.lastIter_    = d.lastWord_ = 0;
	d.qprofStride_ = d.gbarStride_ = 2;
	d.bias_ = 0; // no bias needed for end-to-end alignment; just use subtraction
	// For each reference character A, C, G, T, N ...
	for(size_t refc = 0; refc < ALPHA_SIZE; refc++) {
		// For each segment ...
		for(size_t i = 0; i < seglen; i++) {
			size_t j = i;
			uint8_t *qprofWords =
				reinterpret_cast<uint8_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2));
			uint8_t *gbarWords =
				reinterpret_cast<uint8_t*>(d.profbuf_.ptr() + (refc * seglen * 2) + (i * 2) + 1);
			// For each sub-word (byte) ...
			for(size_t k = 0; k < EEU8_NWORDS_PER_REG; k++) {
				int sc = 0;
				*gbarWords = 0;
				if(j < len) {
					int readc = (*rd)[j];
					int readq = (*qu)[j];
					sc = psc->score(readc, (int)(1 << refc), readq - 33);
					// Make score positive, to fit in an unsigned
					sc = -sc;
					assert_range(0, 255, sc);
					size_t j_from_end = len - j - 1;
					if(j < (size_t)psc->gapbar ||
					   j_from_end < (size_t)psc->gapbar)
					{
						// Inside the gap barrier
						*gbarWords = 0xff;
					}
				}
				if(refc == 0 && j == len-1) {
					// Remember which 128-bit word and which smaller word has
					// the final row
					d.lastIter_ = i;
					d.lastWord_ = k;
				}
				if((size_t)sc > d.maxPen_) {
					d.maxPen_ = (size_t)sc;
				}
				*qprofWords = (uint8_t)sc;
				gbarWords++;
				qprofWords++;
				j += seglen; // update offset into query
			}
		}
	}
}


#ifndef SWSSE_INLINE_ONLY
void SwAligner::buildQueryProfileEnd2EndSseU8(bool fw) {
	bool& done = fw ? sseU8fwBuilt_ : sseU8rcBuilt_;
	if(done) {
		return;
	}
	done = true;
	const BTDnaString* rd = fw ? rdfw_ : rdrc_;
	const BTString* qu = fw ? qufw_ : qurc_;
	auto& d = fw ? sseU8fw_ : sseU8rc_;
	EEU8_buildQueryProfile(rd, qu, sc_, d);
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

template<typename TIdxSize>
inline SSERegI EEU8_alignOne(const TIdxSize iter,
			const size_t colstride,
			const SSERegI *pvScore,
			const SSERegI *pvHLoad, const SSERegI *pvELoad,
			SSERegI *pvHStore, SSERegI *pvEStore, SSERegI *pvFStore,
			const SSERegI rfgapo, const SSERegI rfgape, const SSERegI rdgapo, const SSERegI rdgape) {
		// vhilsw: topmost (least sig) word set to 0xff, all other words=0
		SSERegI vhilsw   = sse_setzero_siall();
		sse_set_low_u8(0xff, vhilsw);	
	
		// Set all cells to low value
		SSERegI vf       = sse_setzero_siall();

		// Load H vector from the final row of the previous column
		SSERegI vh = sse_load_siall(pvHLoad + colstride - ROWSTRIDE);
		// Shift N bytes down so that topmost (least sig) cell gets 0
		vh = sse_slli_u8(vh);
		// Fill topmost (least sig) cell with high value
		vh = sse_or_siall(vh, vhilsw);
		
		// For each character in the reference text:
		for(TIdxSize j = 0; j < iter; j++) {
			SSERegI vs0 = sse_load_siall(pvScore);
                        pvScore++;
			SSERegI vs1 = sse_load_siall(pvScore);
                        pvScore++;
			// Load cells from E, calculated previously
			SSERegI ve = sse_load_siall(pvELoad);
			pvELoad += ROWSTRIDE;
			
			// Store cells in F, calculated previously
			vf = sse_subs_epu8(vf, vs1); // veto some ref gap extensions
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Factor in query profile (matches and mismatches)
			vh = sse_subs_epu8(vh, vs0);
			
			// Update H, factoring in E and F
			vh = sse_max_epu8(vh, ve);
			vh = sse_max_epu8(vh, vf);
			
			// Save the new vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update vE value
			SSERegI vtmp = vh;
			vh = sse_subs_epu8(vh, rdgapo);
			vh = sse_subs_epu8(vh, vs1); // veto some read gap opens
			ve = sse_subs_epu8(ve, rdgape);
			ve = sse_max_epu8(ve, vh);

			// Load the next h value
			vh = sse_load_siall(pvHLoad);
			pvHLoad += ROWSTRIDE;
			
			// Save E values
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			
			// Update vf value
			vtmp = sse_subs_epu8(vtmp, rfgapo);
			vf = sse_subs_epu8(vf, rfgape);
			vf = sse_max_epu8(vf, vtmp);
		}

		return vf;
}

#if NBYTES_PER_REG>1
template<typename TIdxSize>
inline void EEU8_lazyF(const SSERegI vf0,
			const TIdxSize iter, const size_t colstride,
			const SSERegI *pvScore,
			SSERegI *pvHStore, SSERegI *pvEStore, SSERegI *pvFStore,
			const SSERegI rfgape, const SSERegI rdgapo) {
		const SSERegI vzero    = sse_setzero_siall();  // needed by sse_anygt_epu8

		// vf from last row gets shifted down by one to overlay the first row
		// rfgape has already been subtracted from it.
		SSERegI vf = sse_slli_u8(vf0);
		
		pvScore += 1;
	        SSERegI vs1 = sse_load_siall(pvScore);
		SSERegI vtmp = sse_load_siall(pvFStore);
		
		vf = sse_subs_epu8(vf, vs1); // veto some ref gap extensions
		vf = sse_max_epu8(vtmp, vf);
		bool anygt;
		sse_anygt_epu8(vf,vtmp,anygt);
	
		// Load after computing cmp, so the result is ready by the time it is tested in while
		SSERegI vh = sse_load_siall(pvHStore);
		SSERegI ve = sse_load_siall(pvEStore);
		
		// If any element of vtmp is greater than H - gap-open...
		TIdxSize j = 0;
		while(anygt) {
			// Store this vf
			sse_store_siall(pvFStore, vf);
			pvFStore += ROWSTRIDE;
			
			// Update vh w/r/t new vf
			vh = sse_max_epu8(vh, vf);
			
			// Save vH values
			sse_store_siall(pvHStore, vh);
			pvHStore += ROWSTRIDE;
			
			// Update E in case it can be improved using our new vh
			vh = sse_subs_epu8(vh, rdgapo);
			vh = sse_subs_epu8(vh, vs1); // veto some read gap opens
			ve = sse_max_epu8(ve, vh);
			sse_store_siall(pvEStore, ve);
			pvEStore += ROWSTRIDE;
			pvScore += 2;
			
			assert_lt(j, iter);
			if(++j == iter) {
				pvScore -= iter*2;
				j = 0;
				pvFStore -= colstride;
				pvHStore -= colstride;
				pvEStore -= colstride;
				vf = sse_slli_u8(vf);
			}
			vs1 = sse_load_siall(pvScore);
			vtmp = sse_load_siall(pvFStore);   // load next vf ASAP

			// Update F with another gap extension
			vf = sse_subs_epu8(vf, rfgape);
			vf = sse_subs_epu8(vf, vs1); // veto some ref gap extensions
			vf = sse_max_epu8(vtmp, vf);
			sse_anygt_epu8(vf,vtmp,anygt);

			// Load after computing cmp, so the result is afailable by the time it is tested
			vh = sse_load_siall(pvHStore);     // load next vh 
			ve = sse_load_siall(pvEStore);     // load next ve
		}
}
#endif

/**
 * Solve the current alignment problem using SSE instructions that operate on 16
 * unsigned 8-bit values packed into a single 128-bit register.
 *
 * Select intput parameters:
 *   profbuf - buffer for query profile & temp vecs
 *   rf      - reference sequence
 *
 * Select output parameters:
 *   pmat    - SSE matrix for holding all E, F, H vectors
 *   btncand - cells we might backtrace from
 */
template<typename TIdxSize>
inline EEU8_TCScore EEU8_alignNucleotides(const SSERegI profbuf[],
					const char   rf[], const TIdxSize rfd,
					SSERegI pmat[],
					const TIdxSize iter, const size_t colstride, const size_t lastWordIdx,
					const TAlScore minsc, const size_t nrow,
					DpBtCandidate btncand[], TIdxSize& btnfilled_,
					const int8_t refGapOpen, const int8_t refGapExtend, const int8_t readGapOpen, const int8_t readGapExtend) {

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

	assert_gt(refGapOpen, 0);
	sse_fill_i8(refGapOpen, rfgapo);
	
	// Set all elts to reference gap extension penalty
	assert_gt(refGapExtend, 0);
	assert_leq(refGapExtend, refGapOpen);
	sse_fill_i8(refGapExtend, rfgape);

	// Set all elts to read gap open penalty
	assert_gt(readGapOpen, 0);
	sse_fill_i8(readGapOpen, rdgapo);
	
	// Set all elts to read gap extension penalty
	assert_gt(readGapExtend, 0);
	assert_leq(readGapExtend, readGapOpen);
	sse_fill_i8(readGapExtend, rdgape);

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
	  SSERegI *pvHTmp = pmat + SSEMatrixConsts::TMP;
	  SSERegI *pvETmp = pmat + SSEMatrixConsts::E;
	  SSERegI vlo      = sse_setzero_siall();
	
	  for(size_t i = 0; i < iter; i++) {
		sse_store_siall(pvETmp, vlo);
		sse_store_siall(pvHTmp, vlo); // start high in end-to-end mode
		pvETmp += ROWSTRIDE;
		pvHTmp += ROWSTRIDE;
	  }
	}

	// These are swapped just before the innermost loop
	SSERegI *pvHLoad  = pmat + SSEMatrixConsts::TMP;
	SSERegI *pvHStore = pmat + SSEMatrixConsts::H;
	SSERegI *pvELoad  = pmat + SSEMatrixConsts::E;
	SSERegI *pvEStore = pmat + colstride + SSEMatrixConsts::E;
	SSERegI *pvFStore = pmat + SSEMatrixConsts::F;
	
	// Maximum score in final row
	EEU8_TCScore lrmax = MIN_U8;

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
		assert(pvFStore == (pmat + i*colstride + SSEMatrixConsts::F));
		assert(pvHStore == (pmat + i*colstride + SSEMatrixConsts::H));
		
		// Fetch the appropriate query profile.  Note that elements of rf must
		// be numbers, not masks.
		size_t off = (size_t)std::countr_zero( uint8_t(rf[i]) ) * iter * 2;
		// points into the query profile
		const SSERegI *pvScore = profbuf + off; // even elts = query profile, odd = gap barrier
	
		SSERegI vf = EEU8_alignOne(iter,
                        	colstride,
                        	pvScore,
                        	pvHLoad, pvELoad,
                        	pvHStore, pvEStore, pvFStore,
                        	rfgapo, rfgape, rdgapo, rdgape);

		if constexpr(NBYTES_PER_REG>1) {
			EEU8_lazyF(vf,
				iter, colstride,
				pvScore,
				pvHStore, pvEStore, pvFStore,
				rfgape, rdgapo);
		}
		// else, no need for lazyF
		pvHLoad = pvHStore;    // new pvHLoad = pvHStore
		
		// Note: we may not want to extract from the final row
		EEU8_TCScore lr = ((EEU8_TCScore*)(pvHLoad))[lastWordIdx];
		TAlScore sc = (TAlScore)(lr - 0xff);
		if(lr > lrmax) {
			lrmax = lr;
		}
		if(sc >= minsc) {
			// Yes, this is legit
			btncand[btnfilled].init(nrow-1, i, lr);
			btnfilled++;
		}

		// pvHLoad are already where they need to be
		
		// Adjust the load and store vectors here.  
		pvHStore = pvHStore + colstride;
		pvELoad  = pvELoad  + colstride;
		pvEStore = pvEStore + colstride;
		pvFStore = pvFStore + colstride;
	}


	btnfilled_ = btnfilled;  // pass it out
	return lrmax;
}


/**
 * Like the previous alignNucleotides, but the Scalar version allows to forfeit some extra strutures.
 * Note that this version does not return the E, F, H vectors.
 * If those are needed, e.g. when (lrmax - 0xff)>=minsc, i.e. btnfilled>0, use the regular, full version.
 *
 * Select intput parameters:
 *   profbuf - buffer for query profile & temp vecs
 *   rf      - reference sequence
 *
 */
template<typename TIdxSize=uint16_t, uint16_t MAX_ITER=151, uint16_t MAX_RB=5>
inline EEU8_TCScore EEU8_alignNucleotidesScalar(const uint8_t profbuf[],
					const char   rf[], const TIdxSize rfd,
					const size_t nrow,
					const int8_t refGapOpen, const int8_t refGapExtend, const int8_t readGapOpen, const int8_t readGapExtend) {
        class TPackedEH {
	private:
		uint32_t val;
	public:
		constexpr TPackedEH() : val(0) {};
		constexpr TPackedEH(const TPackedEH& other) : val(other.val) {}
		constexpr TPackedEH& operator=(const TPackedEH& other) { val = other.val; return *this;}

		constexpr void set(uint8_t E_even, uint8_t H_even, uint8_t E_odd, uint8_t H_odd) { val =(uint32_t(H_odd)<<24) | (uint32_t(E_odd)<<16) |  (uint32_t(H_even)<<8) | uint32_t(E_even);}
		constexpr void set(uint8_t E_even, uint8_t H_even) { val = (uint32_t(H_even)<<8) | uint32_t(E_even);}

		// preserve odd, update even (low 16 bits)
		constexpr void set_even(uint8_t E, uint8_t H) { val = (val & uint32_t(0xffff0000U)) | ((uint32_t(H)<<8) | uint32_t(E));}
		// preserve even, update odd (high 16 bits)
		constexpr void set_odd(uint8_t E, uint8_t H) { val = (val & uint32_t(0xffffU)) | ((uint32_t(H)<<24) | (uint32_t(E)<<16));}

		constexpr uint8_t get_E_even() const {return uint8_t(val); }
		constexpr uint8_t get_H_even() const {return uint8_t(val>>8); }
		constexpr uint8_t get_E_odd() const {return uint8_t(val>>16); }
		constexpr uint8_t get_H_odd() const {return uint8_t(val>>24); }
	};

        class TPackedScore {
	private:
		uint32_t val;
	public:
		constexpr TPackedScore() : val(0) {};
		constexpr TPackedScore(const uint32_t packed_val) : val(packed_val) {};
		// Note: j must be even
		constexpr TPackedScore(const uint8_t *pvScore, const uint16_t j) { val = ((uint32_t *)(pvScore+(2*j)))[0]; }

		constexpr TPackedScore(const TPackedScore& other) : val(other.val) {}
		constexpr TPackedScore& operator=(const TPackedScore& other) { val = other.val; return *this;}

		constexpr uint8_t get_vs0_even() const {return uint8_t(val); }
		constexpr uint8_t get_vs1_even() const {return uint8_t(val>>8); }
		constexpr uint8_t get_vs0_odd() const {return uint8_t(val>>16); }
		constexpr uint8_t get_vs1_odd() const {return uint8_t(val>>24); }

		constexpr uint32_t operator()() const {return val;}
	};

        class TPackedScoreHalf {
	private:
		uint16_t val;
	public:
		constexpr TPackedScoreHalf() : val(0) {};
		// Note: j must be even
		constexpr TPackedScoreHalf(const uint8_t *pvScore, const uint16_t j) { val = ((uint16_t *)(pvScore+(2*j)))[0]; }

		constexpr TPackedScoreHalf(const TPackedScoreHalf& other) : val(other.val) {}
		constexpr TPackedScoreHalf& operator=(const TPackedScoreHalf& other) { val = other.val; return *this;}

		constexpr uint8_t get_vs0_even() const {return uint8_t(val); }
		constexpr uint8_t get_vs1_even() const {return uint8_t(val>>8); }

		constexpr uint16_t operator()() const {return val;}
	};

        constexpr uint16_t iter = MAX_ITER;

	assert_gt(refGapOpen, 0);
	uint8_t rfgapo = uint8_t(refGapOpen);
	
	// Set all elts to reference gap extension penalty
	assert_gt(refGapExtend, 0);
	assert_leq(refGapExtend, refGapOpen);
	uint8_t rfgape = uint8_t(refGapExtend);

	// Set all elts to read gap open penalty
	assert_gt(readGapOpen, 0);
	uint8_t rdgapo = uint8_t(readGapOpen);
	
	// Set all elts to read gap extension penalty
	assert_gt(readGapExtend, 0);
	assert_leq(readGapExtend, readGapOpen);
	uint8_t rdgape = uint8_t(readGapExtend);

	// Load the procbuf in local memeory
	// as we will read over it many times
	uint32_t loc_procbuf[MAX_RB][(MAX_ITER+1)/2];
	for (int ir=0; ir<MAX_RB; ir++) {
		size_t off = (size_t)( ir ) * iter * 2;
		// points into the query profile
		const uint8_t *pvScore = profbuf + off;
		for(TIdxSize j = 0; j < ((iter-1)/2); j++) {
			const TPackedScore full_score(pvScore,j*2);
			loc_procbuf[ir][j] = full_score();
		}
		if constexpr((iter%2)!=0) {
			// Last Even step
		        const TPackedScoreHalf full_score(pvScore, iter-1);
			loc_procbuf[ir][(iter-1)/2] = full_score();
		}
	}

	//
	//
	// Maximum score in final row
	EEU8_TCScore lrmax = MIN_U8;

	// Set all elts to reference gap open penalty
	// iAutomatic initialize the H and E vectors in the first matrix column
	TPackedEH bufEH[(iter+1)/2];

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
		
		// Fetch the appropriate query profile.  Note that elements of rf must
		// be numbers, not masks.
		const size_t ir = (size_t)std::countr_zero( uint8_t(rf[i]) );
	
		// scalar version of EEU8_alignOne()
		{
		  // vhilsw: topmost (least sig) word set to 0xff, all other words=0
		  uint8_t vhilsw   = 0xff;
	
		  // Set all cells to low value
		  uint8_t vf       = 0;

		  // in scalar mode, we always start high
		  uint8_t vh = vhilsw; //0xff

		  // For each character in the reference text
#ifdef OMPGPU
		  // Full unrolling allows bufEH to be mapped to the GPU register file
		  // Only downsides for the CPUs, that have fewer registers
#pragma omp unroll full
#endif
		  for(TIdxSize j = 0; j < ((iter-1)/2); j++) {
		    // Load cells from E and H, calculated previously
		    TPackedEH full_eh(bufEH[j]);
		    // Load the scores
		    const TPackedScore full_score(loc_procbuf[ir][j]);
		    // Even step
		    {
		  	// Load cells from E and H, calculated previously
			uint8_t ve = full_eh.get_E_even();
		  	uint8_t vh_next = full_eh.get_H_even();

			// Store cells in F, calculated previously
			vf = subs_u8(vf, full_score.get_vs1_even()); // veto some ref gap extensions
		  	
		  	// Factor in query profile (matches and mismatches)
		  	vh = subs_u8(vh, full_score.get_vs0_even());
		  	
		  	// Update H, factoring in E and F
		  	vh = std::max(vh, ve);
		  	vh = std::max(vh, vf);
		  	
		  	// Save the new vH values
		  	uint8_t vtmp = vh;
		  	
		  	// Update vE value
		  	vh = subs_u8(vh, rdgapo);
		  	vh = subs_u8(vh, full_score.get_vs1_even()); // veto some read gap opens
		  	ve = subs_u8(ve, rdgape);
		  	ve = std::max(ve, vh);

		  	// Save E and H values for next i round
		  	full_eh.set_even(ve,vtmp);
		  	
		  	// Update vf value for next round
		  	vtmp = subs_u8(vtmp, rfgapo);
		  	vf = subs_u8(vf, rfgape);
		  	vf = std::max(vf, vtmp);
			vh = vh_next;
		    }
		    // Odd step
		    {
		  	// Load cells from E and H, calculated previously
			uint8_t ve = full_eh.get_E_odd();
		  	uint8_t vh_next = full_eh.get_H_odd();

			// Store cells in F, calculated previously
			vf = subs_u8(vf, full_score.get_vs1_odd()); // veto some ref gap extensions
		  	
		  	// Factor in query profile (matches and mismatches)
		  	vh = subs_u8(vh, full_score.get_vs0_odd());
		  	
		  	// Update H, factoring in E and F
		  	vh = std::max(vh, ve);
		  	vh = std::max(vh, vf);
		  	
		  	// Save the new vH values
		  	uint8_t vtmp = vh;
		  	
		  	// Update vE value
		  	vh = subs_u8(vh, rdgapo);
		  	vh = subs_u8(vh, full_score.get_vs1_odd()); // veto some read gap opens
		  	ve = subs_u8(ve, rdgape);
		  	ve = std::max(ve, vh);

		  	// Save E and H values for next i round
		  	full_eh.set_odd(ve,vtmp);
		  	
		  	// Update vf value for next round
		  	vtmp = subs_u8(vtmp, rfgapo);
		  	vf = subs_u8(vf, rfgape);
		  	vf = std::max(vf, vtmp);
			vh = vh_next;
		    }
		    // save the packed result back for the next i loop
		    bufEH[j] = full_eh;
		  }
		  if constexpr((iter%2)!=0) {
			// Last Even step
		  	// Load cells from E and H, calculated previously
			const TPackedEH full_eh(bufEH[iter/2]);
		        // Load the scores
			// Use the Full Packed class for consistency
		        const TPackedScore full_score(loc_procbuf[ir][(iter-1)/2]);

			uint8_t ve = full_eh.get_E_even();
			// no next round, do not need vh_next

			// Store cells in F, calculated previously
			vf = subs_u8(vf, full_score.get_vs1_even()); // veto some ref gap extensions
		  	
		  	// Factor in query profile (matches and mismatches)
		  	vh = subs_u8(vh, full_score.get_vs0_even());
		  	
		  	// Update H, factoring in E and F
		  	vh = std::max(vh, ve);
		  	vh = std::max(vh, vf);
		  	
		  	// Save the new vH values
		  	uint8_t vtmp = vh;
		  	
		  	// Update vE value
		  	vh = subs_u8(vh, rdgapo);
		  	vh = subs_u8(vh, full_score.get_vs1_even()); // veto some read gap opens
		  	ve = subs_u8(ve, rdgape);
		  	ve = std::max(ve, vh);

		  	// Save E and H values for next i round
		  	bufEH[iter/2].set(ve,vtmp);
		  	
		  	// no next round
		  }
		}

		// Note: we may not want to extract from the final row
		//       EEU8_TCScore == uint8_t == uint8_t
		// Was: EEU8_TCScore lr = bufEH2[iter - 1].get_H();
		EEU8_TCScore lr = ((iter%2)!=0) ? bufEH[iter/2].get_H_even() : bufEH[(iter-1)/2].get_H_odd();
		if(lr > lrmax) {
			lrmax = lr;
		}
	}

	return lrmax;
}

#ifndef SWSSE_INLINE_ONLY
/**
 * Align read 'rd' to reference using read & reference information given
 * last time init() was called.
 */
bool SwAligner::alignEnd2EndSseU8(
	TAlScore& best)    // best alignment score observed in DP matrix
{
#ifdef ENABLE_SSE_METRICS
	constexpr bool debug = false;
#endif

	assert(initedRef() && initedRead());
	assert_eq(STATE_INITED, state_);

	assert_leq(rdlen_, qu_->length());
	assert_eq(rd_->length(), qu_->length());
	assert_geq(sc_->gapbar, 1);
	assert(repOk());
#ifndef NDEBUG
	for(size_t i = 0; i < (size_t)rflen_; i++) {
		assert_range(0, 16, (int)rf_[i]);
	}
#endif

	auto& d = fw_ ? sseU8fw_ : sseU8rc_;
#ifdef ENABLE_SSE_METRICS
	SSEMetrics& met = sseMet_;
	if(!debug) met.dp++;
#endif
	buildQueryProfileEnd2EndSseU8(fw_);
	assert(!d.profbuf_.empty());

	assert_eq(0, d.maxBonus_);
	const size_t iter =
		(dpRows() + (EEU8_NWORDS_PER_REG-1)) / EEU8_NWORDS_PER_REG; // iter = segLen
        const size_t lastWordIdx = EEU8_NWORDS_PER_REG*(d.lastIter_*ROWSTRIDE)+d.lastWord_;


	assert_gt(sc_->refGapOpen(), 0);
	assert_leq(sc_->refGapOpen(), MAX_U8);
	
	// Set all elts to reference gap extension penalty
	assert_gt(sc_->refGapExtend(), 0);
	assert_leq(sc_->refGapExtend(), MAX_U8);
	assert_leq(sc_->refGapExtend(), sc_->refGapOpen());

	// Set all elts to read gap open penalty
	assert_gt(sc_->readGapOpen(), 0);
	assert_leq(sc_->readGapOpen(), MAX_U8);
	
	// Set all elts to read gap extension penalty
	assert_gt(sc_->readGapExtend(), 0);
	assert_leq(sc_->readGapExtend(), MAX_U8);
	assert_leq(sc_->readGapExtend(), sc_->readGapOpen());

	assert_gt(sc_->gapbar, 0);

	colstop_ = rflen_ - 1;
	lastsolcol_ = 0;

	// should never get in here, but just in case
	if (rflen_==0) {
		btncand_.clear();
#ifdef ENABLE_SSE_METRICS
		if(!debug) met.dpfail++;
#endif
		best = MIN_I64;
		return false;
	}

	assert_eq(int(d.mat_.wperf_), int(EEU8_NWORDS_PER_REG));
	d.mat_.init(dpRows(), rflen_);

	assert_leq(iter,   (size_t)MAX_U16);
	assert_leq(rflen_, (size_t)MAX_U16);
	uint16_t btnfilled = 0;
	btncand_.resizeNoCopy(rflen_); // cannot be bigger that this

#ifdef SSE_FAST_SCALAR
	const EEU8_TCScore lrmax = EEU8_alignNucleotidesScalar<uint16_t>(d.profbuf_.ptr(), rf_, rflen_,
					d.mat_.ptr(),
                                        iter, d.mat_.colstride(), lastWordIdx,
					minsc_, dpRows(),
					btncand_.ptr(), btnfilled,
					sc_->refGapOpen(), sc_->refGapExtend(), sc_->readGapOpen(), sc_->readGapExtend());
#else
	const EEU8_TCScore lrmax = EEU8_alignNucleotides<uint16_t>(d.profbuf_.ptr(), rf_, rflen_,
					d.mat_.ptr(),
                                        iter, d.mat_.colstride(), lastWordIdx,
					minsc_, dpRows(),
					btncand_.ptr(), btnfilled,
					sc_->refGapOpen(), sc_->refGapExtend(), sc_->readGapOpen(), sc_->readGapExtend());
#endif
	btncand_.trim(btnfilled);
	state_ = STATE_ALIGNED;
	cural_ = 0;

#ifdef ENABLE_SSE_METRICS
	// Update metrics
	if(!debug) {
		size_t ninner = (rflen_) * iter;
		met.col   += (rflen_);             // DP columns
		met.cell  += (ninner * EEU8_NWORDS_PER_REG); // DP cells
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
	if(lrmax == MIN_U8) {
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
bool SwAligner::backtraceNucleotidesEnd2EndSseU8(
	TAlScore       escore, // in: expected score
	SwResult&      res,    // out: store results (edits and scores) here
	size_t&        off,    // out: store diagonal projection of origin
	size_t&        nbts,   // out: # backtracks
	size_t         row,    // start in this row
	size_t         col)    // start in this column
{
	assert_lt(row, dpRows());
	assert_lt(col, (size_t)(rflen_));
	auto& d = fw_ ? sseU8fw_ : sseU8rc_;
#ifdef ENABLE_SSE_METRICS
	SSEMetrics& met = sseMet_;
	met.bt++;
#endif
	assert(!d.profbuf_.empty());
	assert_lt(row, rd_->length());
	btnstack_.clear(); // empty the backtrack stack
	btcells_.clear();  // empty the cells-so-far list
	AlnScore score;
	score.reset();
	score.basesAligned_ = std::numeric_limits<int>::min();
	score.edits_ = std::numeric_limits<int>::min();
	size_t origCol = col;
	size_t gaps = 0, readGaps = 0, refGaps = 0;
	res.alres.reset();
	auto& ned = res.alres.ned();
	assert(ned.empty());
	assert_gt(dpRows(), row);
	size_t trimEnd = dpRows() - row - 1; 
	size_t trimBeg = 0;
	size_t ct = SSEMatrixConsts::H; // cell type
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
		int readc = (*rd_)[row];
		int refm  = (int)rf_[col];
		int readq = (*qu_)[row];
		assert_leq(col, origCol);
		// Get score in this cell
		bool empty = false;
		bool reportedThru, canMoveThru, branch = false;
		int cur = SSEMatrixConsts::H;
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
				const int offsetsc = -0xff;
				// Move to beginning of column/row
				if(ct == SSEMatrixConsts::E) { // AKA rdgap
					assert_gt(col, 0);
					TAlScore sc_cur = ((EEU8_TCScore*)(cur_vec + SSEMatrixConsts::E))[rowelt] + offsetsc;
					assert(gapsAllowed);
					// Currently in the E matrix; incoming transition must come from the
					// left.  It's either a gap open from the H matrix or a gap extend from
					// the E matrix.
					// TODO: save and restore origMask as well as mask
					int origMask = 0, mask = 0;
					// Get H score of cell to the left
					TAlScore sc_h_left = ((EEU8_TCScore*)(left_vec + SSEMatrixConsts::H))[left_rowelt] + offsetsc;
					if(sc_h_left > floorsc && sc_h_left - sc_->readGapOpen() == sc_cur) {
						mask |= (1 << 0);
					}
					// Get E score of cell to the left
					TAlScore sc_e_left = ((EEU8_TCScore*)(left_vec + SSEMatrixConsts::E))[left_rowelt] + offsetsc;
					if(sc_e_left > floorsc && sc_e_left - sc_->readGapExtend() == sc_cur) {
						mask |= (1 << 1);
					}
					origMask = mask;
					assert(origMask > 0 || sc_cur <= sc_->match());
					if(d.mat_.isEMaskSet(row, col)) {
						mask = (d.mat_.masks(row,col) >> 8) & 3;
					}
					if(mask == 3) {
						// Pick H -> E cell
						cur = SW_BT_OALL_READ_OPEN;
						d.mat_.eMaskSet(row, col, 2); // might choose E later
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
				} else if(ct == SSEMatrixConsts::F) { // AKA rfgap
					assert_gt(row, 0);
					assert(gapsAllowed);
					TAlScore sc_h_up = ((EEU8_TCScore*)(up_vec  + SSEMatrixConsts::H))[up_rowelt] + offsetsc;
					TAlScore sc_f_up = ((EEU8_TCScore*)(up_vec  + SSEMatrixConsts::F))[up_rowelt] + offsetsc;
					TAlScore sc_cur  = ((EEU8_TCScore*)(cur_vec + SSEMatrixConsts::F))[rowelt] + offsetsc;
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
						mask = (d.mat_.masks(row,col) >> 11) & 3;
					}
					if(mask == 3) {
						// I chose the H cell
						cur = SW_BT_OALL_REF_OPEN;
						d.mat_.fMaskSet(row, col, 2); // might choose E later
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
					assert_eq(SSEMatrixConsts::H, ct);
					TAlScore sc_cur      = ((EEU8_TCScore*)(cur_vec + SSEMatrixConsts::H))[rowelt]    + offsetsc;
					TAlScore sc_f_up     = ((EEU8_TCScore*)(up_vec  + SSEMatrixConsts::F))[up_rowelt] + offsetsc;
					TAlScore sc_h_up     = ((EEU8_TCScore*)(up_vec  + SSEMatrixConsts::H))[up_rowelt] + offsetsc;
					TAlScore sc_h_left   = col > 0 ? (((EEU8_TCScore*)(left_vec   + SSEMatrixConsts::H))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_e_left   = col > 0 ? (((EEU8_TCScore*)(left_vec   + SSEMatrixConsts::E))[left_rowelt]   + offsetsc) : floorsc;
					TAlScore sc_h_upleft = col > 0 ? (((EEU8_TCScore*)(upleft_vec + SSEMatrixConsts::H))[upleft_rowelt] + offsetsc) : floorsc;
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
						mask = (d.mat_.masks(row,col) >> 2) & 31;
					}
					assert(gapsAllowed || mask == (1 << 4) || mask == 0);
					int opts = std::popcount(uint8_t(mask));
					int select = -1;
					if(opts == 1) {
						select = std::countr_zero(uint8_t(mask));
						assert_geq(mask, 0);
						d.mat_.hMaskSet(row, col, 0);
					} else if(opts > 1) {
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
				assert(!empty || !canMoveThru || ct == SSEMatrixConsts::H);
			}
		}
		//cerr << "reportedThrough rejected (" << row << ", " << col << ")" << endl;
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
			assert_eq(SSEMatrixConsts::H, ct);
			btcells_.expand();
			btcells_.back().first = row;
			btcells_.back().second = col;
			// This cell is at the end of a legitimate alignment
			trimBeg = row;
			assert_eq(0, trimBeg);
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
				int refNmask = (int)rf_[col];
				assert_gt(refNmask, 0);
				int m = matchesEx(readC, refNmask);
				ct = SSEMatrixConsts::H;
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
				assert_geq((int)(rdlen_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrixConsts::H;
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
				assert_geq((int)(rdlen_-row-1), sc_->gapbar-1);
				row--;
				ct = SSEMatrixConsts::F;
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
					mask2dna[(int)rf_[col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdlen_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrixConsts::H;
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
					mask2dna[(int)rf_[col]],
					'-',
					EDIT_TYPE_READ_GAP);
				assert(e.repOk());
				assert(ned.empty() || ned.back().pos >= row);
				ned.push_back(e);
				assert_geq(row, (size_t)sc_->gapbar);
				assert_geq((int)(rdlen_-row-1), sc_->gapbar-1);
				col--;
				ct = SSEMatrixConsts::E;
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
	assert_eq(SSEMatrixConsts::H, ct);
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
#ifndef NDEBUG
		//assert(!d.mat_.reportedThrough(rw, cl));
		//d.mat_.setReportedThrough(rw, cl);
		assert(d.mat_.reportedThrough(rw, cl));
#endif
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
	int readC = (*rd_)[row];      // get last char in read
	int refNmask = (int)rf_[col]; // get last ref char ref involved in aln
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
	assert_lt(col, (size_t)rflen_);
	score.gaps_ = gaps;
	score.edits_ = (int)ned.size();
	score.basesAligned_ = (int)(rdlen_ - trimBeg - trimEnd - score.edits_);
	res.alres.setScore(score);
	res.alres.setShape(
		refidx_,                  // ref id
		off + rect_->refl, // 0-based ref offset
		reflen_,                  // length of entire reference
		fw_,                      // aligned to Watson?
		rdlen_,                   // read length
		true,                     // pretrim soft?
		0,                        // pretrim 5' end
		0,                        // pretrim 3' end
		true,                     // alignment trim soft?
		fw_ ? trimBeg : trimEnd,  // alignment trim 5' end
		fw_ ? trimEnd : trimBeg); // alignment trim 3' end
	size_t refns = 0;
	for(size_t i = col; i <= origCol; i++) {
		if((int)rf_[i] > 15) {
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
		refstr.append(std::countr_zero(uint8_t(rf_[i])));
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
#endif /* SWSSE_INLINE_ONLY */

