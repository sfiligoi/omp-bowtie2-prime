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

#ifndef ALIGNER_SWSSE_H_
#define ALIGNER_SWSSE_H_

#include "ds.h"
#include "mem_ids.h"
#include "random_source.h"
#include "scoring.h"
#include "mask.h"
#include "sse_util.h"
#include <strings.h>

static constexpr size_t EEU8_NWORDS_PER_REG  = NBYTES_PER_REG;
//static const size_t EEU8_NBITS_PER_WORD  = 8;
static constexpr size_t EEU8_NBYTES_PER_WORD = 1;

static constexpr size_t EEI16_NWORDS_PER_REG  = NBYTES_PER_REG/2;
static constexpr size_t EEI16_NBITS_PER_WORD  = 16;
static constexpr size_t EEI16_NBYTES_PER_WORD = 2;

// Hardcode for now. May want to pass in Makefile
#define SSE_MAX_ROWS 160
#define SSE_MAX_COLS 200


struct SSEMetrics {
	
	SSEMetrics() { reset(); }

	void clear() { reset(); }
	void reset() {
		dp = dpsat = dpfail = dpsucc = 
		col = cell = inner = fixup =
		gathsol = bt = btfail = btsucc = btcell =
		corerej = nrej = 0;
	}

	void merge(const SSEMetrics& o) {
		dp       += o.dp;
		dpsat    += o.dpsat;
		dpfail   += o.dpfail;
		dpsucc   += o.dpsucc;
		col      += o.col;
		cell     += o.cell;
		inner    += o.inner;
		fixup    += o.fixup;
		gathsol  += o.gathsol;
		bt       += o.bt;
		btfail   += o.btfail;
		btsucc   += o.btsucc;
		btcell   += o.btcell;
		corerej  += o.corerej;
		nrej     += o.nrej;
	}

	uint64_t dp;       // DPs tried
	uint64_t dpsat;    // DPs saturated
	uint64_t dpfail;   // DPs failed
	uint64_t dpsucc;   // DPs succeeded
	uint64_t col;      // DP columns
	uint64_t cell;     // DP cells
	uint64_t inner;    // DP inner loop iters
	uint64_t fixup;    // DP fixup loop iters
	uint64_t gathsol;  // DP gather solution cells found
	uint64_t bt;       // DP backtraces
	uint64_t btfail;   // DP backtraces failed
	uint64_t btsucc;   // DP backtraces succeeded
	uint64_t btcell;   // DP backtrace cells traversed
	uint64_t corerej;  // DP backtrace core rejections
	uint64_t nrej;     // DP backtrace N rejections
};

/**
 * Encapsulates matrix information calculated by the SSE aligner.
 *
 * Matrix memory is laid out as follows:
 *
 * - Elements (individual cell scores) are packed into SSERegI vectors
 * - Vectors are packed into quartets, quartet elements correspond to: a vector
 *   from E, one from F, one from H, and one that's "reserved"
 * - Quartets are packed into columns, where the number of quartets is
 *   determined by the number of query characters divided by the number of
 *   elements per vector
 *
 * Regarding the "reserved" element of the vector quartet: we use it for two
 * things.  First, we use the first column of reserved vectors to stage the
 * initial column of H vectors.  Second, we use the "reserved" vectors during
 * the backtrace procedure to store information about (a) which cells have been
 * traversed, (b) whether the cell is "terminal" (in local mode), etc.
 */
struct SSEMatrixConsts {

	// Each matrix element is a quartet of vectors.  These constants are used
	// to identify members of the quartet.
	constexpr static uint16_t E   = 0;
	constexpr static uint16_t F   = 1;
	constexpr static uint16_t H   = 2;
	constexpr static uint16_t TMP = 3;
	constexpr static uint16_t nvecPerCell_ = 4; // # vectors per matrix cell (4)
};

template<bool i16, uint16_t max_rows = SSE_MAX_ROWS, uint16_t max_cols = SSE_MAX_COLS>
struct SSEMatrix {

	// Each matrix element is a quartet of vectors.  These constants are used
	// to identify members of the quartet.
	constexpr static uint16_t E   = SSEMatrixConsts::E;
	constexpr static uint16_t F   = SSEMatrixConsts::F;
	constexpr static uint16_t H   = SSEMatrixConsts::H;
	constexpr static uint16_t TMP = SSEMatrixConsts::TMP;
	constexpr static uint16_t nvecPerCell_ = SSEMatrixConsts::nvecPerCell_;

	SSEMatrix(int cat = 0) : matbuf_(cat) { }

	void set_alloc(BTAllocator *alloc, bool propagate_alloc=true) {
		matbuf_.set_alloc(alloc, propagate_alloc);
		masks_.set_alloc(alloc, propagate_alloc);
		reset_.set_alloc(alloc, propagate_alloc);
	}

	void set_alloc(std::pair<BTAllocator *, bool> arg) {
		set_alloc(arg.first,arg.second);
	}

	/**
	 * Return a pointer to the matrix buffer.
	 */
	inline SSERegI *ptr() {
		assert(inited_);
		return matbuf_.ptr();
	}
	
	/**
	 * Return a pointer to the E vector at the given row and column.  Note:
	 * here row refers to rows of vectors, not rows of elements.
	 */
	inline SSERegI* evec(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + E;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}

	/**
	 * Like evec, but it's allowed to ask for a pointer to one column after the
	 * final one.
	 */
	inline SSERegI* evecUnsafe(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_leq(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + E;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}

	/**
	 * Return a pointer to the F vector at the given row and column.  Note:
	 * here row refers to rows of vectors, not rows of elements.
	 */
	inline SSERegI* fvec(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + F;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}

	/**
	 * Return a pointer to the H vector at the given row and column.  Note:
	 * here row refers to rows of vectors, not rows of elements.
	 */
	inline SSERegI* hvec(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + H;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}

	/**
	 * Return a pointer to the TMP vector at the given row and column.  Note:
	 * here row refers to rows of vectors, not rows of elements.
	 */
	inline SSERegI* tmpvec(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_lt(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + TMP;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}

	/**
	 * Like tmpvec, but it's allowed to ask for a pointer to one column after
	 * the final one.
	 */
	inline SSERegI* tmpvecUnsafe(size_t row, size_t col) {
		assert_lt(row, nvecrow_);
		assert_leq(col, nveccol_);
		size_t elt = row * rowstride() + col * colstride() + TMP;
		assert_lt(elt, matbuf_.size());
		return ptr() + elt;
	}
	
	/**
	 * Given a number of rows (nrow), a number of columns (ncol), and the
	 * number of words to fit inside a single SSERegI vector, initialize the
	 * matrix buffer to accomodate the needed configuration of vectors.
	 */
	void init(
		size_t nrow,
		size_t ncol)
	{
	 nrow_ = nrow;
	 ncol_ = ncol;
	 assert_leq(nrow,max_rows);
	 assert_leq(ncol,max_cols);
	 nvecPerCol_ = get_vecPerCol(nrow);
	 {
		// compute vecshift_=log2(wperv)
		size_t wp = wperv_;
		assert(wp == (NBYTES_PER_REG/2) || wp == NBYTES_PER_REG);
	        uint8_t vecshift = 0;
        	while( wp>>=1 ) vecshift++;
		nvecrow_ = (nrow + (wperv_-1)) >> vecshift;
	 }
	 nveccol_ = ncol;
	 colstride_ = nvecPerCol_ * nvecPerCell_;
	 rowstride_ = nvecPerCell_;
	 inited_ = true;
	}
	
	/**
	 * Return the number of SSERegI's you need to skip over to get from one
	 * cell to the cell one column over from it.
	 */
	inline size_t colstride() const { return colstride_; }

	/**
	 * Return the number of SSERegI's you need to skip over to get from one
	 * cell to the cell one row down from it.
	 */
	inline size_t rowstride() const { return rowstride_; }

	/**
	 * Given a row, col and matrix (i.e. E, F or H), return the corresponding
	 * element.
	 */
	int eltSlow(size_t row, size_t col, size_t mat) const {
		assert_lt(row, nrow_);
		assert_lt(col, ncol_);
		assert_leq(mat, 3);
		// Move to beginning of column/row
		size_t rowelt = row / nvecrow_;
		size_t rowvec = row % nvecrow_;
		size_t eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
		if constexpr(!i16) {
			return (int)((uint8_t*)(matbuf_.ptr() + eltvec))[rowelt];
		} else {
			return (int)((int16_t*)(matbuf_.ptr() + eltvec))[rowelt];
		}
	}
	
	/**
	 * Given a row, col and matrix (i.e. E, F or H), return the corresponding
	 * element.
	 */
	inline int elt(size_t row, size_t col, size_t mat) const {
		assert(inited_);
		assert_lt(row, nrow_);
		assert_lt(col, ncol_);
		assert_lt(mat, 3);
		// Move to beginning of column/row
		size_t rowelt = row / nvecrow_;
		size_t rowvec = row % nvecrow_;
		size_t eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
		assert_lt(eltvec, matbuf_.size());
		if(wperv_ == NBYTES_PER_REG) {
			return (int)((uint8_t*)(matbuf_.ptr() + eltvec))[rowelt];
		} else {
			assert_eq((NBYTES_PER_REG/2), wperv_);
			return (int)((int16_t*)(matbuf_.ptr() + eltvec))[rowelt];
		}
	}

	/**
	 * Return the element in the E matrix at element row, col.
	 */
	inline int eelt(size_t row, size_t col) const {
		return elt(row, col, E);
	}

	/**
	 * Return the element in the F matrix at element row, col.
	 */
	inline int felt(size_t row, size_t col) const {
		return elt(row, col, F);
	}

	/**
	 * Return the element in the H matrix at element row, col.
	 */
	inline int helt(size_t row, size_t col) const {
		return elt(row, col, H);
	}
	
	/**
	 * Return true iff the given cell has its reportedThru bit set.
	 */
	inline bool reportedThrough(
		size_t row,          // current row
		size_t col) const    // current column
	{
		return (masks_[row][col] & (1 << 0)) != 0;
	}

	/**
	 * Set the given cell's reportedThru bit.
	 */
	inline void setReportedThrough(
		size_t row,          // current row
		size_t col)          // current column
	{
		masks_[row][col] |= (1 << 0);
	}

	/**
	 * Return true iff the H mask has been set with a previous call to hMaskSet.
	 */
	bool isHMaskSet(
		size_t row,          // current row
		size_t col) const    // current column
	{
		return (masks_[row][col] & (1 << 1)) != 0;
	}

	/**
	 * Set the given cell's H mask.  This is the mask of remaining legal ways to
	 * backtrack from the H cell at this coordinate.  It's 5 bits long and has
	 * offset=2 into the 16-bit field.
	 */
	void hMaskSet(
		size_t row,          // current row
		size_t col,          // current column
		int mask)
	{
		assert_lt(mask, 32);
		masks_[row][col] &= ~(31 << 1);
		masks_[row][col] |= (1 << 1 | mask << 2);
	}

	/**
	 * Return true iff the E mask has been set with a previous call to eMaskSet.
	 */
	bool isEMaskSet(
		size_t row,          // current row
		size_t col) const    // current column
	{
		return (masks_[row][col] & (1 << 7)) != 0;
	}

	/**
	 * Set the given cell's E mask.  This is the mask of remaining legal ways to
	 * backtrack from the E cell at this coordinate.  It's 2 bits long and has
	 * offset=8 into the 16-bit field.
	 */
	void eMaskSet(
		size_t row,          // current row
		size_t col,          // current column
		int mask)
	{
		assert_lt(mask, 4);
		masks_[row][col] &= ~(7 << 7);
		masks_[row][col] |=  (1 << 7 | mask << 8);
	}
	
	/**
	 * Return true iff the F mask has been set with a previous call to fMaskSet.
	 */
	bool isFMaskSet(
		size_t row,          // current row
		size_t col) const    // current column
	{
		return (masks_[row][col] & (1 << 10)) != 0;
	}

	/**
	 * Set the given cell's F mask.  This is the mask of remaining legal ways to
	 * backtrack from the F cell at this coordinate.  It's 2 bits long and has
	 * offset=11 into the 16-bit field.
	 */
	void fMaskSet(
		size_t row,          // current row
		size_t col,          // current column
		int mask)
	{
		assert_lt(mask, 4);
		masks_[row][col] &= ~(7 << 10);
		masks_[row][col] |=  (1 << 10 | mask << 11);
	}

	/**
	 * Analyze a cell in the SSE-filled dynamic programming matrix.  Determine &
	 * memorize ways that we can backtrack from the cell.  If there is at least one
	 * way to backtrack, select one at random and return the selection.
	 *
	 * There are a few subtleties to keep in mind regarding which cells can be at
	 * the end of a backtrace.  First of all: cells from which we can backtrack
	 * should not be at the end of a backtrace.  But have to distinguish between
	 * cells whose masks eventually become 0 (we shouldn't end at those), from
	 * those whose masks were 0 all along (we can end at those).
	 */
	void analyzeCell(
		size_t row,          // current row
		size_t col,          // current column
		size_t ct,           // current cell type: E/F/H
		int refc,
		int readc,
		int readq,
		const Scoring& sc,   // scoring scheme
		int64_t offsetsc,    // offset to add to each score
		RandomSource& rand,  // rand gen for choosing among equal options
		bool& empty,         // out: =true iff no way to backtrace
		int& cur,            // out: =type of transition
		bool& branch,        // out: =true iff we chose among >1 options
		bool& canMoveThru,   // out: =true iff ...
		bool& reportedThru); // out: =true iff ...

	/**
	 * Initialize the matrix of masks and backtracking flags.
	 */
	void initMasks() {
		assert_gt(nrow_, 0);
		assert_gt(ncol_, 0);
		masks_.resize(nrow_);
		reset_.resizeNoCopy(nrow_);
		reset_.fill(false);
	}

	/**
	 * Return the number of rows in the dynamic programming matrix.
	 */
	size_t nrow() const {
		return nrow_;
	}

	/**
	 * Return the number of columns in the dynamic programming matrix.
	 */
	size_t ncol() const {
		return ncol_;
	}
	
	/**
	 * Prepare a row so we can use it to store masks.
	 */
	void resetRow(size_t i) {
		assert(!reset_[i]);
		masks_[i].resizeNoCopy(ncol_);
		masks_[i].fillZero();
		reset_[i] = true;
	}

	static constexpr uint16_t get_wperv() {
		const uint16_t NWORDS_PER_REG = i16 ? EEI16_NWORDS_PER_REG : EEU8_NWORDS_PER_REG;
		return NWORDS_PER_REG;
	}

	static constexpr uint16_t get_vecPerCol(uint16_t nrows) {
		constexpr uint16_t  wperv = get_wperv();
		return (nrows + (wperv-1)) / wperv;
	}

	static constexpr uint32_t get_matbuf_size() {
		const uint32_t maxVecPerCol = get_vecPerCol(max_rows);
		// The +1 is so that we don't have to special-case the final column;
		// instead, we just write off the end of the useful part of the table
		// with pvEStore.
		return (max_cols+1) * (nvecPerCell_ * maxVecPerCol); // group last 2 to enforce uint32_t
	}

	typedef EList_sse<get_matbuf_size()>		EList_mb;

	bool             inited_;      // initialized?
	size_t           nrow_;        // # rows
	size_t           ncol_;        // # columns
	size_t           nvecrow_;     // # vector rows (<= nrow_)
	size_t           nveccol_;     // # vector columns (<= ncol_)
	static constexpr uint16_t  wperv_ = get_wperv();       // # words per vector
	size_t           nvecPerCol_;  // # vectors per column
	size_t           colstride_;   // # vectors b/t adjacent cells in same row
	size_t           rowstride_;   // # vectors b/t adjacent cells in same col
	EList_mb         matbuf_;      // buffer for holding vectors
	ELList<uint16_t> masks_;       // buffer for masks/backtracking flags
	EList<bool>      reset_;       // true iff row in masks_ has been reset
};

#define ALPHA_SIZE 5

/**
 * All the data associated with the query profile and other data needed for SSE
 * alignment of a query.
 */
template<bool i16, uint16_t max_rows = SSE_MAX_ROWS, uint16_t max_cols = SSE_MAX_COLS>
struct SSEData {
	SSEData(int cat = 0) { }

	static constexpr uint32_t get_sse_seglen(uint32_t rdlen) {
		const uint32_t NWORDS_PER_REG = i16 ? EEI16_NWORDS_PER_REG : EEU8_NWORDS_PER_REG;
		const uint32_t seglen = (rdlen + (NWORDS_PER_REG-1)) / NWORDS_PER_REG;
		return seglen;
	}

	static constexpr uint32_t get_nsses(uint32_t rdlen) {
		const uint32_t seglen = get_sse_seglen(rdlen);
		// How many SSERegI's are needed
		uint32_t nsses =
			64 +                    // slack bytes, for alignment?
			(seglen * ALPHA_SIZE)   // query profile data
			* 2;                    // & gap barrier data
		return nsses;
	}

	typedef EList_sse<get_nsses(max_rows)>		EList_pb;
	typedef SSEMatrix<i16,max_rows,max_cols>	SDMatrix;

	EList_pb       profbuf_;     // buffer for query profile & temp vecs
	size_t         qprofStride_; // stride for query profile
	size_t         gbarStride_;  // gap barrier for query profile
	SDMatrix       mat_;         // SSE matrix for holding all E, F, H vectors
	size_t         maxPen_;      // biggest penalty of all
	size_t         maxBonus_;    // biggest bonus of all
	size_t         lastIter_;    // which N-bit striped word has final row?
	size_t         lastWord_;    // which word within N-word has final row?
	int            bias_;        // all scores shifted up by this for unsigned
};

#define ROWSTRIDE_2COL 4
#define ROWSTRIDE 4
#define ROWSTRIDE_LOG2 2

#endif /*ndef ALIGNER_SWSSE_H_*/
