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

#ifndef ALIGNER_CACHE_H_
#define ALIGNER_CACHE_H_

/**
 * CACHEING
 *
 * By caching the results of some alignment sub-problems, we hope to
 * enable a "fast path" for read alignment whereby answers are mostly
 * looked up rather than calculated from scratch.  This is particularly
 * effective when the input is sorted or otherwise grouped in a way
 * that brings together reads with (at least some) seed sequences in
 * common.
 *
 * But the cache is also where results are held, regardless of whether
 * the results are maintained & re-used across reads.
 *
 * The cache consists of two linked potions:
 *
 * 1. A multimap from seed strings (i.e. read substrings) to reference strings
 *    that are within some edit distance (roughly speaking).  This is the "seed
 *    multimap".
 *
 *    Key:   Read substring (2-bit-per-base encoded + length)
 *    Value: Set of reference substrings (i.e. keys into the suffix
 *           array multimap).
 *
 * 2. A multimap from reference strings to the corresponding elements of the
 *    suffix array.  Elements are filled in with reference-offset info as it's
 *    calculated.  This is the "suffix array multimap"
 *
 *    Key:   Reference substring (2-bit-per-base encoded + length)
 *    Value: (a) top from BWT, (b) length of range, (c) offset of first
 *           range element in
 *
 * For both multimaps, we use a combo Red-Black tree and EList.  The payload in
 * the Red-Black tree nodes points to a range in the EList.
 */

#include <iostream>
#include "ds.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "btypes.h"

#define CACHE_PAGE_SZ (16 * 1024)

class BwtTopBot {
public:
	BwtTopBot() : topf(0), botf(0), topb(0), botb(0) {}

	BwtTopBot(
		TIndexOffU _topf,        // top in BWT
		TIndexOffU _botf,        // bot in BWT
		TIndexOffU _topb,        // top in BWT'
		TIndexOffU _botb)        // bot in BWT'
	: topf(_topf)
	, botf(_botf)
	, topb(_topb)
	, botb(_botb)
	{}

	void set(
		TIndexOffU _topf,        // top in BWT
		TIndexOffU _botf,        // bot in BWT
		TIndexOffU _topb,        // top in BWT'
		TIndexOffU _botb)        // bot in BWT'
	{
		topf = _topf;
		botf = _botf;
		topb = _topb;
		botb = _botb;
	}

	TIndexOffU topf;        // top in BWT
	TIndexOffU botf;        // bot in BWT
	TIndexOffU topb;        // top in BWT'
	TIndexOffU botb;        // bot in BWT'
};

class BwtTopBotFw {
public:
	BwtTopBotFw() : topf(0), botf(0) {}

	BwtTopBotFw(
		TIndexOffU _topf,        // top in BWT
		TIndexOffU _botf)        // bot in BWT
	: topf(_topf)
	, botf(_botf)
	{}

	void set(
		TIndexOffU _topf,        // top in BWT
		TIndexOffU _botf)        // bot in BWT
	{
		topf = _topf;
		botf = _botf;
	}

	TIndexOffU topf;        // top in BWT
	TIndexOffU botf;        // bot in BWT
};

typedef EListSlice<TIndexOffU,128> TSlice;

/**
 * Key for the query multimap: the read substring and its length.
 */
struct QKey {

	/**
	 * Initialize invalid QKey.
	 */
	QKey() { reset(); }

	/**
	 * Initialize QKey from DNA string.
	 */
	QKey(const BTDnaString& s ASSERT_ONLY(, BTDnaString& tmp)) {
		init(s ASSERT_ONLY(, tmp));
	}

	QKey(const char * s, const uint32_t l ASSERT_ONLY(, BTDnaString& tmp)) {
		init(s, l ASSERT_ONLY(, tmp));
	}

	/**
	 * Initialize QKey from DNA string.  Rightmost character is placed in the
	 * least significant bitpair.
	 */
	bool init(
		const char *   s,
		const uint32_t l
		ASSERT_ONLY(, BTDnaString& tmp))
	{
		seq = 0;
		len = l;
		ASSERT_ONLY(tmp.clear());
		if(len > 32) {
			len = 0xffffffff;
			return false; // wasn't cacheable
		} else {
			// Rightmost char of 's' goes in the least significant bitpair
			for(uint32_t i = 0; i < 32 && i < l; i++) {
				int c = (int)s[i];
				assert_range(0, 4, c);
				if(c == 4) {
					len = 0xffffffff;
					return false;
				}
				seq = (seq << 2) | s[i];
			}
			ASSERT_ONLY(toString(tmp));
			assert(sstr_eq(tmp, s));
			assert_leq(len, 32);
			return true; // was cacheable
		}
	}

	bool init(
		const BTDnaString& s
		ASSERT_ONLY(, BTDnaString& tmp))
	{
		init(s.buf(), (uint32_t)s.length() ASSERT_ONLY(, tmp));
	}

	/**
	 * Convert this key to a DNA string.
	 */
	void toString(BTDnaString& s) {
		s.resize(len);
		uint64_t sq = seq;
		for(int i = (len)-1; i >= 0; i--) {
			s.set((uint32_t)(sq & 3), i);
			sq >>= 2;
		}
	}

	/**
	 * Return true iff the read substring is cacheable.
	 */
	bool cacheable() const { return len != 0xffffffff; }

	/**
	 * Reset to uninitialized state.
	 */
	void reset() { seq = 0; len = 0xffffffff; }

	/**
	 * True -> my key is less than the given key.
	 */
	bool operator<(const QKey& o) const {
		return seq < o.seq || (seq == o.seq && len < o.len);
	}

	/**
	 * True -> my key is greater than the given key.
	 */
	bool operator>(const QKey& o) const {
		return !(*this < o || *this == o);
	}

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const QKey& o) const {
		return seq == o.seq && len == o.len;
	}


	/**
	 * True -> my key is not equal to the given key.
	 */
	bool operator!=(const QKey& o) const {
		return !(*this == o);
	}

#ifndef NDEBUG
	/**
	 * Check that this is a valid, initialized QKey.
	 */
	bool repOk() const {
		return len != 0xffffffff;
	}
#endif

	uint64_t seq; // sequence
	uint32_t len; // length of sequence
};

class AlignmentCache;

/**
 * Payload for the query multimap: a range of elements in the reference
 * string list.
 */
class QVal {

public:

	QVal() { reset(); }

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const QVal& o) const {
		return (i_ == o.i_) && (rangen_ == o.rangen_) && (eltn_ == o.eltn_);
	}

	/**
	 * True -> my key is not equal to the given key.
	 */
	bool operator!=(const QVal& o) const {
		return !(*this == o);
	}

	/**
	 * Return the offset of the first reference substring in the qlist.
	 */
	TIndexOffU offset() const { return i_; }

	/**
	 * Return the number of reference substrings associated with a read
	 * substring.
	 */
	TIndexOffU numRanges() const {
		assert(valid());
		return rangen_;
	}

	/**
	 * Return the number of elements associated with all associated
	 * reference substrings.
	 */
	TIndexOffU numElts() const {
		assert(valid());
		return eltn_;
	}

	/**
	 * Return true iff the read substring is not associated with any
	 * reference substrings.
	 */
	bool empty() const {
		assert(valid());
		return numRanges() == 0;
	}

	/**
	 * Return true iff the QVal is valid.
	 */
	bool valid() const { return rangen_ != OFF_MASK; }

	/**
	 * Reset to invalid state.
	 */
	void reset() {
		i_ = 0; rangen_ = eltn_ = OFF_MASK;
	}

	/**
	 * Initialize Qval.
	 */
	void init(TIndexOffU i, TIndexOffU ranges, TIndexOffU elts) {
		i_ = i; rangen_ = ranges; eltn_ = elts;
	}

	/**
	 * Tally another range with given number of elements.
	 */
	void addRange(TIndexOffU numElts) {
		rangen_++;
		eltn_ += numElts;
	}

#ifndef NDEBUG
	/**
	 * Check that this QVal is internally consistent and consistent
	 * with the contents of the given cache.
	 */
	bool repOk(const AlignmentCache& ac) const;
#endif

protected:

	TIndexOffU i_;      // idx of first elt in qlist
	TIndexOffU rangen_; // # ranges (= # associated reference substrings)
	TIndexOffU eltn_;   // # elements (total)
};

/**
 * Key for the suffix array multimap: the reference substring and its
 * length.  Same as QKey so I typedef it.
 */
typedef QKey SAKey;

/**
 * Payload for the suffix array multimap: (a) the top element of the
 * range in BWT, (b) the offset of the first elt in the salist, (c)
 * length of the range.
 */
struct SAVal {

	SAVal() : topf(), topb(), i(), len(OFF_MASK) { }

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const SAVal& o) const {
		return  (topf == o.topf) && (topb == o.topb) && 
			(i == o.i) && (len == o.len);
	}

	/**
	 * True -> my key is not equal to the given key.
	 */
	bool operator!=(const SAVal& o) const {
		return !(*this == o);
	}

	/**
	 * Return true iff the SAVal is valid.
	 */
	bool valid() { return len != OFF_MASK; }

#ifndef NDEBUG
	/**
	 * Check that this SAVal is internally consistent and consistent
	 * with the contents of the given cache.
	 */
	bool repOk(const AlignmentCache& ac) const;
#endif

	/**
	 * Initialize the SAVal.
	 */
	void init(
		TIndexOffU tf,
		TIndexOffU tb,
		TIndexOffU ii,
		TIndexOffU ln)
	{
		topf = tf;
		topb = tb;
		i = ii;
		len = ln;
	}

	TIndexOffU topf;  // top in BWT
	TIndexOffU topb;  // top in BWT'
	TIndexOffU i;     // idx of first elt in salist
	TIndexOffU len;   // length of range
};

/**
 * One data structure that encapsulates all of the cached information
 * associated with a particular reference substring.  This is useful
 * for summarizing what info should be added to the cache for a partial
 * alignment.
 */
class SATuple {

public:

	SATuple() { reset(); };

	SATuple(SAKey k, TIndexOffU tf, TIndexOffU tb, TSlice o) {
		init(k, tf, tb, o);
	}

	void init(SAKey k, TIndexOffU tf, TIndexOffU tb, TSlice o) {
		key = k; topf = tf; topb = tb; offs = o;
	}

	/**
	 * Initialize this SATuple from a subrange of the SATuple 'src'.
	 */
	void init(const SATuple& src, size_t first, size_t last) {
		assert_neq(OFF_MASK, src.topb);
		key = src.key;
		topf = (TIndexOffU)(src.topf + first);
		topb = OFF_MASK; // unknown!
		offs.init(src.offs, first, last);
	}

#ifndef NDEBUG
	/**
	 * Check that this SATuple is internally consistent and that its
	 * PListSlice is consistent with its backing PList.
	 */
	bool repOk() const {
		assert(offs.repOk());
		return true;
	}
#endif

	/**
	 * Function for ordering SATuples.  This is used when prioritizing which to
	 * explore first when extending seed hits into full alignments.  Smaller
	 * ranges get higher priority and we use 'top' to break ties, though any
	 * way of breaking a tie would be fine.
	 */
	bool operator<(const SATuple& o) const {
		if(offs.size() < o.offs.size()) {
			return true;
		}
		if(offs.size() > o.offs.size()) {
			return false;
		}
		return topf < o.topf;
	}
	bool operator>(const SATuple& o) const {
		if(offs.size() < o.offs.size()) {
			return false;
		}
		if(offs.size() > o.offs.size()) {
			return true;
		}
		return topf > o.topf;
	}

	bool operator==(const SATuple& o) const {
		return key == o.key && topf == o.topf && topb == o.topb && offs == o.offs;
	}

	void reset() { topf = topb = OFF_MASK; offs.reset(); }

	/**
	 * Set the length to be at most the original length.
	 */
	void setLength(size_t nlen) {
		assert_leq(nlen, offs.size());
		offs.setLength(nlen);
	}

	/**
	 * Return the number of times this reference substring occurs in the
	 * reference, which is also the size of the 'offs' TSlice.
	 */
	size_t size() const { return offs.size(); }

	// bot/length of SA range equals offs.size()
	SAKey    key;  // sequence key
	TIndexOffU topf;  // top in BWT index
	TIndexOffU topb;  // top in BWT' index
	TSlice   offs; // offsets
};

/**
 * Encapsulate the data structures and routines that constitute a
 * particular cache, i.e., a particular stratum of the cache system,
 * which might comprise many strata.
 *
 * Each thread has a "current-read" AlignmentCache which is used to
 * build and store subproblem results as alignment is performed.  When
 * we're finished with a read, we might copy the cached results for
 * that read (and perhaps a bundle of other recently-aligned reads) to
 * a higher-level "across-read" cache.  Higher-level caches may or may
 * not be shared among threads.
 *
 * A cache consists chiefly of two multimaps, each implemented as a
 * Red-Black tree map backed by an EList.  A 'version' counter is
 * incremented every time the cache is cleared.
 *
 * TODO: Doucment the changes
 */
class AlignmentCache {

	typedef EList<SAKey> TQList;
	typedef EList<TIndexOffU> TSAList;

public:

	AlignmentCache():
		qlist_(CA_CAT),
		samap_(CA_CAT),
		salist_(CA_CAT),
		version_(0) { }

	AlignmentCache(const AlignmentCache& other) = delete;
	AlignmentCache& operator=(const AlignmentCache& other) = delete;

	/**
	 * Given a QVal, populate the given EList of SATuples with records
	 * describing all of the cached information about the QVal's
	 * reference substrings.
	 */
	template <int S>
	void queryQval(
		const QVal& qv,
		EList<SATuple, S>& satups,
		size_t& nrange,
		size_t& nelt,
		bool getLock = true)
	{
			queryQvalImpl(qv, satups, nrange, nelt);
	}

	/**
	 * Add a new association between a read sequnce ('seq') and a
	 * reference sequence ('')
	 */
	bool addOnTheFly(
		QVal& qv,         // qval that points to the range of reference substrings
		const SAKey& sak, // the key holding the reference substring
		TIndexOffU topf,    // top range elt in BWT index
		TIndexOffU botf,    // bottom range elt in BWT index
		bool getLock = true);

	/**
	 * Clear the cache, i.e. turn it over.  All HitGens referring to
	 * ranges in this cache will become invalid and the corresponding
	 * reads will have to be re-aligned.
	 */
	void clear() {
		qlist_.clear();
		samap_.clear();
		salist_.clear();
		version_++;
	}

	/**
	 * Return the number of keys in the suffix array multimap.
	 */
	size_t saNumKeys() const { return samap_.size(); }

	/**
	 * Return the number of elements in the reference substring list.
	 */
	size_t qSize() const { return qlist_.size(); }

	/**
	 * Return the number of elements in the SA range list.
	 */
	size_t saSize() const { return salist_.size(); }

	/**
	 * Return the current "version" of the cache, i.e. the total number
	 * of times it has turned over since its creation.
	 */
	uint32_t version() const { return version_; }

protected:

	TQList                 qlist_;  // list of reference substrings
	ESimpleMap<SAKey, SAVal>  samap_;  // map from reference substrings to SA ranges
	TSAList                salist_; // list of SA ranges

	uint32_t version_; // cache version

private:

	template <int S>
	void queryQvalImpl(
		const QVal& qv,
		EList<SATuple, S>& satups,
		size_t& nrange,
		size_t& nelt)
	{
		assert(qv.repOk(*this));
		const size_t refi = qv.offset();
		const size_t reff = refi + qv.numRanges();
		// For each reference sequence sufficiently similar to the
		// query sequence in the QKey...
		for(size_t i = refi; i < reff; i++) {
			// Get corresponding SAKey, containing similar reference
			// sequence & length
			SAKey sak = qlist_.get(i);
			// Shouldn't have identical keys in qlist_
			assert(i == refi || qlist_.get(i) != qlist_.get(i-1));
			// Get corresponding SAVal
			assert(samap_.lookup(sak) != NULL);
			const SAVal& sav = *(samap_.lookup(sak));
			assert(sav.repOk(*this));
			if(sav.len > 0) {
				nrange++;
				satups.expand();
				satups.back().init(sak, sav.topf, sav.topb, TSlice(salist_, sav.i, sav.len));
				nelt += sav.len;
#ifndef NDEBUG
				// Shouldn't add consecutive identical entries too satups
				if(i > refi) {
					const SATuple b1 = satups.back();
					const SATuple b2 = satups[satups.size()-2];
					assert(b1.key != b2.key || b1.topf != b2.topf || b1.offs != b2.offs);
				}
#endif
			}
		}
	}

	/**
	 * Add a new association between a read sequnce ('seq') and a
	 * reference sequence ('')
	 */
	bool addOnTheFlyImpl(
		QVal& qv,         // qval that points to the range of reference substrings
		const SAKey& sak, // the key holding the reference substring
		TIndexOffU topf,    // top range elt in BWT index
		TIndexOffU botf);    // bottom range elt in BWT index

};

/**
 * Interface used to query and update a pair of caches: one thread-
 * local and unsynchronized, another shared and synchronized.  One or
 * both can be NULL.
 */
class AlignmentCacheIface {

public:

	AlignmentCacheIface(
		AlignmentCache *current):
		rangen_(0),
		eltsn_(0),
		current_(current)
	{
		assert(current_ != NULL);
	}

	~AlignmentCacheIface() {}

	AlignmentCacheIface(const AlignmentCacheIface& other) = delete;
	AlignmentCacheIface& operator=(const AlignmentCacheIface& other) = delete;

	ASSERT_ONLY(BTDnaString tmpdnastr_);

	const QVal& align(
		const SAKey& sak, // the key holding the reference substring
		TIndexOffU topf,              // top in BWT index
		TIndexOffU botf) {            // bot in BWT index
		//beginAlign;
		qvbuf_.reset();
		addOnTheFly(sak, topf, botf);
		// finishAlign
		resetRead();
		return qvbuf_;
	}

	/**
	 * A call to this member indicates that the caller has finished
	 * with the last read (if any) and is ready to work on the next.
	 * This gives the cache a chance to reset some of its state if
	 * necessary.
	 */
	void nextRead() {
		current_->clear();
		resetRead();
		assert(!aligning());
	}


	/**
	 * Clears both the local and shared caches.
	 */
	void clear() {
		if(current_ != NULL) current_->clear();
	}

	/**
	 * Add an alignment to the running list of alignments being
	 * compiled for the current read in the local cache.
	 */
	bool addOnTheFly(
		const SAKey& sak, // the key holding the reference substring
		TIndexOffU topf,            // top in BWT index
		TIndexOffU botf,            // bot in BWT index
		bool getLock = true)      // true -> lock is not held by caller
	{

		assert(aligning());
		assert(repOk());
		//assert(sak.cacheable());
		if(current_->addOnTheFly(qvbuf_, sak, topf, botf, getLock)) {
			rangen_++;
			eltsn_ += (botf-topf);
			return true;
		}
		return false;
	}

	/**
	 * Given a QVal, populate the given EList of SATuples with records
	 * describing all of the cached information about the QVal's
	 * reference substrings.
	 */
	template<int S>
	void queryQval(
		const QVal& qv,
		EList<SATuple, S>& satups,
		size_t& nrange,
		size_t& nelt,
		bool getLock = true)
	{
		current_->queryQval(qv, satups, nrange, nelt, getLock);
	}

	/**
	 * Return a pointer to the current-read cache object.
	 */
	const AlignmentCache* currentCache() const { return current_; }

	size_t curNumRanges() const { return rangen_; }
	size_t curNumElts()   const { return eltsn_;  }

#ifndef NDEBUG
	/**
	 * Check that AlignmentCacheIface is internally consistent.
	 */
	bool repOk() const {
		assert(current_ != NULL);
		assert_geq(eltsn_, rangen_);
		return true;
	}
#endif

	/**
	 * Return the alignment cache for the current read.
	 */
	const AlignmentCache& current() {
		return *current_;
	}

protected:

	/**
	 * Reset fields encoding info about the in-process read.
	 */
	void resetRead() {
		rangen_ = eltsn_ = 0;
	}

	QVal qvbuf_; // buffer for when key is uncacheable but we need a qv

	size_t rangen_; // number of ranges since last alignment job began
	size_t eltsn_;  // number of elements since last alignment job began

	AlignmentCache *current_; // cache dedicated to the current read
};

#endif /*ALIGNER_CACHE_H_*/
