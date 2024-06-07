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

	QVal(const QKey &k, TIndexOffU elts)
	: k_(k), rangen_(1), eltn_(elts) {}

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const QVal& o) const {
		return (k_ == o.k_) && (rangen_ == o.rangen_) && (eltn_ == o.eltn_);
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
	const QKey& key() const { return k_; }

	/**
	 * Return the number of reference substrings associated with a read
	 * substring.
	 */
	uint8_t numRanges() const {
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
	bool valid() const { return rangen_ != 0; }

	/**
	 * Reset to invalid state.
	 */
	void reset() {
		rangen_ = eltn_ = 0;
	}

	/**
	 * Initialize Qval.
	 */
	void init(const QKey &k, TIndexOffU elts) {
		k_ = k; rangen_ = 1; eltn_ = elts;
	}

#ifndef NDEBUG
	/**
	 * Check that this QVal is internally consistent and consistent
	 * with the contents of the given cache.
	 */
	bool repOk(const AlignmentCache& ac) const;
#endif

protected:

	QKey       k_;      // idx of first elt in qlist
	uint8_t    rangen_; // # ranges (= # associated reference substrings)
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

	SAVal() : topf(), i(), len(OFF_MASK) { }

	/**
	 * True -> my key is equal to the given key.
	 */
	bool operator==(const SAVal& o) const {
		return  (topf == o.topf) &&
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
		TIndexOffU ii,
		TIndexOffU ln)
	{
		topf = tf;
		i = ii;
		len = ln;
	}

	TIndexOffU topf;  // top in BWT
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

	SATuple(SAKey k, TIndexOffU tf, TSlice o) {
		init(k, tf, o);
	}

	void init(SAKey k, TIndexOffU tf, TSlice o) {
		key = k; topf = tf; offs = o;
	}

	/**
	 * Initialize this SATuple from a subrange of the SATuple 'src'.
	 */
	void init(const SATuple& src, size_t first, size_t last) {
		key = src.key;
		topf = (TIndexOffU)(src.topf + first);
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
		return key == o.key && topf == o.topf && offs == o.offs;
	}

	void reset() { topf = OFF_MASK; offs.reset(); }

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
	TSlice   offs; // offsets
};

class AlignmentCacheInterface;

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
	friend class AlignmentCacheInterface;

public:

	AlignmentCache():
		samap_(CA_CAT),
		salist_size_(0),
		salist_offs_(0),
		version_(0) { }

	AlignmentCache(const AlignmentCache& other) = delete;
	AlignmentCache& operator=(const AlignmentCache& other) = delete;

	void setSAOffset(size_t off) {salist_offs_=off;}

	/**
	 * Add a new association between a read sequnce ('seq') and a
	 * reference sequence ('')
	 */
	void addOnTheFly(
		QVal& qv,         // qval that points to the range of reference substrings
		const SAKey& sak, // the key holding the reference substring
		TIndexOffU topf,    // top range elt in BWT index
		TIndexOffU botf) {   // bottom range elt in BWT index
		qv.init(sak, botf-topf);
		bool added = !samap_.contains(sak);
		if(added) {
			SAVal sav;
			sav.init(topf, salist_size_, botf - topf);
			salist_size_ += sav.len; // just remember how much we need for now	
			assert(sav.repOk(*this));
			samap_.insert(sak,sav);
		}
	}

	/**
	 * Clear the cache, i.e. turn it over.  All HitGens referring to
	 * ranges in this cache will become invalid and the corresponding
	 * reads will have to be re-aligned.
	 */
	void clear() {
		samap_.clear();
		salist_size_=0;
		salist_offs_=0;
		version_++;
	}

	/**
	 * Return the number of keys in the suffix array multimap.
	 */
	inline size_t saNumKeys() const { return samap_.size(); }

	/**
	 * Return the number of elements in the SA range list.
	 */
	inline size_t saSize() const { return salist_size_; }

	/**
	 * Return the current "version" of the cache, i.e. the total number
	 * of times it has turned over since its creation.
	 */
	inline uint32_t version() const { return version_; }

protected:

	ESimpleMap<SAKey, SAVal>  samap_;  // map from reference substrings to SA ranges
	size_t                 salist_size_; // local salist size
	size_t                 salist_offs_; // global offset

	uint32_t version_; // cache version

private:

	/**
	 * Given a QVal, populate the given EList of SATuples with records
	 * describing all of the cached information about the QVal's
	 * reference substrings.
	 */
	template <int S1, int S2>
	inline void queryQvalImpl(
		const QVal& qv,
		EList<TIndexOffU, S1>& salist,
		EList<SATuple, S2>&    satups,
		size_t& nrange,
		size_t& nelt) const
	{
		assert(qv.repOk(*this));
		if (qv.numRanges()>0) {
			// assume only one element
			// Get corresponding SAKey, containing similar reference
			// sequence & length
			const SAKey& sak = qv.key();
			// Get corresponding SAVal
			assert(samap_.lookup(sak) != NULL);
			const SAVal& sav = *(samap_.lookup(sak));
			assert(sav.repOk(*this));
			if(sav.len > 0) {
				nrange++;
				satups.expand();
				TSlice offs = TSlice(salist, salist_offs_+sav.i, sav.len);
				offs.fill(OFF_MASK);
				satups.back().init(sak, sav.topf, offs);
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

};

/**
 * Implements the cache lookup
 **/
class AlignmentCacheInterface {
public:
	AlignmentCacheInterface(const AlignmentCache& cache, EList<TIndexOffU>& salist)
		:cache_(cache) , salist_(salist) {}

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
		size_t& nelt)
	{
		cache_.queryQvalImpl(qv, salist_, satups, nrange, nelt);
	}
private:
	const AlignmentCache& cache_;
	EList<TIndexOffU>&    salist_;
};

/**
 * Handle a set of caches
 **/
class AlignmentMultiCache {
public:
	AlignmentMultiCache(const uint32_t ncaches)
		: caches_(new AlignmentCache[ncaches])
		, ncaches_(ncaches)
		{}

	~AlignmentMultiCache() {
		delete[] caches_;
	}

	inline uint32_t n_els() const {return ncaches_;}

	inline AlignmentCache &getCache(uint32_t idx) { return caches_[idx];}
	inline const AlignmentCache &getCache(uint32_t idx) const { return caches_[idx];}

	inline AlignmentCacheInterface getCacheInterface(uint32_t idx) {
		AlignmentCacheInterface out(caches_[idx], salists_);
		return out;
	}

	/**
	 * Finalize the cache, initialize salist_
	 **/
	void finalize() {
		size_t nsas = 0;
		for (uint32_t i=0; i<ncaches_; i++) {
			AlignmentCache& cache = caches_[i];
			cache.setSAOffset(nsas);
			nsas += cache.saSize();
		}
		salists_.resizeNoCopy(nsas);
	}

private:
	AlignmentCache   *caches_;
	EList<TIndexOffU> salists_;
	const uint32_t    ncaches_;  // size of caches_
};

#endif /*ALIGNER_CACHE_H_*/
