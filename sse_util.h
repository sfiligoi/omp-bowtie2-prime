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

#ifndef SSE_UTIL_H_
#define SSE_UTIL_H_

#include "assert_helpers.h"
#include "ds.h"
#include "limit.h"
#include <iostream>
#include "sse_wrap.h"

class EList_sse : public BTAllocatorHandler<SSERegI> {
public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList_sse  (int cat = 0) :
		cat_(cat), last_alloc_(NULL), list_(NULL), sz_(0), cur_(0)
	{
		assert_geq(cat, 0);
	}

	/**
	 * Destructor.
	 */
	~EList_sse  () { free(); }

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return number of elements allocated.
	 */
	inline size_t capacity() const { return sz_; }
	
	/**
	 * Ensure that there is sufficient capacity to expand to include
	 * 'thresh' more elements without having to expand.
	 */
	inline void ensure(size_t thresh) {
		if(list_ == NULL) lazyInit();
		expandCopy(cur_ + thresh);
	}

	/**
	 * Ensure that there is sufficient capacity to include 'newsz' elements.
	 * If there isn't enough capacity right now, expand capacity to exactly
	 * equal 'newsz'.
	 */
	inline void reserveExact(size_t newsz) {
		if(list_ == NULL) lazyInitExact(newsz);
		expandCopyExact(newsz);
	}

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }
	
	/**
	 * Return true iff list hasn't been initialized yet.
	 */
	inline bool null() const { return list_ == NULL; }

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandCopy(sz);
		}
		cur_ = sz;
	}
	
	/**
	 * Zero out contents of vector.
	 */
	void zero() {
		if(cur_ > 0) {
			memset(list_, 0, cur_ * sizeof(SSERegI));
		}
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.  Do not copy the elements over.
	 */
	void resizeNoCopy(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandNoCopy(sz);
		}
		cur_ = sz;
	}

	/**
	 * If size is less than requested size, resize up to exactly sz and set
	 * cur_ to requested sz.
	 */
	void resizeExact(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInitExact(sz);
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
	}

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline SSERegI& operator[](size_t i) {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline SSERegI operator[](size_t i) const {
		assert_lt(i, cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline SSERegI& get(size_t i) {
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	inline SSERegI get(size_t i) const {
		return operator[](i);
	}

	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	SSERegI *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const SSERegI *ptr() const { return list_; }

	/**
	 * Return memory category.
	 */
	int cat() const { return cat_; }

private:

	/**
	 * Initialize memory for EList.
	 */
	void lazyInit() {
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}

	/**
	 * Initialize exactly the prescribed number of elements for EList.
	 */
	void lazyInitExact(size_t sz) {
		assert_gt(sz, 0);
		assert(list_ == NULL);
		sz_ = sz;
		list_ = alloc(sz);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	SSERegI *alloc(size_t sz) {
		SSERegI* last_alloc_;
		last_alloc_ = this->allocate(sz + 4);
                this->last_alloc_ = last_alloc_;
		SSERegI* tmp = last_alloc_;
		size_t tmpint = (size_t)tmp;
		// Align it!
		const size_t alignmask = NBYTES_PER_REG-1;
		if((tmpint & alignmask) != 0) {
			tmpint += alignmask;
			tmpint &= (~alignmask);
			tmp = reinterpret_cast<SSERegI*>(tmpint);
		}
		assert_eq(0, (tmpint & alignmask)); // should be NBYTES_PER_REG-byte aligned
		assert(tmp != NULL);
#ifdef USE_MEM_TALLY
		gMemTally.add(cat_, sz);
#endif
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			deallocate(last_alloc_, sz_ + 4);
#ifdef USE_MEM_TALLY
			gMemTally.del(cat_, sz_);
#endif
			list_ = NULL;
			sz_ = cur_ = 0;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.  Size
	 * increases quadratically with number of expansions.  Copy old contents
	 * into new buffer using operator=.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Copy
	 * old contents into new buffer using operator=.
	 */
	void expandCopyExact(size_t newsz) {
		if(newsz <= sz_) return;
                SSERegI* prev_last_alloc = last_alloc_;
		SSERegI* tmp = alloc(newsz);
		assert(tmp != NULL);
		size_t cur = cur_;
		if(list_ != NULL) {
 			for(size_t i = 0; i < cur_; i++) {
				// Note: operator= is used
				tmp[i] = list_[i];
			}
                        SSERegI* current_last_alloc = last_alloc_;
                        last_alloc_ = prev_last_alloc;
			free();
                        last_alloc_ = current_last_alloc;
		}
		list_ = tmp;
		sz_ = newsz;
		cur_ = cur;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Size increases quadratically with number of expansions.  Don't copy old
	 * contents into the new buffer.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandNoCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Don't
	 * copy old contents into the new buffer.
	 */
	void expandNoCopyExact(size_t newsz) {
		assert(list_ != NULL);
		assert_gt(newsz, 0);
		free();
		SSERegI* tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert_gt(sz_, 0);
	}

	int      cat_;        // memory category, for accounting purposes
	SSERegI* last_alloc_; // what new[] originally returns
	SSERegI *list_;       // list ptr, aligned version of what new[] returns
	size_t   sz_;         // capacity
	size_t   cur_;        // occupancy (AKA size)
};

struct  CpQuad {
	CpQuad() { reset(); }
	
	void reset() { sc[0] = sc[1] = sc[2] = sc[3] = 0; }
	
	bool operator==(const CpQuad& o) const {
		return sc[0] == o.sc[0] &&
		       sc[1] == o.sc[1] &&
			   sc[2] == o.sc[2] &&
			   sc[3] == o.sc[3];
	}

	int16_t sc[4];
};

#endif
