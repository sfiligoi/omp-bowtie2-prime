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

// Mimick EList, but with static buffer of max size
template<uint32_t sz_>
class EList_sse {
public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList_sse  (int cat = 0) :
		cur_(0)
	{
	}

	/**
	 * Destructor.
	 */
	~EList_sse  () {}

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return number of elements allocated.
	 */
	static constexpr size_t capacity() { return sz_; }
	
	/**
	 * Ensure that there is sufficient capacity to expand to include
	 * 'thresh' more elements without having to expand.
	 */
	inline void ensure(size_t thresh) {
		assert_lt(cur_ + thresh, sz_);
	}

	/**
	 * Ensure that there is sufficient capacity to include 'newsz' elements.
	 * If there isn't enough capacity right now, expand capacity to exactly
	 * equal 'newsz'.
	 */
	inline void reserveExact(size_t newsz) {
		assert_leq(newsz, sz_);
	}

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }
	
	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		assert_leq(sz, sz_);
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
	void resizeNoCopy(size_t sz) { resize(sz); }

	/**
	 * If size is less than requested size, resize up to exactly sz and set
	 * cur_ to requested sz.
	 */
	void resizeExact(size_t sz) { resize(sz); }

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
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
	int cat() const {return this->cat_; }

private:

	SSERegI  list_[sz_];  // list ptr, aligned version of what new[] returns
	uint32_t cur_;        // occupancy (AKA size)
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
