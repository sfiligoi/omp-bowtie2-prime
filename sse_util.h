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

// Partially mimick EList, with static buffer of max size and always all elements used
template<uint32_t sz_>
class EList_sse {
public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList_sse  (int cat = 0)
	{
	}

	/**
	 * Destructor.
	 */
	~EList_sse  () {}

	/**
	 * Return number of elements.
	 */
	static constexpr size_t size() { return sz_; }

	/**
	 * Return number of elements allocated.
	 */
	static constexpr size_t capacity() { return sz_; }
	
	/**
	 * Return true iff there are no elements.
	 */
	static constexpr bool empty() { return false; }
	
	
	/**
	 * Zero out contents of vector.
	 */
	inline void zero() {
		memset(list_, 0, sz_ * sizeof(SSERegI));
	}

	/**
	 * Return a reference to the ith element.
	 */
	constexpr SSERegI& operator[](size_t i) {
		assert_lt(i, sz_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	constexpr SSERegI operator[](size_t i) const {
		assert_lt(i, sz_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	constexpr SSERegI& get(size_t i) {
		assert_lt(i, sz_);
		return operator[](i);
	}
	
	/**
	 * Return a reference to the ith element.
	 */
	constexpr SSERegI get(size_t i) const {
		assert_lt(i, sz_);
		return operator[](i);
	}

	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	constexpr SSERegI *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	constexpr const SSERegI *ptr() const { return list_; }

public:

	SSERegI  list_[sz_];  // list ptr, aligned version of what new[] returns
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
