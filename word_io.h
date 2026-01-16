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

#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "assert_helpers.h"
#include "btypes.h"

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeU(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

/**
 * Write a 32/64 bit unsigned to an output stream using the native
 * endianness.
 */
template <typename T>
static inline void writeI(std::ostream& out, T x) {
	out.write((const char*)&x, sizeof(T));
}

/**
 * Read a 32/64 bit unsigned from an input stream
 */


template <typename T>
static inline T readU(std::istream& in) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	return x;
}


/**
 * Read a 32/64 bit unsigned from a FILE*
 */
template <typename T>
static inline T readU(FILE* in) {
	T x;
	if(fread((void *)&x, sizeof(T), 1, in) != 1) {
		perror("readU");
		exit(1);
	}
	return x;
}


/**
 * Read a 32/64 bit signed from an input stream
 */
template <typename T>
static inline T readI(std::istream& in) {
	T x;
	in.read((char *)&x, sizeof(T));
	assert_eq(sizeof(T), in.gcount());
	return x;
}


/**
 * Read a 32/64 bit unsigned from a FILE*
 */
template <typename T>
static inline T readI(FILE* in) {
	T x;
	if(fread((void *)&x, sizeof(T), 1, in) != 1) {
		perror("readI");
		exit(1);
	}
	return x;
}



#endif /*WORD_IO_H_*/
