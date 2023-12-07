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

#ifndef LIMIT_H_
#define LIMIT_H_

#include <stdint.h>
#include <cstring>

#define MIN_U8 0
#define MAX_U8 std::numeric_limits<uint8_t>::max()
#define MIN_U16 0
#define MAX_U16 std::numeric_limits<uint16_t>::max()
#define MIN_U32 0
#define MAX_U32 std::numeric_limits<uint32_t>::max()
#define MIN_U64 0
#define MAX_U64 std::numeric_limits<uint64_t>::max()
#define MIN_SIZE_T std::numeric_limits<size_t>::min()
#define MAX_SIZE_T std::numeric_limits<size_t>::max()

#define MIN_I std::numeric_limits<int>::min()
#define MAX_I std::numeric_limits<int>::max()
#define MIN_I8 std::numeric_limits<int8_t>::min()
#define MAX_I8 std::numeric_limits<int8_t>::max()
#define MIN_I16 std::numeric_limits<int16_t>::min()
#define MAX_I16 std::numeric_limits<int16_t>::max()
#define MIN_I32 std::numeric_limits<int32_t>::min()
#define MAX_I32 std::numeric_limits<int32_t>::max()
#define MIN_I64 std::numeric_limits<int64_t>::min()
#define MAX_I64 std::numeric_limits<int64_t>::max()

#endif
