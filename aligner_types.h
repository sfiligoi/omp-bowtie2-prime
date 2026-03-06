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

#ifndef ALIGNER_TYPES_H_
#define ALIGNER_TYPES_H_

#include "cstdint"

typedef int64_t TAlScore;

// Hardcode for now. May want to pass in Makefile
// The cosen values are appropriate for short reads alignment
// since rdlen is typically 150 (as of July 2024)
#define ALN_MAX_ROWS 160
#define ALN_MAX_COLS 200

#endif
