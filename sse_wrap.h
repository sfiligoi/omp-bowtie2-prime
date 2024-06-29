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

/*
 * sse_wrap.h
 *
 * Routines to wrap Streaming SIMD Extensions (SSE) emmintrin.h
 * for an Intel x86 CPU and SIMD Everywhere (simde) for other CPUs.
 */

#ifndef SSE_WRAP_H_
#define SSE_WRAP_H_

#ifdef SSE_AVX2
#include <immintrin.h>
#define NBYTES_PER_REG 32
#define BYTES_LOG2_PER_REG 5
#define SSE_MASK_ALL ((int) 0xffffffff)

typedef __m256i SSEReg;
typedef SSEReg  SSEMem;  // memory and register representation are the same

#define sse_load_siall(x) _mm256_load_si256(x)
#define sse_store_siall(x, y) _mm256_store_si256(x, y)

#define sse_set1_epi8(x) _mm256_set1_epi8(x)
#define sse_subs_epu8(x, y) _mm256_subs_epu8(x, y)
#define sse_setzero_siall() _mm256_setzero_si256()
#define sse_max_epu8(x, y) _mm256_max_epu8(x, y)
#define sse_or_siall(x, y) _mm256_or_si256(x, y)

#define sse_insert_epi16(x, y, z) _mm256_insert_epi16(x, y, z)
#define sse_cmpeq_epi8(x, y) _mm256_cmpeq_epi8(x, y)
#define sse_movemask_epi8(x) _mm256_movemask_epi8(x)
/* AVX2 does not have a native 256-bit shift instruction */
/* Note only works for y<=16, which is OK for this code */
#define sse_slli_siall(x, y) \
	_mm256_alignr_epi8(x, _mm256_permute2x128_si256(x, x, _MM_SHUFFLE(0, 0, 2, 0)), 16-y)

#ifdef ENABLE_I16
#define sse_adds_epi16(x, y) _mm256_adds_epi16(x, y)
#define sse_subs_epi16(x, y) _mm256_subs_epi16(x, y)
#define sse_set1_epi16(x) _mm256_set1_epi16(x)
#define sse_max_epi16(x, y) _mm256_max_epi16(x, y)
#define sse_slli_epi16(x, y) _mm256_slli_epi16(x, y)

#define sse_cmpeq_epi16(x, y) _mm256_cmpeq_epi16(x, y)
#define sse_cmpgt_epi16(x, y) _mm256_cmpgt_epi16(x, y)
#endif

#else /* no SSE_AVX2 */

#ifdef SSE_SW4

#define SSE_DISABLE 1

#define NBYTES_PER_REG 4
#define BYTES_LOG2_PER_REG 2
#define SSE_MASK_ALL 0xff

#include <algorithm>
#include <stdint.h>

typedef struct {
  int16_t el[NBYTES_PER_REG];

  // need to bound to emulate u8
  static inline uint8_t enforce_u8(int16_t one_el) {
    return std::min(std::max(int16_t(0),one_el),int16_t(255));
  }
} SSEReg;  // We will use higher precision work area, so we do not have to worry about overflow during compute

typedef struct {
  uint8_t el[NBYTES_PER_REG];

  // need to bound to emulate u8
  inline void load(const SSEReg val) {
    for (int j=0; j<NBYTES_PER_REG; j++) el[j] = SSEReg::enforce_u8(val.el[j]);
  }
} SSEMem;  // Memory representation is indeed uint8_t

inline SSEReg sse_set1_epi8(const int16_t a) {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = a;
  return out;
};

// not checking for overflow at this point
inline SSEReg sse_subs_epu8(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j] - y.el[j];
  return out;
};

inline SSEReg sse_load_u8(SSEMem const *x) {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x->el[j];
  return out;
};

inline void sse_store_u8(SSEMem *x, const SSEReg y) {
  x->load(y);
};

inline SSEReg sse_setzero_siall() {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = 0;
  return out;
};

inline SSEReg sse_max_epu8(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = std::max(x.el[j],y.el[j]);
  return out;
};

// Note: We are assuming x and y are already properly bound
inline SSEReg sse_or_siall(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j] | y.el[j];
  return out;
};

#define sse_anygt_epu8(val1,val2,outval) { \
        outval = (val1.el[0]>val2.el[0]) | (val1.el[1]>val2.el[1]) |(val1.el[2]>val2.el[2]) |(val1.el[3]>val2.el[3]); }

inline SSEReg sse_slli_u8(const SSEReg x) {
  SSEReg out;
  out.el[0] = 0;
  for (int j=1; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j-1];
  return out;
};

#define sse_set_low_u8(v, outval) outval.el[0] = v

#ifdef ENABLE_I16

//
// TODO: Right now ENABLE_I16 is not supported with SW simulation
//

inline SSEReg sse_adds_epi16(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<8; j++) {
     int32_t x1 = x.i16.el[j];
     int32_t y1 = y.i16.el[j];
     int32_t tmp32 = x1 + y1;
     int16_t tmp = std::min(std::max(-32768,tmp32),32767);
     out.i16.el[j] = tmp;
  }
  return out;
};

inline SSEReg sse_subs_epi16(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<8; j++) {
     int32_t x1 = x.i16.el[j];
     int32_t y1 = y.i16.el[j];
     int32_t tmp32 = x1 - y1;
     int16_t tmp = std::min(std::max(-32768,tmp32),32767);
     out.i16.el[j] = tmp;
  }
  return out;
};

inline SSEReg sse_set1_epi16(const int16_t a) {
  SSEReg out;
  for (int j=0; j<8; j++) out.i16.el[j] = a;
  return out;
};

inline SSEReg sse_max_epi16(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<8; j++) out.i16.el[j] = std::max(x.i16.el[j],y.i16.el[j]);
  return out;
};

inline SSEReg sse_slli_epi16(const SSEReg x, const unsigned int i) {
  SSEReg out;
  for (int j=0; j<8; j++) out.u16.el[j] = x.u16.el[j] << i;
  return out;
};



// TODO: Only used in sse_setall_ff
inline SSEReg sse_cmpeq_epi16(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<8; j++) out.u16.el[j] = (x.u16.el[j] == y.u16.el[j]) ? 0xFFFF : 0;
  return out;
};




inline SSEReg sse_cmpgt_epi16(const SSEReg x, const SSEReg y) {
  SSEReg out;
  for (int j=0; j<8; j++) out.u16.el[j] = (x.i16.el[j] > y.i16.el[j]) ? 0xFFFF : 0;
  return out;
};

inline SSEReg sse_srli_epi16(const SSEReg x, const unsigned int i) {
  SSEReg out;
  for (int j=0; j<8; j++) out.u16.el[j] = x.u16.el[j] >> i;
  return out;
};

inline int32_t sse_extract_epi16(const SSEReg x, const int32_t i) {
  return x.u16.el[i];
};

inline SSEReg nosse_set_low_i16(const uint16_t v) {
  SSEReg out;
  out.u16.el[0] = v;
  for (int j=1; j<8; j++) out.u16.el[j] = 0;
  return out;
};
#define sse_set_low_i16(inval, outval) { outval = nosse_set_low_i16(inval);}

#define sse_cmplt_epi16(x, y) sse_cmpgt_epi16(y, x)

#endif /* ENABLE_I16 */ 

#else /* no SSE_SW4 */

#define NBYTES_PER_REG 16
#define BYTES_LOG2_PER_REG 4
#define SSE_MASK_ALL 0xffff

#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif


#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
typedef simde__m128i SSEReg;
typedef SSEReg  SSEMem;  // memory and register representation are the same

#define sse_load_siall(x) simde_mm_load_si128(x)
#define sse_store_siall(x, y) simde_mm_store_si128(x, y)
#define sse_set1_epi8(x) simde_mm_set1_epi8(x)
#define sse_subs_epu8(x, y) simde_mm_subs_epu8(x, y)
#define sse_setzero_siall() simde_mm_setzero_si128()
#define sse_max_epu8(x, y) simde_mm_max_epu8(x, y)
#define sse_or_siall(x, y) simde_mm_or_si128(x, y)

#define sse_insert_epi16(x, y, z) simde_mm_insert_epi16(x, y, z)
#define sse_cmpeq_epi8(x, y) simde_mm_cmpeq_epi8(x, y)
#define sse_movemask_epi8(x) simde_mm_movemask_epi8(x)
#define sse_slli_siall(x, y) simde_mm_slli_si128(x, y)


#ifdef ENABLE_I16
#define sse_adds_epi16(x, y) simde_mm_adds_epi16(x, y)
#define sse_subs_epi16(x, y) simde_mm_subs_epi16(x, y)
#define sse_set1_epi16(x) simde_mm_set1_epi16(x)
#define sse_max_epi16(x, y) simde_mm_max_epi16(x, y)
#define sse_slli_epi16(x, y) simde_mm_slli_epi16(x, y)

#define sse_cmpeq_epi16(x, y) simde_mm_cmpeq_epi16(x, y)
#define sse_cmpgt_epi16(x, y) simde_mm_cmpgt_epi16(x, y)
#endif

#else
typedef __m128i SSEReg;
typedef SSEReg  SSEMem;  // memory and register representation are the same

#define sse_load_siall(x) _mm_load_si128(x)
#define sse_store_siall(x, y) _mm_store_si128(x, y)
#define sse_set1_epi8(x) _mm_set1_epi8(x)
#define sse_subs_epu8(x, y) _mm_subs_epu8(x, y)
#define sse_setzero_siall() _mm_setzero_si128()
#define sse_max_epu8(x, y) _mm_max_epu8(x, y)
#define sse_or_siall(x, y) _mm_or_si128(x, y)

#define sse_insert_epi16(x, y, z) _mm_insert_epi16(x, y, z)
#define sse_cmpeq_epi8(x, y) _mm_cmpeq_epi8(x, y)
#define sse_movemask_epi8(x) _mm_movemask_epi8(x)
#define sse_slli_siall(x, y) _mm_slli_si128(x, y)


#ifdef ENABLE_I16
#define sse_adds_epi16(x, y) _mm_adds_epi16(x, y)
#define sse_subs_epi16(x, y) _mm_subs_epi16(x, y)
#define sse_set1_epi16(x) _mm_set1_epi16(x)
#define sse_max_epi16(x, y) _mm_max_epi16(x, y)
#define sse_slli_epi16(x, y) _mm_slli_epi16(x, y)

#define sse_cmpeq_epi16(x, y) _mm_cmpeq_epi16(x, y)
#define sse_cmpgt_epi16(x, y) _mm_cmpgt_epi16(x, y)
#endif

#endif

#endif /* not x86 */

#endif /* SSE_AVX2 */

#ifndef SSE_DISABLE

#define sse_load_u8(x) sse_load_siall(x)
#define sse_store_u8(x, y) sse_store_siall(x, y)

#ifdef ENABLE_I16
#define sse_load_i16(x) sse_load_siall(x)
#define sse_store_i16(x, y) sse_store_siall(x, y)
#endif

#endif

/* Fill all elements in outval with inval */
#define sse_fill_i8(inval, outval) outval=sse_set1_epi8(inval)

#ifndef SSE_DISABLE

#define sse_slli_u8(x) sse_slli_siall(x,1)

#define sse_anygt_epu8(val1,val2,outval) { \
	SSEReg s = sse_subs_epu8(val1, val2); \
        s = sse_cmpeq_epi8(s, vzero); \
        outval = (sse_movemask_epi8(s) != SSE_MASK_ALL); }

/* Set the low element with invl, all others to 0 */
#define sse_set_low_i16(inval, outval) { \
	outval = sse_setzero_siall(); \
	outval = sse_insert_epi16(outval, inval, 0); \
}

#define sse_set_low_u8(inval, outval) sse_set_low_i16(inval, outval)

#endif /* SSE_DISABLE */


#ifdef ENABLE_I16

#define sse_slli_i16(x) sse_slli_siall(x, 2)

/* Fill all elements in outval with inval */
/* opt version will check for special ivals that can use shortcuts */
#define sse_fill_i16(inval, outval) outval=sse_set1_epi16(inval)

#define sse_anygt_epi16(val1,val2,outval) { \
	SSEReg s = sse_cmpgt_epi16(val1, val2); \
        outval = (sse_movemask_epi8(s) != 0); }

#define sse_setall_ff(val) val =  sse_cmpeq_epi16(val,val)

#endif

#endif /* SSE_WRAP_H_ */
