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

typedef __m256i SSERegI;

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

#if defined(SSE_SW16) || defined(SSE_SW8) || defined(SSE_SW4) || defined(SSE_SW2)

#define SSE_DISABLE 1

#if defined(SSE_SW16)
#define NBYTES_PER_REG 16
#define BYTES_LOG2_PER_REG 4

#elif defined(SSE_SW8)
#define NBYTES_PER_REG 8
#define BYTES_LOG2_PER_REG 3

#elif defined(SSE_SW4)
#define NBYTES_PER_REG 4
#define BYTES_LOG2_PER_REG 2

#elif defined(SSE_SW2)
#define NBYTES_PER_REG 2
#define BYTES_LOG2_PER_REG 1
#endif

#include <algorithm>
#include <stdint.h>

typedef struct {
  uint8_t el[NBYTES_PER_REG];
} SSERegI;


inline SSERegI sse_set1_epi8(const uint8_t a) {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = a;
  return out;
};

inline SSERegI sse_subs_epu8(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j] - std::min(x.el[j],y.el[j]);
  return out;
};

inline SSERegI sse_load_siall(SSERegI const *x) {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x->el[j];
  return out;
};

inline void sse_store_siall(SSERegI *x, const SSERegI y) {
  for (int j=0; j<NBYTES_PER_REG; j++) x->el[j] = y.el[j];
};

inline SSERegI sse_setzero_siall() {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = 0;
  return out;
};

inline SSERegI sse_max_epu8(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = std::max(x.el[j],y.el[j]);
  return out;
};

inline SSERegI sse_or_siall(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j] | y.el[j];
  return out;
};

#if defined(SSE_SW16)
#define sse_anygt_epu8(val1,val2,outval) { \
        outval = (val1.el[0]>val2.el[0]) | (val1.el[1]>val2.el[1]) |(val1.el[2]>val2.el[2]) |(val1.el[3]>val2.el[3]) || \
                 (val1.el[4]>val2.el[4]) | (val1.el[5]>val2.el[5]) |(val1.el[6]>val2.el[6]) |(val1.el[7]>val2.el[7]) || \
                 (val1.el[8]>val2.el[8]) | (val1.el[9]>val2.el[9]) |(val1.el[10]>val2.el[10]) |(val1.el[11]>val2.el[11]) || \
                 (val1.el[12]>val2.el[12]) | (val1.el[13]>val2.el[13]) |(val1.el[14]>val2.el[14]) |(val1.el[15]>val2.el[15]); }

#elif defined(SSE_SW8)
#define sse_anygt_epu8(val1,val2,outval) { \
        outval = (val1.el[0]>val2.el[0]) | (val1.el[1]>val2.el[1]) |(val1.el[2]>val2.el[2]) |(val1.el[3]>val2.el[3]) || \
                 (val1.el[4]>val2.el[4]) | (val1.el[5]>val2.el[5]) |(val1.el[6]>val2.el[6]) |(val1.el[7]>val2.el[7]); }

#elif defined(SSE_SW4)
#define sse_anygt_epu8(val1,val2,outval) { \
        outval = (val1.el[0]>val2.el[0]) | (val1.el[1]>val2.el[1]) |(val1.el[2]>val2.el[2]) |(val1.el[3]>val2.el[3]); }

#elif defined(SSE_SW2)
#define sse_anygt_epu8(val1,val2,outval) { \
        outval = (val1.el[0]>val2.el[0]) | (val1.el[1]>val2.el[1]); }

#endif

inline SSERegI sse_slli_u8(const SSERegI x) {
  SSERegI out;
  out.el[0] = 0;
  for (int j=1; j<NBYTES_PER_REG; j++) out.el[j] = x.el[j-1];
  return out;
};

#define sse_set_low_u8(v, outval) outval.el[0] = v

#ifdef ENABLE_I16

//
// TODO: Right now ENABLE_I16 is not supported with SW simulation
//

#if 0
inline SSERegI sse_adds_epi16(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) {
     int32_t x1 = x.i16.el[j];
     int32_t y1 = y.i16.el[j];
     int32_t tmp32 = x1 + y1;
     int16_t tmp = std::min(std::max(-32768,tmp32),32767);
     out.i16.el[j] = tmp;
  }
  return out;
};

inline SSERegI sse_subs_epi16(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) {
     int32_t x1 = x.i16.el[j];
     int32_t y1 = y.i16.el[j];
     int32_t tmp32 = x1 - y1;
     int16_t tmp = std::min(std::max(-32768,tmp32),32767);
     out.i16.el[j] = tmp;
  }
  return out;
};

inline SSERegI sse_set1_epi16(const int16_t a) {
  SSERegI out;
  for (int j=0; j<8; j++) out.i16.el[j] = a;
  return out;
};

inline SSERegI sse_max_epi16(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) out.i16.el[j] = std::max(x.i16.el[j],y.i16.el[j]);
  return out;
};

inline SSERegI sse_slli_epi16(const SSERegI x, const unsigned int i) {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = x.u16.el[j] << i;
  return out;
};



// TODO: Only used in sse_setall_ff
inline SSERegI sse_cmpeq_epi16(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = (x.u16.el[j] == y.u16.el[j]) ? 0xFFFF : 0;
  return out;
};




inline SSERegI sse_cmpgt_epi16(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = (x.i16.el[j] > y.i16.el[j]) ? 0xFFFF : 0;
  return out;
};

inline SSERegI sse_srli_epi16(const SSERegI x, const unsigned int i) {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = x.u16.el[j] >> i;
  return out;
};

inline int32_t sse_extract_epi16(const SSERegI x, const int32_t i) {
  return x.u16.el[i];
};

inline SSERegI nosse_set_low_i16(const uint16_t v) {
  SSERegI out;
  out.u16.el[0] = v;
  for (int j=1; j<8; j++) out.u16.el[j] = 0;
  return out;
};
#define sse_set_low_i16(inval, outval) { outval = nosse_set_low_i16(inval);}

#define sse_cmplt_epi16(x, y) sse_cmpgt_epi16(y, x)

#endif

#endif /* ENABLE_I16 */ 

#elif defined(SSE_SCALAR)

#define SSE_DISABLE 1

#define NBYTES_PER_REG 1
#define BYTES_LOG2_PER_REG 0

#include <algorithm>
#include <stdint.h>

typedef uint8_t SSERegI;


inline SSERegI sse_set1_epi8(const uint8_t a) { return a; }

inline SSERegI sse_subs_epu8(const SSERegI x, const SSERegI y) { return x - std::min(x,y); }

inline SSERegI sse_load_siall(SSERegI const *x) { return *x; }

inline void sse_store_siall(SSERegI *x, const SSERegI y) { *x = y; }

inline SSERegI sse_setzero_siall() { return 0; }

inline SSERegI sse_max_epu8(const SSERegI x, const SSERegI y) { return std::max(x,y); }

inline SSERegI sse_or_siall(const SSERegI x, const SSERegI y) { return x | y; }

#define sse_set_low_u8(v, outval) outval = v

#define sse_slli_u8(x) 0

#ifdef ENABLE_I16

//
// TODO: Right now ENABLE_I16 is not supported with SW simulation
//

#endif /* ENABLE_I16 */ 

#else /* no SSE_SCALAR */

#define NBYTES_PER_REG 16
#define BYTES_LOG2_PER_REG 4
#define SSE_MASK_ALL 0xffff

#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif


#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
typedef simde__m128i SSERegI;

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
typedef __m128i SSERegI;

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

/* Fill all elements in outval with inval */
#define sse_fill_i8(inval, outval) outval=sse_set1_epi8(inval)

#ifndef SSE_DISABLE

#define sse_slli_u8(x) sse_slli_siall(x,1)

#define sse_anygt_epu8(val1,val2,outval) { \
	SSERegI s = sse_subs_epu8(val1, val2); \
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
	SSERegI s = sse_cmpgt_epi16(val1, val2); \
        outval = (sse_movemask_epi8(s) != 0); }

#define sse_setall_ff(val) val =  sse_cmpeq_epi16(val,val)

#endif

#endif /* SSE_WRAP_H_ */
