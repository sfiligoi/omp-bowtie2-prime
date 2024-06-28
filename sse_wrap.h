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


#define sse_adds_epi16(x, y) _mm256_adds_epi16(x, y)
#define sse_subs_epi16(x, y) _mm256_subs_epi16(x, y)
#define sse_set1_epi16(x) _mm256_set1_epi16(x)
#define sse_max_epi16(x, y) _mm256_max_epi16(x, y)

#define sse_cmpeq_epi16(x, y) _mm256_cmpeq_epi16(x, y)
#define sse_slli_epi16(x, y) _mm256_slli_epi16(x, y)

#else /* no SSE_AVX2 */

#define NBYTES_PER_REG 16
#define BYTES_LOG2_PER_REG 4
#define SSE_MASK_ALL 0xffff

#ifdef SSE_DISABLE

#include <algorithm>
#include <stdint.h>

typedef struct {
  int16_t el[8];
} SSERegI16;

typedef struct {
  uint16_t el[8];
} SSERegU16;

typedef struct {
  int8_t el[16];
} SSERegI8;

typedef struct {
  uint8_t el[16];
} SSERegU8;

typedef union {
  SSERegI16 i16;
  SSERegU16 u16;
  SSERegI8  i8;
  SSERegU8  u8;
} SSERegI;

inline SSERegI sse_set1_epi8(const int8_t a) {
  SSERegI out;
  for (int j=0; j<16; j++) out.i8.el[j] = a;
  return out;
};

inline SSERegI sse_subs_epu8(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<16; j++) {
     int16_t x1 = x.u8.el[j];
     int16_t y1 = y.u8.el[j];
     int16_t tmp16 = x1 - y1;
     uint8_t tmp = std::min(std::max(int16_t(0),tmp16),int16_t(255));
     out.u8.el[j] = tmp;
  }
  return out;
};

inline SSERegI sse_load_siall(SSERegI const *x) {
  return x[0];
};

inline void sse_store_siall(SSERegI *x, const SSERegI y) {
  x[0] = y;
};

inline SSERegI sse_setzero_siall() {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = 0;
  return out;
};

inline SSERegI sse_max_epu8(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<16; j++) out.u8.el[j] = std::max(x.u8.el[j],y.u8.el[j]);
  return out;
};

inline SSERegI sse_or_siall(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<8; j++) out.u16.el[j] = x.u16.el[j] | y.u16.el[j];
  return out;
};


// TODO: We really just need sse_anygt_epu8
inline SSERegI sse_cmpeq_epi8(const SSERegI x, const SSERegI y) {
  SSERegI out;
  for (int j=0; j<16; j++) out.u8.el[j] = (x.u8.el[j] == y.u8.el[j]) ? 0xFF : 0;
  return out;
};

inline uint16_t sse_movemask_epi8(const SSERegI x) {
  uint16_t out = 0;
  for (int j=0; j<16; j++) out |= ((uint16_t)(x.u8.el[j]>>7)) << j;
  return out;
};

inline SSERegI sse_slli_siall(const SSERegI x, const unsigned int i) {
  SSERegI out;
  if (i==2) {
     out.u16.el[0] = 0;
     for (int j=1; j<8; j++) out.u16.el[j] = x.u16.el[j-1];
  } else if (i==1) {
     out.u8.el[0] = 0;
     for (int j=1; j<16; j++) out.u8.el[j] = x.u8.el[j-1];
  } else {
     throw 1; // unsupported
  }
  return out;
};


// Note: We are not using saturation, not needed

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

#else /* no SSE_DISABLE */

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


#define sse_adds_epi16(x, y) simde_mm_adds_epi16(x, y)
#define sse_subs_epi16(x, y) simde_mm_subs_epi16(x, y)
#define sse_set1_epi16(x) simde_mm_set1_epi16(x)
#define sse_max_epi16(x, y) simde_mm_max_epi16(x, y)

#define sse_cmpeq_epi16(x, y) simde_mm_cmpeq_epi16(x, y)
#define sse_slli_epi16(x, y) simde_mm_slli_epi16(x, y)

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


#define sse_adds_epi16(x, y) _mm_adds_epi16(x, y)
#define sse_subs_epi16(x, y) _mm_subs_epi16(x, y)
#define sse_set1_epi16(x) _mm_set1_epi16(x)
#define sse_max_epi16(x, y) _mm_max_epi16(x, y)

#define sse_cmpeq_epi16(x, y) _mm_cmpeq_epi16(x, y)
#define sse_slli_epi16(x, y) _mm_slli_epi16(x, y)

#endif

#endif /* SSE_DISABLE */

#endif /* SSE_AVX2 */

/* Fill all elements in outval with inval */
#define sse_fill_i8(inval, outval) outval=sse_set1_epi8(inval)

#define sse_anygt_epu8(val1,val2,outval) { \
	SSERegI s = sse_subs_epu8(val1, val2); \
        s = sse_cmpeq_epi8(s, vzero); \
        outval = (sse_movemask_epi8(s) != SSE_MASK_ALL); }

#define sse_slli_u8(x) sse_slli_siall(x,1)

#ifndef SSE_DISABLE

/* Set the low element with invl, all others to 0 */
#define sse_set_low_i16(inval, outval) { \
	outval = sse_setzero_siall(); \
	outval = sse_insert_epi16(outval, inval, 0); \
}

#endif /* SSE_DISABLE */

#define sse_set_low_u8(inval, outval) sse_set_low_i16(inval, outval)



#define sse_slli_i16(x) sse_slli_epi16(vf,1)

/* Fill all elements in outval with inval */
/* opt version will check for special ivals that can use shortcuts */
#define sse_fill_i16(inval, outval) outval=sse_set1_epi16(inval)

#define sse_anygt_epi16(val1,val2,outval) { \
	SSERegI s = sse_cmpgt_epi16(val1, val2); \
        outval = (sse_movemask_epi8(s) != 0); }

#define sse_setall_ff(val) val =  sse_cmpeq_epi16(val,val)

#endif /* SSE_WRAP_H_ */
