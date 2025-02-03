/**
 *
 * Optimized Implementation of LESS.
 *
 * @version 1.2 (May 2025)
 * @author Floyd Zweydinger <zweydfg8+github@rub.de>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#pragma once

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <immintrin.h>

#include "fq_arith.h"

typedef __m256i vec256_t;
typedef __m128i vec128_t;


/// number of Fq elements per vector register
#define LESS_WSZ 32u

// number of vector register for N bytes
#define NW ((NEXT_MULTIPLE(N, LESS_WSZ))/LESS_WSZ)


// c <- src
#define vload256(c, src) c = _mm256_loadu_si256(src);
#define vload128(c, src) c = _mm_loadu_si128(src);

// src <- c
#define vstore256(src, c) _mm256_storeu_si256(src, c);
#define vstore128(src, c) _mm_storeu_si128(src, c);
// #define vstore(src, c) _mm256_store_si256(src, c);

// c = a + b
#define vadd8(c, a, b)  c = _mm256_add_epi8(a, b);
#define vadd16(c, a, b) c = _mm256_add_epi16(a, b);
#define vadd64(c, a, b) c = _mm256_add_epi64(a, b);

// c = a - b
#define vsub8(c, a, b) c = _mm256_sub_epi8(a, b);
#define vsub16(c, a, b) c = _mm256_sub_epi16(a, b);

// c = a * b
#define vmul_lo16(c, a, b) c = _mm256_mullo_epi16(a, b);
#define vmul_hi16(c, a, b) c = _mm256_mulhi_epu16(a, b);

// c = a >> n
#define vsr16(c, a, n) c = _mm256_srai_epi16(a, n);
#define vsr32(c, a, n) c = _mm256_srli_epi32(a, n);

// c = a << n
#define vsl16(c, a, n) c = _mm256_slli_epi16(a, n);

// c = a & b
#define vand(c, a, b) c = _mm256_and_si256(a, b);

// c = a ^ b
#define vxor(c, a, b) c = _mm256_xor_si256(a, b);

// c = a | b
#define vor(c, a, b) c = _mm256_or_si256(a, b);

// c[0..16] = n
#define vset8(c, n) c = _mm256_set1_epi8((char)n);
#define vset17(c, n) c = _mm256_set1_epi16((short)n);

// c = a == b
#define vcmp8(c, a, b) c = _mm256_cmpeq_epi8(a, b);

/// move the msb of each 8bit limb into an integer
#define vmovemask8(a) _mm256_movemask_epi8(a);

// Unpack 8-bit low: a[0] | b[0] ... a[7] | b[7]
#define vunpackl8(c, a, b) c = _mm256_unpacklo_epi8(a, b);
// Unpack 8-bit hi: a[8] | b[8] ... a[15] | b[15]
#define vunpackh8(c, a, b) c = _mm256_unpackhi_epi8(a, b);

// Extend 8-bit unsigned from vec128 to
// 16-bit unsigned from vec256
#define vextend8_16(c, a) c = _mm256_cvtepu8_epi16(a);

// Select low 128-bit of vec256
#define vget_lo(c, a) c = _mm256_extracti128_si256(a, 0);
// Select high 128-bit of vec256
#define vget_hi(c, a) c = _mm256_extracti128_si256(a, 1);

// Set vector to 0
#define vzero(c) c = _mm256_setzero_si256();

// Shuffle vector a according to vector b
#define vshuffle8(c, a, b) c = _mm256_shuffle_epi8(a, b);

// Permute 128-bit, combine 2 128-bit low from
// a and b to c
// c = a[0..127] | b[0..127]
#define vpermute(c, a, b) c = _mm256_permute2x128_si256(a, b, 0x20);
#define vpermute2(c, a, b, d) c = _mm256_permute2x128_si256(a, b, d);

#define vpermute_4x64(c, a, b) c = _mm256_permute4x64_epi64(a, b);

/*
 * Fix width 16-bit Barrett modulo reduction Q = 127
 * c = a % q
 */
#define barrett_red16(c, a, t, c127, c516, c1)                 \
    t = _mm256_add_epi16(a, c1);     /* t = (a + 1) */         \
    t = _mm256_mulhi_epu16(t, c516); /* t = (a * 516) >> 16 */ \
    t = _mm256_mullo_epi16(t, c127); /* t = (t * Q) */         \
    a = _mm256_sub_epi16(a, t);      /* a = (a - t)*/

/*
 * Fix width 8-bit Barrett modulo reduction Q = 127
 * c = a % q
 */
#define barrett_red8(a, t, c127, c1) \
    t = _mm256_srli_epi16(a, 7);     \
    t = _mm256_and_si256(t, c1);     \
    a = _mm256_add_epi8(a, t);       \
    a = _mm256_and_si256(a, c127);

#define barrett_red8_half(a, t, c127, c1) \
    t = _mm_srli_epi16(a, 7);             \
    t = _mm_and_si128(t, c1);             \
    a = _mm_add_epi8(a, t);               \
    a = _mm_and_si128(a, c127);

/*
 * t = tmp register
 * Fix width 16-bit Barrett multiplication Q = 127
 * c = (a * b) % q
 */
#define barrett_mul_u16(c, a, b, t)          \
    vmul_lo16(a, a, b); /* lo = (a * b)  */  \
    vsr16(t, a, 7);     /* hi = (lo >> 7) */ \
    vadd16(a, a, t);    /* lo = (lo + hi) */ \
    vsl16(t, t, 7);     /* hi = (hi << 7) */ \
    vsub16(c, a, t);    /* c  = (lo - hi) */

/// original reduction formula
#define W_RED127(x)                                                            \
  {                                                                            \
    t = _mm256_srli_epi16(x, 7);                                               \
    t = _mm256_and_si256(t, c01);                                              \
    x = _mm256_add_epi8(x, t);                                                 \
    x = _mm256_and_si256(x, c7f);                                              \
  }

/// new reduction formula, catches the case where input is q=127
#define W_RED127_(x)                                                           \
  x = _mm256_and_si256(                                                        \
      _mm256_add_epi8(_mm256_and_si256(                                        \
                          _mm256_srli_epi16(_mm256_add_epi8(x, c01), 7), c01), \
                      x),                                                      \
      c7f);

// Extend from 8-bit to 16-bit type
extern const uint8_t shuff_low_half[32];

void print256_num(vec256_t var, const char *string);

/// \return in[0] + in[1] + ... + in[31] % q
static inline uint8_t vhadd8(const __m256i in) {
    vec256_t c01, c7f;
    vset8(c01, 0x01);
    vset8(c7f, 0x7F);

    __m256i a = _mm256_srli_epi16(in, 8);
    __m256i t = _mm256_add_epi8(a, in);
    W_RED127_(t)

    a = _mm256_srli_epi32(t, 16);
    t = _mm256_add_epi8(a, t);
    W_RED127_(t)

    a = _mm256_srli_epi64(t, 32);
    t = _mm256_add_epi8(a, t);
    W_RED127_(t)

    a = _mm256_srli_si256(t, 8);
    t = _mm256_add_epi8(a, t);
    W_RED127_(t)

    a = _mm256_permute2x128_si256(t, t, 1);
    t = _mm256_add_epi8(a, t);
    W_RED127_(t)

    return _mm256_extract_epi8(t, 0);
}
