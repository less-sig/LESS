/**
 *
 * Optimized Implementation of LESS.
 *
 * @version 1.2 (May 2025)
 *
 * @author Floyd Zweydinger
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

#if !defined(__APPLE__) && !defined(__aarch64__) && !defined(_M_ARM64)
#error "only available on aarch64"
#endif

#include <_types/_uint32_t.h>
#include <_types/_uint64_t.h>
#include <arm_neon.h>

#ifndef __clang__
#include <arm_bf16.h>
#include <arm_fp16.h>
#include <stdint.h>
typedef __Uint8x8_t  __uint8x8_t;
typedef __Uint8x16_t __uint8x16_t;
typedef __Uint16x4_t __uint16x4_t;
typedef __Uint16x8_t __uint16x8_t;
typedef __Uint32x2_t __uint32x2_t;
typedef __Uint32x4_t __uint32x4_t;
typedef __Uint64x1_t __uint64x1_t;
typedef __Uint64x2_t __uint64x2_t;
#else
typedef __attribute__((neon_vector_type(8))) uint8_t  __uint8x8_t;
typedef __attribute__((neon_vector_type(16))) uint8_t __uint8x16_t;
typedef __attribute__((neon_vector_type(4))) uint16_t __uint16x4_t;
typedef __attribute__((neon_vector_type(8))) uint16_t __uint16x8_t;
typedef __attribute__((neon_vector_type(2))) uint32_t __uint32x2_t;
typedef __attribute__((neon_vector_type(4))) uint32_t __uint32x4_t;
typedef __attribute__((neon_vector_type(1))) uint64_t __uint64x1_t;
typedef __attribute__((neon_vector_type(2))) uint64_t __uint64x2_t;
#endif

typedef union {
	uint8_t  v8 [32];
	uint16_t v16[16];
	uint32_t v32[8];
	uint64_t v64[4];
	__uint8x16_t v[2];
} vec256_t;

typedef union {
	uint8_t  v8 [16];
	uint16_t v16[8];
	uint32_t v32[4];
	uint64_t v64[2];
	__uint8x16_t v;
} vec128_t;


/// NOTE: ptr must be (v256_t *)
// c <- src
#define vload256(c, src) c.v[0] = vld1q_u8((uint8_t *)src); c.v[1] = vld1q_u8(((uint8_t *)src) + 16u);
#define vload128(c, src) c.v[0] = vld1q_u8((uint8_t *)src);

// src <- c
#define vstore256(src, c) vst1q_u8((uint8_t *)src, c.v[0]); vst1q_u8(((uint8_t *)src) + 16, c.v[1]);
#define vstore128(src, c) vst1q_u8((uint8_t *)src, c.v);

// c = a + b
#define vadd8(c, a, b)   c.v[0] = vaddq_u8(a.v[0], b.v[0]); c.v[1] = vaddq_u8(a.v[1], b.v[1]);
#define vadd16(c, a, b)  c.v[0] = (__uint8x16_t)vaddq_u16((__uint16x8_t)a.v[0], (__uint16x8_t)b.v[0]); c.v[1] = (__uint8x16_t)vaddq_u16((__uint16x8_t)a.v[1], (__uint16x8_t)b.v[1]);
#define vadd64(c, a, b)  c.v[0] = vaddq_u64(a.v[0], b.v[0]); c.v[1] = vaddq_u64(a.v[1], b.v[1]);

// c = a - b
#define vsub8(c, a, b)   c.v[0] = vsubq_u8(a.v[0], b.v[0]); c.v[1] = vsubq_u8(a.v[1], b.v[1]);
#define vsub16(c, a, b)  c.v[0] = (__uint8x16_t)vsubq_u16((__uint16x8_t)a.v[0], (__uint16x8_t)b.v[0]); c.v[1] = (__uint8x16_t)vsubq_u16((__uint16x8_t)a.v[1], (__uint16x8_t)b.v[1]);

// c = a * b
//#define vmul_hi16(c, a, b) c = _mm256_mulhi_epu16(a, b);
#define vmul_lo16(c, a, b)  c.v[0] = (__uint8x16_t)vmulq_u16((__uint16x8_t)a.v[0], (__uint16x8_t)b.v[0]); c.v[1] = (__uint8x16_t)vmulq_u16((__uint16x8_t)a.v[1], (__uint16x8_t)b.v[1]);

// c = a >> n
#define vsr16(c, a, b)  c.v[0] = (__uint8x16_t)vshrq_n_u16((__uint16x8_t)a.v[0], b); c.v[1] = (__uint8x16_t)vshrq_n_u16((__uint16x8_t)a.v[1], b);

// c = a << n
#define vsl16(c, a, b)  c.v[0] = (__uint8x16_t)vshlq_n_u16((__uint16x8_t)a.v[0], b); c.v[1] = (__uint8x16_t)vshlq_n_u16((__uint16x8_t)a.v[1], b);

// c = a & b
#define vand(c, a, b)   c.v[0] = vandq_u8(a.v[0], b.v[0]); c.v[1] = vandq_u8(a.v[1], b.v[1]);

// c = a ^ u
#define vxor(c, a, b)   c.v[0] = veorq_u8(a.v[0], b.v[0]); c.v[1] = veorq_u8(a.v[1], b.v[1]);

// c = a | u
#define vor(c, a, b)    c.v[0] = vorrq_u8(a.v[0], b.v[0]); c.v[1] = vorrq_u8(a.v[1], b.v[1]);

// c[0..16] = n
#define vset8(c, n)   c.v[0] = vdupq_n_u8((uint8_t)n); c.v[1] = vdupq_n_u8((uint8_t)n);
#define vset16(c, n)  c.v[0] = (__uint16x8_t)vdupq_n_u16(n); c.v[1] = (__uint16x8_t)vdupq_n_u16(n);

// c = a == b
#define vcmp8(c, a, b) c.v[0] = vceqq_u8(a.v[0], b.v[0]); c.v[1] = vceqq_u8(a.v[1], b.v[1]);

// moves the msb of each uint8_t limb int a single bit
static inline uint32_t vmovemask8(const vec256_t a) {
	uint16x8_t high_bits0 = vreinterpretq_u16_u8(vshrq_n_u8(a.v[0], 7));
	uint16x8_t high_bits1 = vreinterpretq_u16_u8(vshrq_n_u8(a.v[1], 7));
	uint32x4_t paired160  = vreinterpretq_u32_u16(vsraq_n_u16(high_bits0, high_bits0, 7));
	uint32x4_t paired161  = vreinterpretq_u32_u16(vsraq_n_u16(high_bits1, high_bits1, 7));
	uint64x2_t paired320  = vreinterpretq_u64_u32(vsraq_n_u32(paired160, paired160, 14));
	uint64x2_t paired321  = vreinterpretq_u64_u32(vsraq_n_u32(paired161, paired161, 14));
	uint8x16_t paired640  = vreinterpretq_u8_u64(vsraq_n_u64(paired320, paired320, 28));
	uint8x16_t paired641  = vreinterpretq_u8_u64(vsraq_n_u64(paired321, paired321, 28));
	return (uint32_t)((vgetq_lane_u8(paired640, 0) <<  0) | ((int) vgetq_lane_u8(paired640, 8) <<  8) |
		              (vgetq_lane_u8(paired641, 0) << 16) | ((int) vgetq_lane_u8(paired641, 8) << 24));
}

// Unpack 8-bit low: a[0] | b[0] ... a[7] | b[7]
#define vunpackl8(c, a, b) c.v[0] = vzip1q_u8(a.v[0], b.v[0]); c.v[1] = vzip1q_u8(a.v[1], b.v[1]);
// Unpack 8-bit hi: a[8] | b[8] ... a[15] | b[15]
#define vunpackh8(c, a, b) c.v[0] = vzip2q_u8(a.v[0], b.v[0]); c.v[1] = vzip2q_u8(a.v[1], b.v[1]);

// Extend 8-bit unsigned from vec128 to
// 16-bit unsigned from vec256
#define vextend8_16(c, a) c.v[0] = vmovl_u8(vget_low_u8(a.v)); c.v[1] = vmovl_u8(vget_high_u8(a.v));

// Select low 128-bit of vec256
#define vget_lo(c, a) c.v = a.v[0];
// Select high 128-bit of vec256
#define vget_hi(c, a) c.v = a.v[1];

// Set vector to 0
#define vzero(c) c.v[0] = 0; v.v[1] = 0;

// Shuffle vector a according to vector b
#define vshuffle8(c, a, b) c.v[0] = vqtbl1q_s8(a.v[0], b.v[0]); c.v[1] = vqtbl1q_s8(a.v[1], b.v[1]);

// Permute 128-bit, combine 2 128-bit low from
// a and b to c
// c = a[0..127] | b[0..127]
#define vpermute(c, a, b) c.v[0] = a.v[0]; c.v[1] = b.v[0];
#define vpermute2(c, a, b, d) 							\
	if (d & (1u << 3u)) {								\
		c.v[0] = vdupq_n_u8(0);							\
	} else { 											\
		c.v[0] = (d&3u) <= 1u ? a.v[d&1u] : b.v[(d>>2u)&1u];		\
	}													\
	if (d & (1u << 7u)) {								\
		c.v[1] = vdupq_n_u8(0);							\
	} else { 											\
		c.v[1] = ((d>>4)&3u) <= 1 ? a.v[d>>4u] : b.v[d>>6u];\
	}

#define vpermute_4x64(c, a, b) 												\
	{																		\
		vec256_t vpermute_4x64_tmp;											\
		for (uint32_t ijk = 0; ijk < 4; ijk++) {							\
			vpermute_4x64_tmp.v64[ijk] = a.v64[(b >> (2*ijk)) & 0x3];		\
		}																	\
		c.v[0] = vpermute_4x64_tmp.v[0];									\
		c.v[1] = vpermute_4x64_tmp.v[1];									\
	}

/*
 * Fix width 8-bit Barrett modulo reduction Q = 127
 * c = a % q
 */
#define barrett_red8(a, t, c127, c1) 	\
	vsr16(t, a, 7)						\
	vand(t, t, c1)						\
	vadd8(a, a, t)						\
	vand(a, a, c127)



/*
 * Fix width 16-bit Barrett multiplication Q = 127
 * c = (a * b) % q
 */
#define barrett_mul_u16(c, a, b, t)        \
    vmul_lo16(a, a, b); /* lo = (a * b)  */  \
    vsr16(t, a, 7);     /* hi = (lo >> 7) */ \
    vadd16(a, a, t);    /* lo = (lo + hi) */ \
    vsl16(t, t, 7);     /* hi = (hi << 7) */ \
    vsub16(c, a, t);    /* c  = (lo - hi) */

/// number of 8 bit elements in an avx register
#define LESS_WSZ 32

/// number of avx registers needed for a full row in a generator matrix
#define NW ((N + LESS_WSZ - 1) / LESS_WSZ)

/// original reduction formula
#define W_RED127(x)     \
	vsr16(t, x, 7)		\
	vand(t, t, c01)		\
	vadd8(x, x, t)		\
	vand(x, x, c7f)


/// new reduction formula, catches the case where input is q=127
#define W_RED127_(x)     \
	x.v[0] = vandq_u8(vaddq_u8(vandq_u8(vshrq_n_u16(vaddq_u8(x.v[0], c01.v[0]), 7), c01.v[0]), x.v[0]), c7f.v[0]); \
	x.v[1] = vandq_u8(vaddq_u8(vandq_u8(vshrq_n_u16(vaddq_u8(x.v[1], c01.v[1]), 7), c01.v[1]), x.v[1]), c7f.v[1]);

// Extend from 8-bit to 16-bit type
extern const uint8_t shuff_low_half[32];

void print256_num(vec256_t var, const char *string);

static inline uint8_t vhadd8(const vec256_t t) {
	vec128_t v;
	v.v = vaddq_u8(t.v[0], t.v[1]);
	uint8_t r = fq_red(v.v8[0]);
	for (uint32_t i = 1; i < 16; i++) {
		r = fq_add(r, fq_red(v.v8[i]));
	}

	return r;
}
