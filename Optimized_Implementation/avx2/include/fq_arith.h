/**
 *
 * Optimized Implementation of LESS.
 *
 * @version 1.2 (May 2023)
 *
 * @author Duc Tri Nguyen <dnguye69@gmu.edu>
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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

#include "parameters.h"
#include "lookup_table.h"
#include "rng.h"

#define NUM_BITS_Q (BITS_TO_REPRESENT(Q))


#define DEF_RAND_STATE(FUNC_NAME, EL_T, MINV, MAXV) \
static inline void FUNC_NAME(SHAKE_STATE_STRUCT *shake_monomial_state, EL_T *buffer, size_t num_elements) { \
   typedef uint64_t WORD_T; \
   static const EL_T MIN_VALUE = (MINV);\
   static const EL_T MAX_VALUE = (MAXV); \
   static const EL_T SPAN = MAX_VALUE - MIN_VALUE; \
   static const size_t REQ_BITS = BITS_TO_REPRESENT(SPAN); \
   static const EL_T EL_MASK = ((EL_T) 1 << REQ_BITS) - 1; \
   WORD_T word; \
   size_t count = 0; \
   do { \
      csprng_randombytes((unsigned char *) &word, sizeof(WORD_T), shake_monomial_state); \
      for (uint32_t i = 0; i < ((sizeof(WORD_T)*8u) / REQ_BITS); i++) { \
         EL_T rnd_value = word & EL_MASK; \
         if (rnd_value <= SPAN) buffer[count++] = rnd_value + MIN_VALUE; \
         if (count >= num_elements) return; \
         word >>= REQ_BITS; \
      } \
   } while (1); }

#define DEF_RAND(FUNC_NAME, EL_T, MINV, MAXV) \
static inline void FUNC_NAME(EL_T *buffer, size_t num_elements) { \
   typedef uint64_t WORD_T; \
   static const EL_T MIN_VALUE = (MINV); \
   static const EL_T MAX_VALUE = (MAXV); \
   static const EL_T SPAN = MAX_VALUE - MIN_VALUE; \
   static const size_t REQ_BITS = BITS_TO_REPRESENT(SPAN); \
   static const EL_T EL_MASK = ((EL_T) 1 << REQ_BITS) - 1; \
   WORD_T word; \
   size_t count = 0; \
   do { \
      randombytes((unsigned char *) &word, sizeof(WORD_T)); \
      for (uint32_t i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
         EL_T rnd_value = word & EL_MASK; \
         if (rnd_value <= SPAN) buffer[count++] = rnd_value + MIN_VALUE; \
         if (count >= num_elements) return; \
         word >>= REQ_BITS; \
      } \
   } while (1); }


/* GCC actually inlines and vectorizes Barrett's reduction already.
 * Backup implementation for less aggressive compilers follows */
static inline
FQ_ELEM fq_cond_sub(const FQ_ELEM x) {
    // equivalent to: (x >= Q) ? (x - Q) : x
    // likely to be ~ constant-time (a "smart" compiler might turn this into conditionals though)
    FQ_ELEM sub_q = x - Q;
    FQ_ELEM mask = -(sub_q >> NUM_BITS_Q);
    return (mask & Q) + sub_q;
}

static inline
FQ_ELEM fq_red(const FQ_DOUBLEPREC x) {
    return fq_cond_sub((x >> NUM_BITS_Q) + ((FQ_ELEM) x & Q));
}

static inline
FQ_ELEM fq_sub(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + Q - y);
}

static inline
FQ_ELEM fq_mul(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_red(((FQ_DOUBLEPREC)x) *(FQ_DOUBLEPREC)y);
}

static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
      return fq_cond_sub(x + y);
}

/// NOTE: non constant-time. Dont use for anything important
/// \param x[in]: < 127
/// \param y[in]: < 127
/// \return x * y;
static inline
FQ_ELEM fq_mul_non_ct(const FQ_ELEM x, const FQ_ELEM y) {
    return __fq127_lookup_table[x*128 + y];
}

/// NOTE: maybe dont use it for sensitive data
static const uint8_t __fq127_inv_table[128] __attribute__((aligned(64))) = { 0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};

// just to a look up
static inline
FQ_ELEM fq_inv(FQ_ELEM x) {
   return __fq127_inv_table[x];
} /* end fq_inv */

/* Sampling functions from the global TRNG state */

DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)

DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

/* Sampling functions from the taking the PRNG state as a parameter*/
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)

DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)

#include "macro.h"

/// \retuns ptr *b
static inline __m256i avx_mul(const uint8_t *ptr, const __m256i b) {
    vec256_t v1, v2, w1, w2;
    const __m256i c7f = _mm256_set1_epi8((short)127);
    const __m256i c7f2 = _mm256_set1_epi16((short)127);
    const __m128i t1 = _mm_loadu_si128((const __m128i *)(ptr +  0));
    const __m128i t2 = _mm_loadu_si128((const __m128i *)(ptr + 16));

    __m256i a_lo = _mm256_cvtepu8_epi16(t1);
    __m256i a_hi = _mm256_cvtepu8_epi16(t2);

    vmul_lo16(a_lo, a_lo, b);
    vmul_lo16(a_hi, a_hi, b);
    vsr16(v1, a_lo, 7);
    vsr16(v2, a_hi, 7);
    w1 = _mm256_packs_epi16(a_lo&c7f2, a_hi&c7f2);
    w2 = _mm256_packs_epi16(v1, v2);

    vadd8(v1, w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0b11011000);
    w1 = _mm256_sub_epi8(v2, c7f);
    w2 = _mm256_blendv_epi8(w1, v2, w1);

    return w2;
}

/// NOTE: these functions are outsourced to this file, to make the
/// optimizied implementation as easy as possible.
/// accumulates a row
/// \param d
/// \return sum(d) for _ in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *d) {
    vec256_t s, t, c01, c7f;
    vset8(s, 0);
    vset8(c01, 0x01);
    vset8(c7f, 0x7F);

    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(t, (const vec256_t *)(d + col));
        vadd8(s, s, t);
        //barrett_red8(s, t, c7f, c01);
        W_RED127_(s);
	 }

    uint32_t k = vhadd8(s);
    return fq_red(k);
}

/// accumulates the inverse of a row
/// \param d
/// \return sum(d[i]**-1) for i in range(N-K)
static inline FQ_ELEM row_acc_inv(const FQ_ELEM *d) {
    // NOTE: actually only the last pos need to be 0
    static FQ_ELEM inv_data[N_K_pad] = {0};
    for (uint32_t col = 0; col < (N-K); col++) {
        inv_data[col] = fq_inv(d[col]);
    }

    return row_acc(inv_data);
}

/// accumulates the inverse of a row
/// \param d
/// \return sum(d[i]**-1) for i in range(N-K)
//static inline
//FQ_ELEM row_acc_inv2(const FQ_ELEM *d) {
//    const __m256i perm = _mm256_setr_epi32(0,1,4,5,2,3,6,7);
//    // NOTE: actually only the last pos need to be 0
//    static FQ_ELEM inv_data[N_K_pad] = {0}; 
//    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
//        const __m128i f1 = _mm_loadu_si128((const __m128i *)(d + col +  0));
//        const __m128i f2 = _mm_loadu_si128((const __m128i *)(d + col + 16));
//        const __m256i t1 = _mm256_cvtepi8_epi32(f1);
//        const __m256i t3 = _mm256_cvtepi8_epi32(f2);
//        const __m256i t2 = _mm256_cvtepi8_epi32(_mm_bsrli_si128(f1, 8));
//        const __m256i t4 = _mm256_cvtepi8_epi32(_mm_bsrli_si128(f2, 8));
//
//        const __m256i l1 = _mm256_i32gather_epi32(__fq127_inv_table, t1, 4);
//        const __m256i l2 = _mm256_i32gather_epi32(__fq127_inv_table, t2, 4);
//        const __m256i l3 = _mm256_i32gather_epi32(__fq127_inv_table, t3, 4);
//        const __m256i l4 = _mm256_i32gather_epi32(__fq127_inv_table, t4, 4);
//
//        const __m256i a1 = _mm256_packs_epi32(l1, l2);
//        const __m256i a2 = _mm256_packs_epi32(l3, l4);
//        const __m256i a3 = _mm256_permutevar8x32_epi32(a1, perm);
//        const __m256i a4 = _mm256_permutevar8x32_epi32(a2, perm);
//
//        const __m256i b1 = _mm256_packus_epi16(a3, a4);
//        const __m256i b2 = _mm256_permute4x64_epi64(b1, 0b11011000);
//        _mm256_storeu_si256((__m256i *)(inv_data + col), b2);
//	}
//
//    return row_acc(inv_data);
//}

/// scalar multiplication of a row
/// NOTE: not a full reduction
/// \param row[in/out] *= s for _ in range(N-K)
/// \param s
static inline
void row_mul(FQ_ELEM *row, const FQ_ELEM s) {
    // precompute b
    const __m256i b = _mm256_set1_epi16(s);
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        const __m256i t = avx_mul(row + col, b);
        vstore256((vec256_t *)(row + col), t);
    }
}

/// scalar multiplication of a row
/// NOTE: not a full reduction
/// \param row[in/out] *= s for _ in range(N-K)
/// \param s
static inline
void row_mul_(FQ_ELEM *row, const FQ_ELEM s) {
    vec256_t shuffle, t, c7f, c01, b, a, a_lo, a_hi, b_lo, b_hi;
    vec128_t tmp;

    vload256(shuffle, (vec256_t *) shuff_low_half);
    vset8(c7f, 127);
    vset8(c01, 1);

    // precompute b
    vset8(b, s);
    vget_lo(tmp, b);
    vextend8_16(b_lo, tmp);
    vget_hi(tmp, b);
    vextend8_16(b_hi, tmp);

    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(row + col));

        vget_lo(tmp, a);
        vextend8_16(a_lo, tmp);
        vget_hi(tmp, a);
        vextend8_16(a_hi, tmp);

        barrett_mul_u16(a_lo, a_lo, b_lo, t);
        barrett_mul_u16(a_hi, a_hi, b_hi, t);

        vshuffle8(a_lo, a_lo, shuffle);
        vshuffle8(a_hi, a_hi, shuffle);

        vpermute_4x64(a_lo, a_lo, 0xd8);
        vpermute_4x64(a_hi, a_hi, 0xd8);

        vpermute2(t, a_lo, a_hi, 0x20);

        W_RED127_(t);
        vstore256((vec256_t *)(row + col), t);
    }
}


/// scalar multiplication of a row
/// NOTE: not a full reduction
/// \param row[in/out] *= s for _ in range(N-K)
/// \param s
//static inline
//void row_mul_gather(FQ_ELEM *row, const FQ_ELEM s) {
//    const __m256i perm = _mm256_setr_epi32(0,1,4,5,2,3,6,7);
//    const uint32_t *ptr = __fq127_lookup_table + s*128;
//    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
//        const __m128i f1 = _mm_loadu_si128((const __m128i *)(row + col +  0));
//        const __m128i f2 = _mm_loadu_si128((const __m128i *)(row + col + 16));
//        const __m256i t1 = _mm256_cvtepi8_epi32(f1);
//        const __m256i t3 = _mm256_cvtepi8_epi32(f2);
//        const __m256i t2 = _mm256_cvtepi8_epi32(_mm_bsrli_si128(f1, 8));
//        const __m256i t4 = _mm256_cvtepi8_epi32(_mm_bsrli_si128(f2, 8));
//
//        const __m256i l1 = _mm256_i32gather_epi32(ptr, t1, 4);
//        const __m256i l2 = _mm256_i32gather_epi32(ptr, t2, 4);
//        const __m256i l3 = _mm256_i32gather_epi32(ptr, t3, 4);
//        const __m256i l4 = _mm256_i32gather_epi32(ptr, t4, 4);
//
//        const __m256i a1 = _mm256_packs_epi32(l1, l2);
//        const __m256i a2 = _mm256_packs_epi32(l3, l4);
//        const __m256i a3 = _mm256_permutevar8x32_epi32(a1, perm);
//        const __m256i a4 = _mm256_permutevar8x32_epi32(a2, perm);
//
//        const __m256i b1 = _mm256_packus_epi16(a3, a4);
//        const __m256i b2 = _mm256_permute4x64_epi64(b1, 0b11011000);
//        _mm256_storeu_si256((__m256i *)(row + col), b2);
//	}
//}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in
/// \param s
static inline
void row_mul2(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    vec256_t shuffle, t, c7f, c01, b, a, a_lo, a_hi, b_lo, b_hi;
    vec128_t tmp;

    vload256(shuffle, (vec256_t *) shuff_low_half);
    vset8(c7f, 127);
    vset8(c01, 1);

    // precompute b
    vset8(b, s);
    vget_lo(tmp, b);
    vextend8_16(b_lo, tmp);
    vget_hi(tmp, b);
    vextend8_16(b_hi, tmp);

    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in + col));

        vget_lo(tmp, a);
        vextend8_16(a_lo, tmp);
        vget_hi(tmp, a);
        vextend8_16(a_hi, tmp);

        barrett_mul_u16(a_lo, a_lo, b_lo, t);
        barrett_mul_u16(a_hi, a_hi, b_hi, t);

        vshuffle8(a_lo, a_lo, shuffle);
        vshuffle8(a_hi, a_hi, shuffle);

        vpermute_4x64(a_lo, a_lo, 0xd8);
        vpermute_4x64(a_hi, a_hi, 0xd8);

        vpermute2(t, a_lo, a_hi, 0x20);

        // barrett_red8(t, r, c7f, c01);
        W_RED127_(t);
        vstore256((vec256_t *)(out + col), t);
    }
}

///
/// \param out = in1[i]*in2[i] for i in range(N-K)
/// \param in1
/// \param in2
static inline
void row_mul3(FQ_ELEM *out, const FQ_ELEM *in1, const FQ_ELEM *in2) {
    vec256_t shuffle, t, c7f, c01, a, a_lo, a_hi, b, b_lo, b_hi;
    vec128_t tmp;

    vload256(shuffle, (vec256_t *) shuff_low_half);
    vset8(c7f, 127);
    vset8(c01, 1);

    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in1 + col));
        vload256(b, (vec256_t *)(in2 + col));

        vget_lo(tmp, a);
        vextend8_16(a_lo, tmp);
        vget_hi(tmp, a);
        vextend8_16(a_hi, tmp);
        vget_lo(tmp, b);
        vextend8_16(b_lo, tmp);
        vget_hi(tmp, b);
        vextend8_16(b_hi, tmp);

        barrett_mul_u16(a_lo, a_lo, b_lo, t);
        barrett_mul_u16(a_hi, a_hi, b_hi, t);

        vshuffle8(a_lo, a_lo, shuffle);
        vshuffle8(a_hi, a_hi, shuffle);

        vpermute_4x64(a_lo, a_lo, 0xd8);
        vpermute_4x64(a_hi, a_hi, 0xd8);

        vpermute2(t, a_lo, a_hi, 0x20);

        // barrett_red8(t, r, c7f, c01);
        W_RED127_(t);
        vstore256((vec256_t *)(out + col), t);
    }
}

/// invert a row
/// \param out[out]: in[i]**-1 for i in range(N-K)
/// \param in [in]
static inline
void row_inv2(FQ_ELEM *out, const FQ_ELEM *in) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_inv(in[col]);
    }
}

/// NOTE: avx512 optimized version (it has a special instruction for stuff like this)
/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
    vec256_t t1, t2, acc;
    vset8(acc, -1);
    vset8(t2, in[0]);

    uint32_t col = 0;
    for (; col < N_K_pad-32; col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, t2);
        vand(acc, acc, t1);
    }

    const uint32_t t3 = (uint32_t)vmovemask8(acc);
    for (;col < N-K; col++) {
        if (in[col-1] != in[col]) {
            return 0;
        }
    }

    return t3 == -1u;
}

/// \param in[in] row
/// \return 1 if a zero was found
///         0 else
static inline
uint32_t row_contains_zero(const FQ_ELEM *in) {
    vec256_t t1, t2, acc;
    vset8(t2, 0);
    vset8(acc, 0);
    uint32_t col = 0;
    for (; col < (N_K_pad-32); col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, t2);
        vor(acc, acc, t1);
    }
    
    const uint32_t t3 = (uint32_t)vmovemask8(acc);
    if (t3 != 0ul) { return 1; }

    for (;col < N-K; col++) {
        if (in[col] == 0) {
            return 1;
        }
    }
    return 0;
}

/// \param in[in]: vector of length N-K
/// \return the number of zeros in the input vector
static inline
uint32_t row_count_zero(const FQ_ELEM *in) {
    vec256_t t1, zero, acc, mask;
    vset8(zero, 0);
    vset8(acc, 0);
    vset8(mask, 1);
    uint32_t col = 0;
    for (; col < (N_K_pad-32); col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, zero);
        vand(t1, t1, mask);
        vadd8(acc, acc, t1);
    }
  
    uint32_t a = vhadd8(acc);
    for (;col < N-K; col++) {
        a += (in[col] == 0);
    }
    return a;
}

