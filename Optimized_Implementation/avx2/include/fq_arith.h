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
#include "utils.h"
#include "lookup_table.h"
#include "rng.h"

// number of needed to represent q = 7
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


/// Sampling functions from the global TRNG state
DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)
DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

/// Sampling functions from the taking the PRNG state as a parameter
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)
DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)


/// constant time implementation
/// \param x[in]: input number < 2*127
/// \return x mod 127
static inline
FQ_ELEM fq_cond_sub(const FQ_ELEM x) {
    // equivalent to: (x >= Q) ? (x - Q) : x
    // likely to be ~ constant-time (a "smart" compiler might turn this into conditionals though)
    FQ_ELEM sub_q = x - Q;
    FQ_ELEM mask = -(sub_q >> NUM_BITS_Q);
    return (mask & Q) + sub_q;
}

/// constant time implementation
/// \param x[in]: input number < 127**2
/// \return x mod 127
static inline
FQ_ELEM fq_red(const FQ_DOUBLEPREC x) {
    return fq_cond_sub((x >> NUM_BITS_Q) + ((FQ_ELEM) x & Q));
}

/// constant time implementation
/// \param x[in]: minuend < 127
/// \param y[in]: subtrahend < 127
/// \return difference x-y mod 127
static inline
FQ_ELEM fq_sub(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + Q - y);
}

/// constant time implementation
/// \param x[in]: first summand < 127
/// \param y[in]: second summand < 127
/// \return sum x+y mod 127
static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + y);
}

/// constant time implementation
/// \param x[in]: first factor < 127
/// \param y[in]: second factor < 127
/// \return product x*y mod 127
static inline
FQ_ELEM fq_mul(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_red((FQ_DOUBLEPREC) x * (FQ_DOUBLEPREC) y);
}

/// NOTE: non constant-time. Don't use it for anything important
/// \param x[in]: first factor < 127
/// \param y[in]: second factor < 127
/// \return product x*y mod 127;
static inline
FQ_ELEM fq_mul_non_ct(const FQ_ELEM x, const FQ_ELEM y) {
    return __fq127_lookup_table[x*128 + y];
}

/// NOTE: non constant-time. Don't use it for anything important
static const uint8_t fq_inv_table[128] __attribute__((aligned(64))) = {
   0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};

/// NOTE: non constant-time. Don't use for anything important.
/// NOTE: input must be reduced.
/// \param x[in]: input value < 127
/// \return x^{-1} mod 127
static inline
FQ_ELEM fq_inv_non_ct(const FQ_ELEM x) {
   return fq_inv_table[x];
}

/// NOTE: constant time
/// NOTE: input must be reduced.
/// \param x[in]: input value < 127
/// \return x^{-1} mod 127
static inline
FQ_ELEM fq_inv(const FQ_ELEM x) {
#if 1
    FQ_ELEM ret = 0;
    for (uint64_t i = 0; i < 127; i++) {
        const FQ_ELEM mask = COMPUTE_CT_MASK(i, x);
        const FQ_ELEM val = fq_inv_table[i] & mask;
        ret ^= val;
    }
    return ret;
#else
    /// Fermat's method for inversion employing r-t-l square and multiply,
    /// unrolled for actual parameters */
    FQ_DOUBLEPREC xlift = x;
    FQ_DOUBLEPREC accum = 1;
    // No need for square and mult always, Q-2 is public
    uint32_t exp = Q-2;
    while(exp) {
        if(exp & 1) {
            accum = fq_red(accum*xlift);
        }
        xlift = fq_red(xlift*xlift);
        exp >>= 1;
    }
    return fq_red(accum);
#endif
} /* end fq_inv */

/// load avx2 macros
#include "macro.h"


/// \param in[in]: avx2 register
/// \return in[0] + in[1] + ... + in[31] % q
static inline 
uint8_t vhadd8(const __m256i in) {
    vec256_t c7f, tmp;
    vset8(c7f, 0x7F);

    __m256i a = _mm256_srli_epi16(in, 8);
    __m256i t = _mm256_add_epi8(a, in);
    vred8(t,tmp,c7f)

    a = _mm256_srli_epi32(t, 16);
    t = _mm256_add_epi8(a, t);
    vred8(t,tmp,c7f)

    a = _mm256_srli_epi64(t, 32);
    t = _mm256_add_epi8(a, t);
    vred8(t,tmp,c7f)

    a = _mm256_srli_si256(t, 8);
    t = _mm256_add_epi8(a, t);
    vred8(t,tmp,c7f)

    a = _mm256_permute2x128_si256(t, t, 1);
    t = _mm256_add_epi8(a, t);
    vred8(t,tmp,c7f)

    return _mm256_extract_epi8(t, 0);
} /* end vhadd8 */

/// \retuns ptr *b
/// \param ptr[in]: pointer to 32 byte of memory, can be unaligned.
///     NOTE: the memory does not need to be aligned
/// \param b[in]: avx2 register containing 16 already zero extended Fq elements.
///     [b_0, ..., b_15], b_i \in Fq and b_i = uint16_t
static inline __m256i avx_mul(const uint8_t *ptr,
                              const __m256i b) {
    vec256_t v1, v2;
    const __m256i c7f = _mm256_set1_epi8((short)0x7f);
    const __m256i c7f2 = _mm256_set1_epi16((short)0x7f);
    const __m128i t1 = _mm_loadu_si128((const __m128i *)(ptr +  0));
    const __m128i t2 = _mm_loadu_si128((const __m128i *)(ptr + 16));

    __m256i a_lo = _mm256_cvtepu8_epi16(t1);
    __m256i a_hi = _mm256_cvtepu8_epi16(t2);

    vmul_lo16(a_lo, a_lo, b);
    vmul_lo16(a_hi, a_hi, b);
    vsr16(v1, a_lo, 7);
    vsr16(v2, a_hi, 7);
    __m256i w1 = _mm256_packs_epi16(_mm256_and_si256(a_lo, c7f2), _mm256_and_si256(a_hi, c7f2));
    __m256i w2 = _mm256_packs_epi16(v1, v2);

    vadd8(v1, w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0xd8); // 0b11011000
    w1 = _mm256_sub_epi8(v2, c7f);
#if defined(LESS_USE_BLEND_IN_ARITH)
    w2 = _mm256_blendv_epi8(w1, v2, w1);
#else
    w2 = _mm256_min_epu8(v2, w1);
#endif
    return w2;
} /* end avx_mul */

/// each of the input SSE registers is in a 8bit limb
/// \param a1[in]: [a_0, ..., a_15], a_i \in Fq, lower values
/// \param a2[in]: [a_0, ..., a_15], a_i \in Fq, upper values
/// \param b1[in]: [b_0, ..., b_15], b_i \in Fq, lower values
/// \param b2[in]: [b_0, ..., b_15], b_i \in Fq, upper values
///         0         a1*b1          15, 16      a2*b2        31 (limb)
/// \return [ a_0 * b_0, ..., a_15*b_15, a_0*b_0, ..., a_15*b_15]
static inline __m256i avx_mul_full(const __m128i a1, const __m128i a2,
                                   const __m128i b1, const __m128i b2) {
    vec256_t v1, v2;
    const __m256i c7f = _mm256_set1_epi8((short)127);
    const __m256i c7f2 = _mm256_set1_epi16((short)127);

    __m256i a_lo = _mm256_cvtepu8_epi16(a1);
    __m256i a_hi = _mm256_cvtepu8_epi16(a2);
    const __m256i b_lo = _mm256_cvtepu8_epi16(b1);
    const __m256i b_hi = _mm256_cvtepu8_epi16(b2);

    vmul_lo16(a_lo, a_lo, b_lo);
    vmul_lo16(a_hi, a_hi, b_hi);
    vsr16(v1, a_lo, 7);
    vsr16(v2, a_hi, 7);
    __m256i w1 = _mm256_packs_epi16(_mm256_and_si256(a_lo, c7f2), _mm256_and_si256(a_hi, c7f2));
    __m256i w2 = _mm256_packs_epi16(v1, v2);

    vadd8(v1, w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0xd8); // 0b11011000
    w1 = _mm256_sub_epi8(v2, c7f);
#if defined(LESS_USE_BLEND_IN_ARITH)
    w2 = _mm256_blendv_epi8(w1, v2, w1);
#else
    w2 = _mm256_min_epu8(v2, w1);
#endif
    return w2;
} /* end avx_mul_full */

/// NOTE assumes each Fq element is in a 16bit limb
/// \param a1[in]: [a_0, ..., a_15], a_i \in Fq, a_i \in 16bit limb lower values
/// \param a2[in]: [a_0, ..., a_15], a_i \in Fq, a_i \in 16bit limb upper values
/// \param b1[in]: [b_0, ..., b_15], b_i \in Fq, b_i \in 16bit limb lower values
/// \param b2[in]: [b_0, ..., b_15], b_i \in Fq, b_i \in 16bit limb upper values
/// NOTE: the return value is compress from 16bits down to 8bits.
///         0         a1*b1          15, 16      a2*b2        31 (limb)
/// \return [ a_0 * b_0, ..., a_15*b_15, a_0*b_0, ..., a_15*b_15]
static inline __m256i avx_mul_full256(const __m256i a1, const __m256i a2,
                                      const __m256i b1, const __m256i b2) {
    vec256_t v1, v2;
    const __m256i c7f = _mm256_set1_epi8((short)127);
    const __m256i c7f2 = _mm256_set1_epi16((short)127);

    __m256i a_lo = a1;
    __m256i a_hi = a2;
    __m256i b_lo = b1;
    __m256i b_hi = b2;

    vmul_lo16(a_lo, a_lo, b_lo);
    vmul_lo16(a_hi, a_hi, b_hi);
    vsr16(v1, a_lo, 7);
    vsr16(v2, a_hi, 7);
    __m256i w1 = _mm256_packs_epi16(_mm256_and_si256(a_lo, c7f2), _mm256_and_si256(a_hi, c7f2));
    __m256i w2 = _mm256_packs_epi16(v1, v2);

    vadd8(v1, w1, w2);
    v2 = _mm256_permute4x64_epi64(v1, 0xd8); // 0b11011000
    w1 = _mm256_sub_epi8(v2, c7f);
#if defined(LESS_USE_BLEND_IN_ARITH) 
    w2 = _mm256_blendv_epi8(w1, v2, w1);
#else
    w2 = _mm256_min_epu8(v2, w1);
#endif
    return w2;
} /* end avx_mul_full256 */

/// NOTE: these functions are outsourced to this file, to make the optimized
/// implementation as easy as possible 
/// NOTE: the length of all input rows must be multiple of 32.
/// accumulates the input row
/// \param row[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \return sum(d[i]) for i in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *row) {
    vec256_t s, t, c7f;
    vset8(s, 0);
    vset8(c7f, 0x7F);
    
    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(t, (const vec256_t *)(row + col));
        vadd8(s, s, t);
        vred8(s, t, c7f);
	 }

    uint32_t k = vhadd8(s);
    return fq_red(k); // TODO remove 
} /* end row_acc */

/// accumulates the inverse to a row of length N_pad
/// \param row[in]: pointer to a row
/// \return sum(d[i]**-1) for i in range(N-K)
static inline FQ_ELEM row_acc_inv(const FQ_ELEM *row) {
    // NOTE: actually only the last pos need to be 0
    static FQ_ELEM inv_data[N_K_pad] = {0};
    for (uint32_t col = 0; col < (N-K); col++) {
        inv_data[col] = fq_inv_non_ct(row[col]);
    }

    return row_acc(inv_data);
} /* end row_acc_inv */

/// TODO rename to scalar mul
/// scalar multiplication of a row
/// \param row[in/out] *= s for _ in range(N-K), pointer to a row of length ROUND_UP(N-K, 32)
/// \param s[in]: scalar value
static inline
void row_mul(FQ_ELEM *row,
             const FQ_ELEM s) {
    const __m256i b = _mm256_set1_epi16(s);
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        const __m256i t = avx_mul(row + col, b);
        vstore256((vec256_t *)(row + col), t);
    }
} /* end row_mul */

/// scalar multiplication of a row
/// \param out[out]: = s*in[i] for i in range(N-K)
/// \param in[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \param s[in]: scalar value
static inline
void row_mul2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in,
              const FQ_ELEM s) {
    const __m256i b = _mm256_set1_epi16(s);
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        const __m256i t = avx_mul(in + col, b);
        vstore256((vec256_t *)(out + col), t);
    }
} /* end row_mul2 */

/// scalar multiplication of a row
/// NOTE: constant time implementation needed for `blind`
/// \param out[out]: = s*in[i] for i in range(N-K)
/// \param in[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \param s[in]: input scalar
static inline
void row_mul2_ct(FQ_ELEM *__restrict__ out,
                 const FQ_ELEM *__restrict__ in,
                 const FQ_ELEM s) {
    row_mul2(out, in, s);
} /* end row_mul2_ct */

/// full element wise multiplication of two rows
/// \param out[out]: = in1[i]*in2[i] for i in range(N-K)
/// \param in1[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \param in2[in]: pointer to a row of length ROUND_UP(N-K, 32)
static inline
void row_mul3(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in1,
              const FQ_ELEM *__restrict__ in2) {
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        const __m128i a1 = _mm_loadu_si128((const __m128i *)(in1 + col +  0));
        const __m128i a2 = _mm_loadu_si128((const __m128i *)(in1 + col + 16));
        const __m128i b1 = _mm_loadu_si128((const __m128i *)(in2 + col +  0));
        const __m128i b2 = _mm_loadu_si128((const __m128i *)(in2 + col + 16));
        const __m256i t = avx_mul_full(a1, a2, b1, b2);
        vstore256((vec256_t *)(out + col), t);
    }
} /* end row_mul3 */

/// NOTE: this is not constant time
/// NOTE: this is not improvable via AVX2
/// invert a row
/// \param out[out]: in[i]**-1 for i in range(N-K)
/// \param in [in]: pointer to a row of length N-K
static inline
void row_inv2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_inv_non_ct(in[col]);
    }
} /* end row_inv2 */

/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
#if 1
    for (uint64_t col = 0;col < N-K; col++) {
        if (in[col-1] != in[col]) {
            return 0;
        }
    }
    return 1;
#else
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
#endif
} /* end row_all_same */

/// NOTE: not CT, but close.
/// \param in[in]: row
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
} /* end row_contains_zero */

/// NOTE: eventhough the input vectors are of length ROUND_UP(N-K, 32), only 
/// N-K bytes are read. As the remaining alignment bytes are zero, which would  
/// alter the result.
/// \param in[in]: vector of length N-K
/// \return: the number of zeros in `in`
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
} /* end row_count_zero */
