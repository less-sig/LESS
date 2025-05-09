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
#include <immintrin.h>

#include "parameters.h"
#include "rng.h"
#include "lookup_table.h"

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


/// constant itme implementation
/// \param x[in] < 256
/// \return x mod 127
static inline
FQ_ELEM fq_cond_sub(const FQ_ELEM x) {
    // equivalent to: (x >= Q) ? (x - Q) : x
    // likely to be ~ constant-time (a "smart" compiler might turn this into conditionals though)
    FQ_ELEM sub_q = x - Q;
    FQ_ELEM mask = -(sub_q >> NUM_BITS_Q);
    return (mask & Q) + sub_q;
} /* end fq_cond_sub */

/// constant itme implementation
/// \param x[in] < 127**2
/// \return x mod 127
static inline
FQ_ELEM fq_red(const FQ_DOUBLEPREC x) {
    return fq_cond_sub((x >> NUM_BITS_Q) + ((FQ_ELEM) x & Q));
} /* end fq_red */

/// constant itme implementation
/// \param x[in] < 127
/// \param y[in] < 127
/// \return x - y mod 127
static inline
FQ_ELEM fq_sub(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + Q - y);
} /* end fq_sub */

/// constant itme implementation
/// \param x[in] < 127
/// \param y[in] < 127
/// \return x * y mod 127
static inline
FQ_ELEM fq_mul(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_red(((FQ_DOUBLEPREC)x) *(FQ_DOUBLEPREC)y);
} /* end fq_mul */

/// constant itme implementation
/// \param x[in] < 127
/// \param y[in] < 127
/// \return x + y mod 127
static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + y);

} /* end fq_add */

/// NOTE: non constant-time. Don't use for anything important
/// \param x[in]: < 127
/// \param y[in]: < 127
/// \return x * y;
static inline
FQ_ELEM fq_mul_non_ct(const FQ_ELEM x, const FQ_ELEM y) {
    return __gf127_lookuptable[x*128 + y];
} /* end fq_mul_non_ct */

/// NOTE: maybe dont use it for sensitive data
static const uint8_t fq_inv_table[128] __attribute__((aligned(64))) = {
   0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};


/// NOTE: non constant-time. Dont use for anything important.
/// NOTE: input must be reduced.
/// \param x[in]: input value < 127
/// \return x^{-1} mod 127
static inline
FQ_ELEM fq_inv(FQ_ELEM x) {
    return fq_inv_table[x];
} /* end fq_inv */

/// Sampling functions from the global TRNG state
DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)
DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

// Sampling functions from the taking the PRNG state as a parameter
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)
DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)

// load avx2 and avx512 macros
#include "macro.h"


/// \param ret[out]: two __m512i register, which will contain a lookup table.
///     This table can be efficently processed via `_mm512_permutex2var_epi8`, 
///     which will then contain the result of the scalar multiplication of each 
///     input with the scalar `a`.
/// \param a[in]: input scalar
static inline void gf127v_scalar_u512_compute_table(__m512i *ret,
                                                    const uint8_t a) {
    ret[0] = _mm512_load_si512((const __m512i *)(__gf127_lookuptable + 128 * a +  0));
    ret[1] = _mm512_load_si512((const __m512i *)(__gf127_lookuptable + 128 * a + 64));
}

/// \param in[in]: input register: containing 64 Fq elements:
///     [in_0, ..., in_63]
/// \param table1[in]: lower 64 byte of the lookup table computed by `gf127v_scalar_u512_compute_table`
/// \param table2[in]: upper 64 byte of the lookup table computed by `gf127v_scalar_u512_compute_table`
/// \return in*a = [in_0*a, ..., in_63*a], where `a` is the scalar input of 
///     `gf127v_scalar_u512_compute_table`
static inline
__m512i gf127v_scalar_table_u512(const __m512i in,
                                 const __m512i table1,
                                 const __m512i table2) {
    return _mm512_permutex2var_epi8(table1, in, table2);
}

/// \param in[in]: avx2 register
/// \return (in[0] + in[1] + ... + in[31]) % q
static inline 
uint8_t vhadd8(const __m256i in) {
    // TODO replace  with: _mm256_reduce_add_epi16
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

/// \param in[in]: avx512 register
/// \return (in[0] + in[1] + ... + in[63]) % q
static inline 
uint8_t vhadd8_512(const __m512i in) {
    vec256_t x, t, c7f;
    vset8(c7f, 0x7F);
    const __m256i low  = _mm512_castsi512_si256   (in);
    const __m256i high = _mm512_extracti32x8_epi32(in, 1);
    vadd8(x, low, high);
    vred8(x, t, c7f);
    return vhadd8(x);
} /* end vhadd8_512 */

/// \param aa[in]: aa[i] < 127 in 16 bit limb for i in range(32)
/// \param bb[in]: bb[i] < 127 in 16 bit limb for i in range(32)
/// \return aa[i] * bb[i] for all i in range(32)
static inline __m256i avx_mul_full512(const __m512i aa,
                                      const __m512i bb) {

    __m512i c01, c7f, c516, tmp;
    vset16_512(c01, 0x01);
    vset16_512(c7f, 0x7F);
    vset16_512(c516, 516);

    __m512i acc = _mm512_mullo_epi16(aa, bb);
    tmp = _mm512_add_epi16(acc, c01);
    tmp = _mm512_mulhi_epu16(tmp, c516);
    tmp = _mm512_mullo_epi16(tmp, c7f);
    acc = _mm512_sub_epi16(acc, tmp);
    const __m256i t = _mm512_cvtepi16_epi8(acc);
    return t;
}

/// accumulates a row
/// \param row[in]: pointer to a row with ROUND_UP(N-K, 32) elements
/// \return sum(d[i]) for i in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *row) {
    __m512i c01, c7f, c516;
    __m256i t;
    vset16_512(c01, 0x01);
    vset16_512(c7f, 0x7F);
    vset16_512(c516, 516);
    __m512i acc = _mm512_setzero_si512();
    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(t, (const vec256_t *)(row + col));
        const __m512i a = _mm512_cvtepi8_epi16(t);
        acc = _mm512_add_epi16(acc, a);
	}

    __m512i tmp;
    tmp = _mm512_add_epi16(acc, c01);
    tmp = _mm512_mulhi_epu16(tmp, c516);
    tmp = _mm512_mullo_epi16(tmp, c7f);
    acc = _mm512_sub_epi16(acc, tmp);

    const __m256i a2 = _mm512_extracti64x4_epi64(acc, 1);
    const __m256i a1 = _mm512_castsi512_si256(acc) ;
    const __m256i a = _mm256_add_epi16(a1, a2);
    const short k = _mm256_reduce_add_epi16(a);
    return fq_red(k);
}

/// accumulates the inverse of a row
/// \param row[in]: pointer to a row with ROUND_UP(N-K, 32) elements
/// \return sum(d[i]**-1) for i in range(N-K)
static inline
FQ_ELEM row_acc_inv(const FQ_ELEM *row) {
    const __m512i t1 = _mm512_load_si512((const __m512i *)(fq_inv_table +  0));
    const __m512i t2 = _mm512_load_si512((const __m512i *)(fq_inv_table + 64));

    static FQ_ELEM inv_data[N_K_pad] __attribute__((aligned(64))) = {0}; 

    uint32_t col = 0;
    for (; (col+64) <= N_K_pad; col += 64) {
        const __m512i a = _mm512_loadu_si512((const __m512i *)(row + col));
        const __m512i k = _mm512_permutex2var_epi8(t1, a, t2);
        _mm512_store_si512((__m512i *)(inv_data + col), k);
	}

    // tail mngt...
    for (; (col+32) <= N_K_pad; col += 32) {
        const __m256i b = _mm256_loadu_si256((const __m256i *)(row + col));
        const __m512i a = _mm512_castsi256_si512(b);
        const __m512i c = gf127v_scalar_table_u512(a, t1, t2); 
        const __m256i d = _mm512_castsi512_si256(c);
        _mm256_storeu_si256((__m256i *)(inv_data + col), d);
    }

    return row_acc(inv_data);
}

/// scalar multiplication of a row
/// \param row[in/out] *= s for _ in range(N-K)
/// \param s[in]: input scalar
static inline
void row_mul(FQ_ELEM *row,
             const FQ_ELEM s) {
    __m512i table[2];
    gf127v_scalar_u512_compute_table(table, s);

    uint32_t col = 0;
    for (; (col+64) <= N_K_pad; col += 64) {
        const __m512i a = _mm512_loadu_si512((const __m512i *)(row + col));
        const __m512i b = gf127v_scalar_table_u512(a, table[0], table[1]); 
        _mm512_storeu_si512((__m512i *)(row + col), b);
	}

    for (; (col+32) <= N_K_pad; col += 32) {
        const __m256i b = _mm256_loadu_si256((const __m256i *)(row + col));
        const __m512i a = _mm512_castsi256_si512(b);
        const __m512i c = gf127v_scalar_table_u512(a, table[0], table[1]); 
        const __m256i d = _mm512_castsi512_si256(c);
        _mm256_storeu_si256((__m256i *)(row + col), d);
    }
}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in[in]: input row with ROUND_UP(N-K, 32) elements
/// \param s[in]: input scalar
static inline
void row_mul2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in,
              const FQ_ELEM s) {
    __m512i table[2];
    gf127v_scalar_u512_compute_table(table, s);

    uint32_t col = 0;
    for (; (col+64) <= N_K_pad; col += 64) {
        const __m512i a = _mm512_loadu_si512((const __m512i *)(in + col));
        const __m512i b = gf127v_scalar_table_u512(a, table[0], table[1]); 
        _mm512_storeu_si512((__m512i *)(out + col), b);
    }

    for (; (col+32) <= N_K_pad; col += 32) {
        const __m256i b = _mm256_loadu_si256((const __m256i *)(in + col));
        const __m512i a = _mm512_castsi256_si512(b);
        const __m512i c = gf127v_scalar_table_u512(a, table[0], table[1]); 
        const __m256i d = _mm512_castsi512_si256(c);
        _mm256_storeu_si256((__m256i *)(out + col), d);
    }
}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in[in]: input row with ROUND_UP(N-K, 32) elements
/// \param s[in]: input scalar
static inline
void row_mul2_ct(FQ_ELEM *__restrict__ out,
                 const FQ_ELEM *__restrict__ in, 
                 const FQ_ELEM s) {
    __m512i c01, c7f, c516, tmp, bb;
    __m256i a;
    vset16_512(c01, 0x01);
    vset16_512(c7f, 0x7F);
    vset16_512(c516, 516);
    vset16_512(bb, s);

    // TODO optimize
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in + col));
        const __m512i aa = _mm512_cvtepi8_epi16(a);
        __m512i acc = _mm512_mullo_epi16(aa, bb);
        tmp = _mm512_add_epi16(acc, c01);
        tmp = _mm512_mulhi_epu16(tmp, c516);
        tmp = _mm512_mullo_epi16(tmp, c7f);
        acc = _mm512_sub_epi16(acc, tmp);
        const __m256i t = _mm512_cvtepi16_epi8(acc);
        vstore256((__m256i *)(out + col), t);
    }
}

/// \param out = in1[i]*in2[i] for i in range(N-K)
/// \param in1[in]: input row with ROUND_UP(N-K, 32) elements
/// \param in2[in]: input row with ROUND_UP(N-K, 32) elements
static inline
void row_mul3(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in1,
              const FQ_ELEM *__restrict__ in2) {
    __m512i c01, c7f, c516, tmp;
    __m256i a, b;
    vset16_512(c01, 0x01);
    vset16_512(c7f, 0x7F);
    vset16_512(c516, 516);

    // TODO optimize
    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in1 + col));
        vload256(b, (vec256_t *)(in2 + col));
        const __m512i aa = _mm512_cvtepi8_epi16(a);
        const __m512i bb = _mm512_cvtepi8_epi16(b);
        __m512i acc = _mm512_mullo_epi16(aa, bb);
        tmp = _mm512_add_epi16(acc, c01);
        tmp = _mm512_mulhi_epu16(tmp, c516);
        tmp = _mm512_mullo_epi16(tmp, c7f);
        acc = _mm512_sub_epi16(acc, tmp);
        const __m256i t = _mm512_cvtepi16_epi8(acc);
        vstore256((__m256i *)(out + col), t);
    }
}

/// invert a row
/// \param out[out]: in[i]**-1 for i in range(N-K)
/// \param in[in]: input row with ROUND_UP(N-K, 32) elements
static inline
void row_inv2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in) {
    const __m512i t1 = _mm512_load_si512((const __m512i *)(fq_inv_table +  0));
    const __m512i t2 = _mm512_load_si512((const __m512i *)(fq_inv_table + 64));

    uint32_t col = 0;
    for (; (col+64) <= N_K_pad; col += 64) {
        const __m512i a = _mm512_loadu_si512((const __m512i *)(in + col));
        const __m512i k = _mm512_permutex2var_epi8(t1, a, t2);

        _mm512_storeu_si512((__m512i *)(out + col), k);
	}

    for (; (col+32) <= N_K_pad; col += 32) {
        const __m256i b = _mm256_loadu_si256((const __m256i *)(in + col));
        const __m512i a = _mm512_castsi256_si512(b);
        const __m512i k = _mm512_permutex2var_epi8(t1, a, t2);
        const __m256i c = _mm512_castsi512_si256(k);
        _mm256_storeu_si256((__m256i *)(out + col), c);
    }
}

/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
    vec256_t t1, t2, acc;
    vset8(acc, -1u);
    vset8(t2, in[0]);

    // TODO make 64 byte wide
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

/// \param in[in]: row of length N-K
/// \return 1 if a zero was found
///         0 else
static inline
uint32_t row_contains_zero(const FQ_ELEM *in) {
    const __m512i zero = _mm512_setzero_si512();
    uint32_t col = 0;
    __mmask64 acc = 0;

    for (; (col+64) <= N_K_pad; col += 64) {
        const __m512i a = _mm512_loadu_si512((const __m512i *)(in + col));
        const __mmask64 b = _mm512_cmpeq_epu8_mask(a, zero);
        acc |= b;
	}

    for (; (col+32) <= N_K_pad; col += 32) {
        const __m256i zero2 = _mm512_castsi512_si256(zero);
        const __m256i a = _mm256_loadu_si256((const __m256i *)(in + col));
        const __mmask64 b = _mm256_cmpeq_epu8_mask(a, zero2);
        acc |= b;
    }

    return acc != 0;
}

/// \param in[in]: row of length N-K
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

    acc =_mm256_sad_epu8(acc, zero);
    t1 = _mm256_srli_si256(acc, 8);
    acc = _mm256_add_epi64(acc, t1);

    t1 = _mm256_permute2x128_si256(acc, acc, 1);
    acc = _mm256_add_epi8(acc, t1);

    uint32_t a = _mm256_extract_epi32(acc, 0);
    for (;col < N-K; col++) {
        a += (in[col] == 0);
    }
    return a;
}

// TODO use this function and catch the case that +=32
// static inline
// uint32_t row_count_zero(const FQ_ELEM *in) {
//     const __m512i mask = _mm512_set1_epi8(1);
//     const __m512i zero = _mm512_setzero_si512();
//     __m512i acc = _mm512_setzero_si512();
//
//     uint32_t col = 0;
//     for (; (col+64) <= N_K_pad; col += 64) {
//         const __m512i a = _mm512_loadu_si512((const __m512i *)(in + col));
//         const __mmask64 b = _mm512_cmpeq_epu8_mask(a, zero);
//         acc = _mm512_maskz_add_epi8(b, acc, mask);
// 	}
//
//     const __m512i t = _mm512_sad_epu8(acc, zero);
//     const int64_t tt = _mm512_reduce_add_epi64(t);
//     return tt;
// }
