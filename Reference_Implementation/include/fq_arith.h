/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (March 2025)
 *
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

#include "parameters.h"
#include "rng.h"
#include "utils.h"
#include "lookup_table.h"

/// number of bits needed to represent Q = 7
#define NUM_BITS_Q (BITS_TO_REPRESENT(Q))

/// macro generating functions to generate random numbers
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
      for (unsigned i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
         EL_T rnd_value = word & EL_MASK; \
         if (rnd_value <= SPAN) buffer[count++] = rnd_value + MIN_VALUE; \
         if (count >= num_elements) return; \
         word >>= REQ_BITS; \
      } \
   } while (1); }

/// macro generating functions to generate random numbers
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
      for (unsigned i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
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

/// NOTE: these functions are outsourced to this file, to make the
/// optimized implementation as easy as possible.
/// accumulates a row
/// \param d[in]: input row
/// \return sum(d)
static inline
FQ_ELEM row_acc(const FQ_ELEM *d) {
    FQ_ELEM s = d[0];
    for (uint32_t col = 1; col < N-K; col++) {
        s = fq_add(s, d[col]);
	 }

    return s;
}

/// NOTE: not constant-time, as fq_inv_non_ct is used.
/// accumulates the inverse of a row
/// \param d[in]: input row
/// \return sum([d[i]^{-1} for i in range(N-K)])
static inline
FQ_ELEM row_acc_inv(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < N-K; col++) {
        s = fq_add(s, fq_inv_non_ct(d[col]));
	 }

    return s;
}

/// NOTE: constant time implementation
/// scalar multiplication of a row
/// \param row[in/out] *= s for i in range(N-K)
/// \param s[in]: scalar value < 127
static inline
void row_mul(FQ_ELEM *row, const FQ_ELEM s) {
    for (uint32_t col = 0; col < N-K; col++) {
        row[col] = fq_mul(s, row[col]);
    }
}

/// NOTE: constant time implementation
/// scalar multiplication of a row
/// \param out[out] = s*in[i] for i in range(N-K)
/// \param in[in]: input row
/// \param s[in]: scalar value < 127
static inline
void row_mul2(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    for (uint32_t col = 0; col < N-K; col++) {
        out[col] = fq_mul(s, in[col]);
    }
}

/// NOTE: constant time implementation
/// \param out[out] = in1[i]*in2[i] for i in range(N-K)
/// \param in1[in]: first input row 
/// \param in2[in]: second input row
static inline
void row_mul3(FQ_ELEM *out, const FQ_ELEM *in1, const FQ_ELEM *in2) {
    for (uint32_t col = 0; col < N-K; col++) {
        out[col] = fq_mul(in1[col], in2[col]);
    }
}

/// NOTE: non ct, as fq_inv_non_ct is used.
/// invert a row
/// \param out[out]: = in[i]**{-1} for i in range(N-K)
/// \param in[in]: input row of length N-K
static inline
void row_inv2(FQ_ELEM *out, const FQ_ELEM *in) {
    for (uint32_t col = 0; col < N-K; col++) {
        out[col] = fq_inv_non_ct(in[col]);
    }
}

/// NOTE: not constant time.
/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
    for (uint32_t col = 1; col < N-K; col++) {
        if (in[col] != in[col - 1]) {
            return 0;
        }
    }
    return 1;
}

/// NOTE: not constant time.
/// \param in[in]: vector of length N-K
/// \return 0 if no zero was found
///         1 if the row contains at least a single 0
static inline
uint32_t row_contains_zero(const FQ_ELEM *in) {
    for (uint32_t col = 0; col < N-K; col++) {
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
    uint32_t r = 0;
    for (uint32_t col = 0; col < N-K; col++) {
        r += in[col] == 0;
    }
    return r;
}
