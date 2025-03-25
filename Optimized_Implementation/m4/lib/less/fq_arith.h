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
#include "lookup_table.h"

#define NUM_BITS_Q (BITS_TO_REPRESENT(QQ))


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


static inline
FQ_ELEM fq_cond_sub(const FQ_ELEM x) {
    // equivalent to: (x >= Q) ? (x - Q) : x
    // likely to be ~ constant-time (a "smart" compiler might turn this into conditionals though)
    FQ_ELEM sub_q = x - QQ;
    FQ_ELEM mask = -(sub_q >> NUM_BITS_Q);
    return (mask & QQ) + sub_q;
}

static inline
FQ_ELEM fq_red(const FQ_DOUBLEPREC x) {
    return fq_cond_sub((x >> NUM_BITS_Q) + ((FQ_ELEM) x & QQ));
}

static inline
FQ_ELEM fq_sub(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + QQ - y);
}

static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_cond_sub(x + y);
}

static inline
FQ_ELEM fq_mul(const FQ_ELEM x, const FQ_ELEM y) {
    return fq_red((FQ_DOUBLEPREC) x * (FQ_DOUBLEPREC) y);
}

inline static uint32_t fq_red_u32(uint32_t Z) {
    const uint32_t mask = 0x7F7F7F7F;
    const uint32_t one = 0x01010101;
    Z = (Z & mask) + ((Z & ~mask) >> 7u);
    uint32_t C  = ((Z+one) & ~mask) ;
    return Z + (C >> 7u) - C ;
}

inline static uint32_t fq_sub_u32(const uint32_t a,
                                  const uint32_t b) {
    return fq_red_u32(a - b + 0x7F7F7F7F);
}

inline static uint32_t fq_add_u32(const uint32_t a,
                                  const uint32_t b) {
    return fq_red_u32(a + b);
}

/// @brief  
/// @param a 
/// @param b 
/// @return a-b*c
inline static uint32_t fq_scalar_sub_u32(const uint32_t a,
                                             const uint32_t b,
                                             const uint8_t c) {
    const uint32_t mask  = 0x7F7F7F7F;
    const uint32_t one   = 0x01010101;
    const uint32_t mask1 = 0x00FF00FF;
    const uint32_t mask2 = 0xFF00FF00;
    const uint32_t mask3 = 0x007F007F;
    const uint32_t mask4 = 0x7F007F00;

    const uint32_t t1 = (b&mask1) * c;
    const uint64_t t2 = ((uint64_t)b&mask2) * c;

    const uint32_t v1 = ((t1 >> 7) + (t1 & mask3)) & mask1;
    const uint32_t v2 = ((t2 >> 7) + (t2 & mask4)) & mask2;

    uint32_t Z = v1 ^ v2;
    Z = (Z & mask) + ((Z & ~mask) >> 7);
    Z = mask - Z;
    Z = Z + a;

    Z = (Z & mask) + ((Z & ~mask) >> 7);
    uint32_t C  = ((Z+one) & ~mask) ;
    return Z + (C>>7) - C ;
}


/// NOTE: non constant-time. Dont use for anything important
/// \param x[in]: < 127
/// \param y[in]: < 127
/// \return x * y;
static inline
FQ_ELEM fq_mul_non_ct(const FQ_ELEM x, const FQ_ELEM y) {
    return __fq127_lookup_table[x*128 + y];
}

/// NOTE: maybe don't use it for sensetive data
static const uint8_t fq_inv_table[128] __attribute__((aligned(64))) = {
   0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};

/// NOTE: input must be reduced, and must not be secret.
static inline
FQ_ELEM fq_inv(const FQ_ELEM x) {
   return fq_inv_table[x];
}

/* Sampling functions from the global TRNG state */
DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, QQ-1)

DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, QQ-1)

/* Sampling functions from the taking the PRNG state as a parameter*/
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, QQ-1)

DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, QQ-1)


/// NOTE: these functions are outsourced to this file, to make the
/// optimizied implementation as easy as possible.
/// accumulates a row
/// \param d
/// \return sum(d) for _ in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *d) {
    uint32_t s = 0;

    uint32_t col = 0;
    for (; (col + 4) <= (NN-K); col+=4) {
        s = fq_add_u32(s, *((uint32_t *)(d + col)));
	}

    uint8_t tmp[4];
    *((uint32_t *)tmp) = s;

    uint8_t ss = tmp[0];
    for (uint32_t i = 1; i < 4; i++) {
        ss = fq_add(ss, tmp[i]);
    }
    for (; col < (NN-K); col++) {
        ss = fq_add(ss, d[col]);
	}

    return ss;
}

/// accumulates the inverse of a row
/// \param d
/// \return sum(d) for _ in range(N-K)
static inline
FQ_ELEM row_acc_inv(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < (NN-K); col++) {
        s = fq_add(s, fq_inv(d[col]));
	 }

    return s;
}

/// scalar multiplication of a row
/// /param row[in/out] *= s for _ in range(N-K)
/// /param s
static inline
void row_mul(FQ_ELEM *row, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (NN-K); col++) {
        row[col] = fq_mul_non_ct(s, row[col]);
    }
}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in
/// \param s
static inline
void row_mul2(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (NN-K); col++) {
        out[col] = fq_mul_non_ct(s, in[col]);
    }
}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in
/// \param s
static inline
void row_mul2_ct(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    row_mul2(out, in, s);
}

/// \param out = in1[i]*in2[i] for i in range(N-K)
/// \param in1
/// \param in2
static inline
void row_mul3(FQ_ELEM *out, const FQ_ELEM *in1, const FQ_ELEM *in2) {
    for (uint32_t col = 0; col < (NN-K); col++) {
        out[col] = fq_mul_non_ct(in1[col], in2[col]);
    }
}

/// invert a row
/// \param out[out]: = in[i]**-1 for i in range(N-K)
/// \param in[in]: vector of length N-K
static inline
void row_inv2(FQ_ELEM *out, const FQ_ELEM *in) {
    for (uint32_t col = 0; col < (NN-K); col++) {
        out[col] = fq_inv(in[col]);
    }
}

/// NOTE: not ct.
/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
    for (uint32_t col = 1; col < NN-K; col++) {
        if (in[col] != in[col - 1]) {
            return 0;
        }
    }
    return 1;
}

/// NOTE: not ct.
/// \param in[in]: vector of length N-K
/// \return 0 if no zero was found
///         1 if the row contains at least a single 0
static inline
uint32_t row_contains_zero(const FQ_ELEM *in) {
    for (uint32_t col = 0; col < NN-K; col++) {
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
    for (uint32_t col = 0; col < NN-K; col++) {
        r += (in[col] == 0);
    }
    return r;
}
