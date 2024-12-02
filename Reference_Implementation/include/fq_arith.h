/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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

#define NUM_BITS_Q (BITS_TO_REPRESENT(Q))

/*
 * assertion _only_ in debug builds (CMAKE_BUILD_TYPE=Debug)
*/
#if !defined(NDEBUG) || defined(DEBUG)
#include <assert.h>
#define ASSERT(X) assert(X);
#else
#define ASSERT(X)  {(void) (X);}
#endif

#define DEF_RAND_STATE(FUNC_NAME, EL_T, MINV, MAXV) \
static inline void FUNC_NAME(SHAKE_STATE_STRUCT *shake_monomial_state, EL_T *buffer, size_t num_elements) { \
   typedef uint64_t WORD_T; \
   static const EL_T MIN_VALUE = (MINV);\
   static const EL_T MAX_VALUE = (MAXV); \
   ASSERT(MIN_VALUE <= MAX_VALUE); \
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
   ASSERT(MIN_VALUE <= MAX_VALUE); \
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


/* GCC actually inlines and vectorizes Barrett's reduction already.
 * Backup implementation for less aggressive compilers follows */
#if 0
#define BARRETT_MU  (((uint32_t)1<<(2*NUM_BITS_Q))/Q)
#define BARRETT_MASK ( ((FQ_DOUBLEPREC)1 << (NUM_BITS_Q+3))-1 )

static inline
FQ_ELEM fq_red(FQ_DOUBLEPREC a)
{
   FQ_DOUBLEPREC q_1, q_2, q_3;
   q_1 = a >> (NUM_BITS_Q);
   q_2 = q_1 * BARRETT_MU;
   q_3 = q_2 >> (NUM_BITS_Q);
   FQ_DOUBLEPREC r_1;
   r_1 = (a & BARRETT_MASK) - ( (q_3*Q) & BARRETT_MASK);
   r_1 = r_1 & BARRETT_MASK;
   FQ_ELEM r_2;
   FQ_DOUBLEPREC need_to_red;
   need_to_red = r_1 >= Q;
   r_1 = r_1-Q*need_to_red; // not needed for 127
   need_to_red = r_1 >= Q;
   r_2 = r_1-Q*need_to_red;
   return r_2;
}
#endif

#if 0
/* Fast Mersenne prime reduction is actually slower than Barrett's */
static inline
FQ_ELEM fq_red(FQ_DOUBLEPREC x)
{
   while (x>=Q) {
      x = ((FQ_DOUBLEPREC) 0x7f & x) + (x>>7);
   }
   return x;
}
#endif

static inline
FQ_ELEM fq_red(FQ_DOUBLEPREC x) {
   return ((FQ_DOUBLEPREC) Q+x) % (FQ_DOUBLEPREC) Q;
}

/// super generic a \in [0, 2**16)
//static inline
//uint8_t fq_red(const uint16_t a) {
//   const uint64_t lowbits = MMM * a;
//   const uint32_t highbits = ((__uint128_t)lowbits * Q) >> 64;
//   return highbits - (Qm1 & (a >> 31));
//}

static inline
FQ_ELEM fq_mul(const FQ_ELEM x, const FQ_ELEM y) {
   return ((FQ_DOUBLEPREC)x * (FQ_DOUBLEPREC)y) % Q;
}

//static inline
//FQ_ELEM fq_mul(const FQ_ELEM x,
//               const FQ_ELEM y) {
//   FQ_DOUBLEPREC z = x*y;
//   return fq_red(z);
//}

//static inline
//FQ_ELEM fq_mul(FQ_ELEM a, FQ_ELEM b){
//   FQ_DOUBLEPREC lo, hi;
//   lo = a * b;
//   hi = lo >> 7;
//   lo += hi;
//   hi <<= 7;
//   return (lo - hi) & Q;
//}

// // Barrett reduction for uint8_t with prime Q = 127
// static inline
// FQ_ELEM fq_red(FQ_ELEM a) {
//    FQ_ELEM t;
//    t = a >> 7;
//    t &= Q;
//    a += t;
//    a &= Q;
//    return a;
// }

static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
      return (x + y) % Q;
}

/// NOTE: maybe dont use it for sensetive data
static const uint8_t fq_inv_table[128] __attribute__((aligned(64))) = {
   0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};


/* Fermat's method for inversion employing r-t-l square and multiply,
 * unrolled for actual parameters */
static inline
FQ_ELEM fq_inv(const FQ_ELEM x) {
    const FQ_ELEM t = x%Q; // TODO, the problem is, that the pseudo reduction, does not fully reduce
   return fq_inv_table[t];
} /* end fq_inv */


static inline
FQ_ELEM fq_pow(FQ_ELEM x, FQ_ELEM exp) {
   FQ_DOUBLEPREC xlift;
   xlift = x;
   FQ_DOUBLEPREC accum = 1;
   /* No need for square and mult always, Q-2 is public*/
   while(exp) {
      if(exp & 1) {
         accum = fq_red(accum*xlift);
      }
      xlift = fq_red(xlift*xlift);
      exp >>= 1;
   }
   return fq_red(accum);
} /* end fq_pow */

/* Sampling functions from the global TRNG state */

DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)

DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

/* Sampling functions from the taking the PRNG state as a parameter*/
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)

DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)

/// NOTE: these functions are outsourced to this file, to make the
/// optimizied implementation as easy as possible.
/// accumulates a row
/// @param d
/// @return sum(d) for _ in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < (N-K); col++) {
        s = fq_add(s, d[col]);
	 }

    return s;
}

/// accumulates the inverse of a row
/// @param d
/// @return sum(d) for _ in range(N-K)
static inline
FQ_ELEM row_acc_inv(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < (N-K); col++) {
        s = fq_add(s, fq_inv(d[col]));
	 }

    return s;
}

/// scalar multiplication of a row
/// @param row[in/out] *= s for _ in range(N-K)
/// @param s
static inline
void row_mul(FQ_ELEM *row, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (N-K); col++) {
        row[col] = fq_mul(s, row[col]);
    }
}

/// scalar multiplication of a row
/// @param out = s*in[i] for i in range(N-K)
/// @param in
/// @param s
static inline
void row_mul2(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_mul(s, in[col]);
    }
}

///
/// @param out = in1[i]*in2[i] for i in range(N-K)
/// @param in1
/// @param in2
static inline
void row_mul3(FQ_ELEM *out, const FQ_ELEM *in1, const FQ_ELEM *in2) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_mul(in1[col], in2[col]);
    }
}

/// invert a row
/// @param out = in[i]**-1 for i in range(N-K)
/// @param in
static inline
void row_inv2(FQ_ELEM *out, const FQ_ELEM *in) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_inv(in[col]);
    }
}

/// @param in
/// @return 1 if all elements are the same
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
