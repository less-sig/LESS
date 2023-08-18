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
      for (int i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
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
      for (int i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
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
FQ_ELEM fq_red(FQ_DOUBLEPREC x)
{
   return ((FQ_DOUBLEPREC) Q+x) % (FQ_DOUBLEPREC) Q;
}

/* Fermat's method for inversion employing r-t-l square and multiply,
 * unrolled for actual parameters */
static inline
FQ_ELEM fq_inv(FQ_ELEM x)
{

   FQ_DOUBLEPREC xlift;
   xlift = x;
   FQ_DOUBLEPREC accum = 1;
   /* No need for square and mult always, Q-2 is public*/
   uint32_t exp = Q-2;
   while(exp) {
      if(exp & 1) {
         accum = fq_red(accum*xlift);
      }
      xlift = fq_red(xlift*xlift);
      exp >>= 1;
   }
   return fq_red(accum);
} /* end fq_inv */

/* Sampling functions from the global TRNG state */

DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)

DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

/* Sampling functions from the taking the PRNG state as a parameter*/
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)

DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)
