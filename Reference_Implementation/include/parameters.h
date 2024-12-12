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
#include <stdint.h>

#define Qm1 (126)
#define MMM 0x204081020408103ull
/* Seed tree max size is computed according to Parameter Generation Script in Utilities folder */

/********************************* Category 1 *********************************/

#if defined(CATEGORY_0)
// NOTE: this is the dev mode. only use it if you know what you are doing
#warning "DEV MODE"

#define SEED_LENGTH_BYTES (16)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#define POSITION_T uint8_t
#define SEED_TREE_LABEL_T uint8_t

#define   N (16)
#define   K (8)
#define   Q (127)
#define NUM_KEYPAIRS (  2)
#define   T (247)
#define   W ( 30)
#define SEED_TREE_MAX_PUBLISHED_BYTES (1472)

#elif defined(CATEGORY_1)
#define SEED_LENGTH_BYTES (16)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#define POSITION_T uint8_t
#define SEED_TREE_LABEL_T uint8_t

#if defined(BALANCED)
#define   N (252)
#define   K (126)
#define   Q (127)
#define NUM_KEYPAIRS (  2)
#define   T (247)
#define   W ( 30)
#define SEED_TREE_MAX_PUBLISHED_BYTES (1472)

#elif defined(INTERMEDIATE)
#define   N (252)
#define   K (126)
#define   Q (127)
#define NUM_KEYPAIRS (  4)
#define   T (244)
#define   W ( 20)
#define SEED_TREE_MAX_PUBLISHED_BYTES (1152)

#elif defined(SHORT_SIG)
#define   N (252)
#define   K (126)
#define   Q (127)
#define NUM_KEYPAIRS (  8)
#define   T (198)
#define   W ( 17)
#define SEED_TREE_MAX_PUBLISHED_BYTES (1056)

#else
#error define optimization corner in parameters.h
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define SEED_LENGTH_BYTES (24)

#if defined(BALANCED)
#define   N (400)
#define   K (200)
#define   Q (127)
#define NUM_KEYPAIRS (  2)
#define   T (759)
#define   W ( 33)
#define SEED_TREE_MAX_PUBLISHED_BYTES (3912)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint16_t
#define POSITION_T uint16_t

#elif defined(SHORT_SIG)
#define   N (400)
#define   K (200)
#define   Q (127)
#define NUM_KEYPAIRS (  3)
#define   T (895)
#define   W ( 26)
#define SEED_TREE_MAX_PUBLISHED_BYTES (3264)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#define POSITION_T uint16_t
#else
#error define optimization corner in parameters.h
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define SEED_LENGTH_BYTES (32)

#if defined(BALANCED)
#define   N (548)
#define   K (274)
#define   Q (127)
#define NUM_KEYPAIRS (  2)
#define   T (1352)
#define   W ( 40)
#define SEED_TREE_MAX_PUBLISHED_BYTES (7168)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#define POSITION_T uint16_t

#elif defined(SHORT_SIG)
#define   N (548)
#define   K (274)
#define   Q (127)
#define NUM_KEYPAIRS (  3)
#define   T (907)
#define   W ( 37)
#define SEED_TREE_MAX_PUBLISHED_BYTES (5600)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t
#define POSITION_T uint16_t
#else
#error define optimization corner in parameters.h
#endif
#else
#error define category for parameters
#endif

#define SIGN_PIVOT_REUSE_LIMIT   (K/2)-32
#define VERIFY_PIVOT_REUSE_LIMIT K 

// #define SIGN_PIVOT_REUSE_LIMIT   0
// #define VERIFY_PIVOT_REUSE_LIMIT 0

/* */
#define K8 ((K+7u)/8u)
#define N8 ((N+7u)/8u)

/// TODO
#define N_K_pad (N-K)

/***************** Derived parameters *****************************************/

#define MAX_PUBLISHED_SEEDS (SEED_TREE_MAX_PUBLISHED_BYTES/SEED_LENGTH_BYTES)

/*length of the output of the cryptographic hash, in bytes */
#define HASH_DIGEST_LENGTH (2*SEED_LENGTH_BYTES)

/* length of the private key seed doubled to avoid multikey attacks */
#define PRIVATE_KEY_SEED_LENGTH_BYTES (2*SEED_LENGTH_BYTES)

#define MASK_Q ((1 << BITS_TO_REPRESENT(Q)) - 1)
#define MASK_N ((1 << BITS_TO_REPRESENT(N)) - 1)


#define IS_REPRESENTABLE_IN_D_BITS(D, N)                \
  (((unsigned long) N>=(1UL << (D-1)) && (unsigned long) N<(1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (15                                     \
                 + IS_REPRESENTABLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(16, N)    \
                 )                                      \
   )

#define NUM_LEAVES_OF_SEED_TREE (1UL << BITS_TO_REPRESENT(T) )
#define NUM_NODES_OF_SEED_TREE (2*NUM_LEAVES_OF_SEED_TREE-1 )

#define RREF_MAT_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*(N-K)*K + 7)/8 + (N + 7)/8)
#define RREF_IS_COLUMNS_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*(N-K)*K + 7)/8)

// #define MONO_ACTION_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*K+7)/8 + (BITS_TO_REPRESENT(N)*K + 7) / 8)

#define LESS_CRYPTO_PUBLICKEYBYTES (NUM_KEYPAIRS*RREF_MAT_PACKEDBYTES)
#define LESS_CRYPTO_SECRETKEYBYTES ((NUM_KEYPAIRS-1)*SEED_LENGTH_BYTES + RREF_MAT_PACKEDBYTES)

// returns the maximum bytes the signature can occupy
#define LESS_CRYPTO_MAX_BYTES (HASH_DIGEST_LENGTH*2 + N8*W + SEED_TREE_MAX_PUBLISHED_BYTES + 1)
#define LESS_CRYPTO_BYTES(NR_LEAVES) (HASH_DIGEST_LENGTH*2 + N8*W + NR_LEAVES*SEED_LENGTH_BYTES + 1)

// TODO
//#define LESS_REUSE_PIVOTS
