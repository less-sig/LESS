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

#define MMM 0x204081020408103ull
/* Seed tree max size is computed according to Parameter Generation Script in Utilities folder */

/***************************** Common Parameters ******************************/
#define Q (127)
#define Qm1 (Q-1)
#define FQ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define POSITION_T uint8_t
#define SEED_TREE_LABEL_T uint8_t

/********************************* Category 1 *********************************/
#if defined(CATEGORY_1)
#define N (252)
#define K (126)
#define SEED_LENGTH_BYTES (16)
#define SIGN_PIVOT_REUSE_LIMIT (25) // Ensures probability of non-CT operation is < 2^-64

#if defined(BALANCED)
// Whenever the SEED_TREE macro is absent for a parameter set, no trees will be used
#define SEED_TREE
#define NUM_KEYPAIRS (2)
#define T (192)
#define W (36)
#define TREE_OFFSETS {0, 0, 0, 0, 0, 0, 0, 0, 128}
#define TREE_NODES_PER_LEVEL {1, 2, 4, 8, 16, 32, 64, 128, 128}
#define TREE_LEAVES_PER_LEVEL {0, 0, 0, 0, 0, 0, 0, 64, 128}
#define TREE_SUBROOTS 2
#define TREE_LEAVES_START_INDICES {255, 191}
#define TREE_CONSECUTIVE_LEAVES {128, 64}
#define MAX_PUBLISHED_SEEDS 388

#elif defined(INTERMEDIATE)
#define NUM_KEYPAIRS (4)
#define T (68)
#define W (42)

#elif defined(SHORT_SIG)
#define NUM_KEYPAIRS (8)
#define T (45)
#define W (34)

#else
#error define parameters in parameters.h
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define N (400)
#define K (200)
#define SEED_LENGTH_BYTES (24)
#define SIGN_PIVOT_REUSE_LIMIT (51) // Ensures probability of non-CT operation is < 2^-64
#define POSITION_T uint16_t

#if defined(BALANCED)
#define NUM_KEYPAIRS (2)
#define T (220)
#define W (68)

#elif defined(SHORT_SIG)
#define NUM_KEYPAIRS (3)
#define T (102)
#define W (61)
#else
#error define optimization corner in parameters.h
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define N (548)
#define K (274)
#define SEED_LENGTH_BYTES (32)
#define SIGN_PIVOT_REUSE_LIMIT (79) // Ensures probability of non-CT operation is < 2^-64
#define POSITION_T uint16_t

#if defined(BALANCED)
#define NUM_KEYPAIRS (2)
#define T (345)
#define W (75)

#elif defined(SHORT_SIG)
#define NUM_KEYPAIRS (3)
#define T (137)
#define W (79)
#else
#error define optimization corner in parameters.h
#endif

#else
#error define category for parameters
#endif


#define VERIFY_PIVOT_REUSE_LIMIT K

// #define SIGN_PIVOT_REUSE_LIMIT   0
// #define VERIFY_PIVOT_REUSE_LIMIT 0

/* */
#define K8 ((K+7u)/8u)
#define N8 ((N+7u)/8u)

// TODO
#define N_K_pad (N-K)

/***************** Derived parameters *****************************************/



/*length of the output of the cryptographic hash, in bytes */
#define HASH_DIGEST_LENGTH (2*SEED_LENGTH_BYTES)
#define SALT_LENGTH_BYTES HASH_DIGEST_LENGTH
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
#define LOG2(L) ( (BITS_TO_REPRESENT(L) > BITS_TO_REPRESENT(L-1)) ? (BITS_TO_REPRESENT(L-1)) : (BITS_TO_REPRESENT(L)) )

#define NUM_LEAVES_SEED_TREE (T)
#define NUM_NODES_SEED_TREE (2*NUM_LEAVES_SEED_TREE-1)

#define RREF_MAT_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*(N-K)*K + 7)/8 + (N + 7)/8)
#define RREF_IS_COLUMNS_PACKEDBYTES ((BITS_TO_REPRESENT(Q)*(N-K)*K + 7)/8)

#define LESS_CRYPTO_PUBLICKEYBYTES (NUM_KEYPAIRS*RREF_MAT_PACKEDBYTES)
#define LESS_CRYPTO_SECRETKEYBYTES ((NUM_KEYPAIRS-1)*SEED_LENGTH_BYTES + RREF_MAT_PACKEDBYTES)

#if defined(SEED_TREE)
// returns the maximum bytes the signature can occupy
#define SEED_TREE_MAX_PUBLISHED_BYTES (MAX_PUBLISHED_SEEDS*SEED_LENGTH_BYTES)
#define LESS_CRYPTO_MAX_BYTES        (HASH_DIGEST_LENGTH*2 + N8*W + SEED_TREE_MAX_PUBLISHED_BYTES + 1)
#else
// returns the maximum bytes the signature can occupy
#define LESS_CRYPTO_BYTES     (HASH_DIGEST_LENGTH + (N8*W) + ((W-T)*SEED_LENGTH_BYTES))
#endif

// TODO doc
#define LESS_REUSE_PIVOTS_VY
#define LESS_REUSE_PIVOTS_SG
