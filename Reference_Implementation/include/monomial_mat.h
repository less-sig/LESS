/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
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
#include "fq_arith.h"

/* Structure representing a monomial matrix */
/* consider the example n=3
 *     [ 2  0  0 ]   [             ]
 * M = [ 0  0  1 ] = [ m_0 m_1 m_2 ]
 *     [ 0  3  0 ]   [             ]
 *
 * Computing GM shuffles the columns of G and rescales them. Assume G is:
 *     [             ]            [                   ]
 * G = [ g_0 g_1 g_2 ]  then GM = [ 2*g_0 3*g_2 1*g_1 ]
 *     [             ]            [                   ]
 *
 * M is stored in compact form as two arrays, the first one stores the scaling
 * coefficients of the columns (i.e., [2 3 1] in the example), while the second
 * stores, for each element, the position of the column placed at its index
 * after the permutation.
 * In the example we have g_0 -> 0, g_1 -> 2, g_2 -> 1 hence we obtain [0 2 1]
 *
 */

typedef struct {
   /* coefficients listed in order of appearance column-wise */
   FQ_ELEM coefficients[N];
   /* considering the product GQ, permutation[...] stores into the cell with
    * index 0, the position of the DESTINATION of column 0 in G after the
    * computation of GQ.
    */
   POSITION_T permutation[N];
} monomial_t;

typedef struct {
    // NOTE: as we are now computing canonical forms, we do not need
    // a full monomial matrix anymore. But only the permutation. So
    // the name of the struct is maybe a little bit wrong.
    POSITION_T permutation[K];
} monomial_action_IS_t;

typedef struct {
   unsigned char value[SEED_LENGTH_BYTES];
} monomial_seed_t;

void yt_shuffle_state_limit(SHAKE_STATE_STRUCT *shake_monomial_state,
                            POSITION_T *permutation,
                            const uint32_t n);
void yt_shuffle_state(SHAKE_STATE_STRUCT *shake_monomial_state,
                      POSITION_T permutation[N]);
void yt_shuffle(POSITION_T permutation[N]);

/* multiplies two monomial matrices */
void monomial_mat_mul(monomial_t *res,
                      const monomial_t *const A,
                      const monomial_t *const B);

/* computes the inverse of the monomial matrix */
void monomial_inv(monomial_t *res,
                  const monomial_t *const to_invert);

/* expands a monomial matrix, given a PRNG seed and a salt (used for ephemeral
 * monomial matrices */
void monomial_sample_salt(monomial_t *res,
                          const unsigned char seed[SEED_LENGTH_BYTES],
                          const unsigned char salt[HASH_DIGEST_LENGTH],
                          const uint16_t round_index);

/* expands a monomial matrix, given a double length PRNG seed (used to prevent
 * multikey attacks) */
void monomial_sample_prikey(monomial_t *res,
                            const unsigned char seed[PRIVATE_KEY_SEED_LENGTH_BYTES]);

///
void monomial_mat_seed_expand_rnd(monomial_t *res,
                                  const unsigned char seed[SEED_LENGTH_BYTES],
                                  const uint16_t round_index);

/* composes a compressed action on an IS with the action of a monomial
 * matrix */
void monomial_compose_action(monomial_action_IS_t * out, 
                             const monomial_t * to_compose, 
                             const monomial_action_IS_t * in);

void CosetRep(uint8_t *b, const monomial_action_IS_t *Q_star);

/* Decompress byte array to MonomialAction object */
void expand_to_monom_action(monomial_action_IS_t *mono,
                            const uint8_t *compressed);

/* Validate MonomialAction object */
int CheckCanonicalAction(const uint8_t* const mono);
