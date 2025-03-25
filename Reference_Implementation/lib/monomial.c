/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Floyd Zweydinger <zweydfg+github@rub.de>
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
#include "monomial_mat.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <string.h>

///
#define POS_BITS BITS_TO_REPRESENT(N-1)

///
#define POS_MASK (((POSITION_T) 1 << POS_BITS) - 1)

/// applies a random permutation between [0, n-1] on the
/// input permutation
/// \param shake_monomial_state[in/out]:
/// \param permutation[in/out]: random permutation. Must be initialized with
///         [0,....,n-1]
/// \param n <= N: number of elements in the permutation
void yt_shuffle_state_limit(SHAKE_STATE_STRUCT *shake_monomial_state,
                            POSITION_T *permutation,
                            const uint32_t n) {
    uint32_t rand_u32[N] = {0};
    POSITION_T tmp;

    csprng_randombytes((unsigned char *) &rand_u32, sizeof(uint32_t)*n, shake_monomial_state);
    for (size_t i = 0; i < n - 1; ++i) {
        rand_u32[i] = i + rand_u32[i] % (n - i);
    }

    for (size_t i = 0; i < n - 1; ++i) {
        tmp = permutation[i];
        permutation[i] = permutation[rand_u32[i]];
        permutation[rand_u32[i]] = tmp;
    }
}

/// \param shake_monomial_state[in/out]:
/// \param permutation[in/out]: random permutation. Must be initialized with
///         [0,....,n-1]
void yt_shuffle_state(SHAKE_STATE_STRUCT *shake_monomial_state, POSITION_T permutation[N]) {
   uint64_t rand_u64;
   POSITION_T tmp;
   POSITION_T x;
   int c;

   csprng_randombytes((unsigned char *) &rand_u64,
                             sizeof(rand_u64),
                             shake_monomial_state);
   c = 0;

   for (int i = 0; i < N; i++) {
      do {
         if (c == (64/POS_BITS)-1) {
            csprng_randombytes((unsigned char *) &rand_u64,
                                sizeof(rand_u64),
                                shake_monomial_state);
            c = 0;
         }
         x = rand_u64 & (POS_MASK);
         rand_u64 = rand_u64 >> POS_BITS;
         c = c + 1;
      } while (x >= N);

      tmp = permutation[i];
      permutation[i] = permutation[x];
      permutation[x] = tmp;
   } 
}
/* FY shuffle on the permutation, sampling from the global TRNG state */
void yt_shuffle(POSITION_T permutation[N]) {
    yt_shuffle_state(&platform_csprng_state, permutation);
}

/* expands a monomial matrix, given a PRNG seed and a salt (used for ephemeral
 * monomial matrices */
void monomial_sample_salt(monomial_t *res,
                          const unsigned char seed[SEED_LENGTH_BYTES],
                          const unsigned char salt[HASH_DIGEST_LENGTH],
                          const uint16_t round_index) {
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    const int shake_buffer_len = SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH + sizeof(uint16_t);
    uint8_t shake_input_buffer[shake_buffer_len];
    memcpy(shake_input_buffer, seed, SEED_LENGTH_BYTES);
    memcpy(shake_input_buffer + SEED_LENGTH_BYTES, salt, HASH_DIGEST_LENGTH);
    memcpy(shake_input_buffer + SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH, &round_index, sizeof(uint16_t));

    initialize_csprng(&shake_monomial_state, shake_input_buffer, shake_buffer_len);
    fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
    for (uint32_t i = 0; i < N; i++) {
        res->permutation[i] = i;
    }

    /* FY shuffle on the permutation */
    yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */

/// expands a monomial matrix, given a double length PRNG seed (used to prevent
/// multikey attacks)
/// \param res[out]: the randomly sampled monomial matrix
/// \param seed[in]: the seed
void monomial_sample_prikey(monomial_t *res,
                            const unsigned char seed[PRIVATE_KEY_SEED_LENGTH_BYTES]) {
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    initialize_csprng(&shake_monomial_state, seed, PRIVATE_KEY_SEED_LENGTH_BYTES);
    fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
    for (uint32_t i = 0; i < N; i++) {
        res->permutation[i] = i;
    }
    /* FY shuffle on the permutation */
    yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */

/// \param res[out]: = to_invert**-1
/// \param to_invert[in]:
void monomial_inv(monomial_t *res,
                  const monomial_t *const to_invert) {
    for(uint32_t i = 0; i < N; i++) {
        res->permutation[to_invert->permutation[i]] = i;
        res->coefficients[to_invert->permutation[i]] =
            fq_inv(to_invert->coefficients[i]);
    }
} /* end monomial_inv */

/* composes a compactly stored action of a monomial on an IS with a regular
 * monomial.
 * NOTE: Only the permutation is computed, as this is the only thing we need
 * since the adaption of canonical forms.
 */
void monomial_compose_action(monomial_action_IS_t *out,
                             const monomial_t *Q_in,
                             const monomial_action_IS_t *in) {
    /* to compose with monomial_action_IS_t, reverse the convention
     * for Q storage: store in permutation[i] the idx of the source column landing
     * as the i-th after the GQ product, and in coefficients[i] the coefficient
     * by which the column is multiplied upon landing */
    monomial_t reverse_Q;
    for (uint32_t i = 0; i < N; i++) {
        reverse_Q.permutation[Q_in->permutation[i]] = i;
    }
    /* compose actions out = Q_in*in */
    for (uint32_t i = 0; i < K; i++) {
        out->permutation[i] = reverse_Q.permutation[in->permutation[i]];
    }
}

/// cf type5 compression
/// \param b[out]: N bits in which K bits will be sed
/// \param Q_star[in]: canonical action matrix
void CosetRep(uint8_t *b,
              const monomial_action_IS_t *Q_star) {
    memset(b, 0, N8);
    for (uint32_t i = 0; i < K; i++) {
        const uint32_t limb = (Q_star->permutation[i])/8;
        const uint32_t pos  = (Q_star->permutation[i])%8;
        b[limb] ^= 1u << pos;
    }
}

/// checks if the given (N+7/8) bytes are a valid
/// canonical form action.
/// \param b[in]: compressed canonical form output from
///         `CosetRep`
/// \return true: if the weight is K
///         false: if the weight is not K
int CheckCanonicalAction(const uint8_t* const b) {
    uint32_t w = 0;
    for (uint32_t i = 0; i < N8; i++) {
        w += (uint32_t)__builtin_popcount(b[i]);
    }

    return w == K;
}
