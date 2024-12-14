/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
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

#include <string.h> // memcpy, memset
#include <stdlib.h>
#include <stdio.h>

#include "LESS.h"
#include "canonical.h"
#include "seedtree.h"
#include "rng.h"
#include "utils.h"
#include "fips202.h"
#include "keccakf1600.h"
#include "sha3.h"

#if !defined(CATEGORY_1)
/// \param SK
/// \param m
/// \param mlen
/// \param sig
/// \return TODO
size_t LESS_sign(const prikey_t *SK,
                   const char *const m,
                   const uint64_t mlen,
                   sign_t *sig) {
    uint8_t g0_initial_pivot_flags [N];

    /*         Private key expansion        */
    /* expand sequence of seeds for private inverse-monomial matrices */
    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state, SK->compressed_sk, PRIVATE_KEY_SEED_LENGTH_BYTES);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    unsigned char private_monomial_seeds[NUM_KEYPAIRS - 1][PRIVATE_KEY_SEED_LENGTH_BYTES];
    for (uint32_t i = 0; i < NUM_KEYPAIRS - 1; i++) {
        csprng_randombytes(private_monomial_seeds[i],
                           PRIVATE_KEY_SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }

    /*         Ephemeral monomial generation        */
    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    randombytes(ephem_monomials_seed, SEED_LENGTH_BYTES);

    // TODO salt is missing?
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    initialize_csprng(&shake_monomial_state, ephem_monomials_seed, SEED_LENGTH_BYTES);

    /*         Public G_0 expansion                  */
    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, SK->G_0_seed);
    generator_get_pivot_flags (&G0_rref, g0_initial_pivot_flags);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0, &G0_rref);

    monomial_t Q_tilde;
    // TODO: remove the values, as we only are interested in the permutation
    monomial_action_IS_t Q_bar[T];
    normalized_IS_t V_array;

    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    for (uint32_t i = 0; i < T; i++) {
        monomial_mat_seed_expand_rnd(&Q_tilde, &shake_monomial_state);

#if defined(LESS_REUSE_PIVOTS)
        // prepare_digest_input_pivot_reuse(&V_array, &Q_bar[i], &full_G0, &Q_tilde, g0_initial_pivot_flags, SIGN_PIVOT_REUSE_LIMIT);
#else
        prepare_digest_input(&V_array, &Q_bar[i], &full_G0, &Q_tilde);
#endif

        const int t = cf5(&V_array);
        if (t == 0) {
            i -= 1;
        } else {
            LESS_SHA3_INC_ABSORB(&state, (uint8_t *)&V_array, sizeof(normalized_IS_t));
        }
    }

    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *)m, mlen);

    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(sig->digest, &state);

    // (x_0, ..., x_{t-1})
    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string, sig->digest);

    monomial_action_IS_t mono_action;
    for (uint32_t i = 0; i < T; i++) {
        monomial_t Q_to_multiply;
        const int sk_monom_seed_to_expand_idx = fixed_weight_string[i];

        monomial_mat_seed_expand_prikey(&Q_to_multiply,
                                        private_monomial_seeds[sk_monom_seed_to_expand_idx - 1]);
        // NOTE: this function is simplify. We do not need the full monomial matrix.
        monomial_compose_action(&mono_action, &Q_to_multiply, &Q_bar[i]);

        cf_compress_monomial_IS_action(sig->cf_monom_actions[i], &mono_action);
    }

    // always return ok
    return 0;
} /* end LESS_sign */

/// \param PK
/// \param m
/// \param mlen
/// \param sig
/// \return
int LESS_verify(const pubkey_t *const PK,
                const char *const m,
                const uint64_t mlen,
                const sign_t *const sig) {
    uint8_t g_initial_pivot_flags [N];
    uint8_t g_permuted_pivot_flags [N];

    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string, sig->digest);

    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, PK->G_0_seed);

    uint8_t is_pivot_column[N] = {0};
    generator_mat_t tmp_full_G;
    generator_mat_t G_hat;
    normalized_IS_t V_array;
    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    for (uint32_t i = 0; i < T; i++) {
        expand_to_rref(&tmp_full_G, PK->SF_G[fixed_weight_string[i] - 1], g_initial_pivot_flags);
        if (!is_cf_monom_action_valid(sig->cf_monom_actions[i])) {
            return 0;
        }

#if defined(LESS_REUSE_PIVOTS)
        apply_cf_action_to_G_with_pivots(&G_hat,
                                         &tmp_full_G,
                                         sig->cf_monom_actions[employed_monoms],
                                         g_initial_pivot_flags,
                                         g_permuted_pivot_flags);
        const int ret = generator_RREF_pivot_reuse(&G_hat, is_pivot_column,
                                                   g_permuted_pivot_flags,
                                                   VERIFY_PIVOT_REUSE_LIMIT);
#else
        apply_cf_action_to_G(&G_hat, &tmp_full_G, sig->cf_monom_actions[i]);
        const int ret = generator_RREF(&G_hat, is_pivot_column);
#endif
        if(ret != 1) {
            return 0;
        }

        // TODO not correct if more than 1 col is not a pivot column. Somehow merge with the loop just below
        // just copy the non IS
        uint32_t ctr = 0, offset = K;
        for(uint32_t j = 0; j < N-K; j++) {
            if (is_pivot_column[j+K]) {
                ctr += 1;
                offset = K - ctr;
            }

            for (uint32_t t = 0; t < K; t++) {
                V_array.values[t][j] = G_hat.values[t][j + offset];
            }

            offset = K;
        }

        const int r = cf5_nonct(&V_array);
        if (r == 0) { return 0; }
        LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) &V_array.values, sizeof(normalized_IS_t));
    }

    uint8_t recomputed_digest[HASH_DIGEST_LENGTH] = {0};
    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) m, mlen);

    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(recomputed_digest, &state);
    return (verify(recomputed_digest, sig->digest, HASH_DIGEST_LENGTH) == 0);
} /* end LESS_verify */

#endif