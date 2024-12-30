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

void LESS_keygen(prikey_t *SK,
                 pubkey_t *PK) {
    /* generating private key from a single seed */
    randombytes(SK->compressed_sk, PRIVATE_KEY_SEED_LENGTH_BYTES);
    /* expanding it onto private (inverse) monomial seeds */
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

    /* Generating public code G_0 */
    randombytes(SK->G_0_seed, SEED_LENGTH_BYTES);
    memcpy(PK->G_0_seed, SK->G_0_seed, SEED_LENGTH_BYTES);

    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, SK->G_0_seed);

    generator_mat_t tmp_full_G;
    generator_rref_expand(&tmp_full_G, &G0_rref);
    uint8_t is_pivot_column[N];

    /* note that the first "keypair" is just the public generator G_0, stored
     * as a seed and the identity matrix (not stored) */
    for (uint32_t i = 0; i < NUM_KEYPAIRS - 1; i++) {
        /* expand inverse monomial from seed */
        monomial_t private_Q;
        monomial_t private_Q_inv;
        monomial_mat_seed_expand_prikey(&private_Q_inv, private_monomial_seeds[i]);
        monomial_mat_inv(&private_Q, &private_Q_inv);

        generator_mat_t result_G;
        generator_monomial_mul(&result_G,
                               &tmp_full_G,
                               &private_Q);
        memset(is_pivot_column, 0, sizeof(is_pivot_column));
        generator_RREF(&result_G, is_pivot_column);
        /* note that the result is stored at i-1 as the first
         * public key element is just a seed */
        compress_rref(PK->SF_G[i],
                      &result_G,
                      is_pivot_column);
    }
} /* end LESS_keygen */

// TODO, quick and dirty hack
#ifdef SEED_TREE

/// returns the number of opened seeds in the tree.
/// \param SK
/// \param m
/// \param mlen
/// \param sig
/// \return
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
    randombytes(sig->tree_salt, HASH_DIGEST_LENGTH);

    /* create the prng for the "blinding" monomials for the canonical form computation */
    uint8_t cf_seed[SEED_LENGTH_BYTES];
    randombytes(cf_seed, SEED_LENGTH_BYTES);
    SHAKE_STATE_STRUCT cf_shake_state;
    initialize_csprng(&cf_shake_state, cf_seed, SEED_LENGTH_BYTES);

    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    generate_seed_tree_from_root(seed_tree, ephem_monomials_seed, sig->tree_salt);
    unsigned char *ephem_monomial_seeds = seed_tree +
                                          SEED_LENGTH_BYTES * (NUM_LEAVES_OF_SEED_TREE - 1);

    /*         Public G_0 expansion                  */
    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, SK->G_0_seed);
    generator_get_pivot_flags (&G0_rref, g0_initial_pivot_flags);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0, &G0_rref);

    monomial_t Q_tilde;
    monomial_action_IS_t Q_bar[T];
    normalized_IS_t V_array;

    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    for (uint32_t i = 0; i < T; i++) {
        monomial_mat_seed_expand_salt_rnd(&Q_tilde,
                                          ephem_monomial_seeds + i * SEED_LENGTH_BYTES,
                                          sig->tree_salt,
                                          i);

#if defined(LESS_REUSE_PIVOTS_SG)
            // TODO half of these operations within this function can be removed
            prepare_digest_input_pivot_reuse(&V_array,
                                             &Q_bar[i],
                                             &full_G0,
                                             &Q_tilde,
                                             g0_initial_pivot_flags,
                                             SIGN_PIVOT_REUSE_LIMIT);
#else
        prepare_digest_input(&V_array, &Q_bar[i], &full_G0, &Q_tilde);
#endif
        // blind(&V_array, &cf_shake_state);
        normalized_IS_t V_array2;
        normalized_copy(&V_array2, &V_array);
        const int t = cf5_nonct(&V_array2);
        if (t == 0) {
            *(ephem_monomial_seeds + i*SEED_LENGTH_BYTES) += 1;
            i -= 1;
            // TODO, what happens in this case?
            normalized_pretty_print(&V_array);
            printf("\n");
            printf("\n");
            normalized_pretty_print(&V_array2);
            printf("\n");
            printf("cf5 failed\n");
        } else {
            LESS_SHA3_INC_ABSORB(&state, (uint8_t *)&V_array2, sizeof(normalized_IS_t));
        }
    }

    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *)m, mlen);
    LESS_SHA3_INC_ABSORB(&state, sig->tree_salt, HASH_DIGEST_LENGTH);

    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(sig->digest, &state);

    // (x_0, ..., x_{t-1})
    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string, sig->digest);

    uint8_t indices_to_publish[T];
    for (uint32_t i = 0; i < T; i++) {
        indices_to_publish[i] = !!(fixed_weight_string[i]);
    }

    int emitted_monoms = 0;
    memset(&sig->seed_storage, 0, SEED_TREE_MAX_PUBLISHED_BYTES);

    const uint32_t num_seeds_published =
            seed_tree_path(seed_tree,
                           indices_to_publish,
                           (unsigned char *) &sig->seed_storage);

    monomial_action_IS_t mono_action;
    for (uint32_t i = 0; i < T; i++) {
        monomial_t Q_to_multiply;
        if (fixed_weight_string[i] != 0) {
            const int sk_monom_seed_to_expand_idx = fixed_weight_string[i];

            monomial_mat_seed_expand_prikey(&Q_to_multiply,
                                            private_monomial_seeds[sk_monom_seed_to_expand_idx - 1]);
            // NOTE: this function is simplify. We do not need the full monomial matrix.
            monomial_compose_action(&mono_action, &Q_to_multiply, &Q_bar[i]);

            cf_compress_monomial_IS_action(sig->cf_monom_actions[emitted_monoms], &mono_action);
            emitted_monoms++;
        }
    }

    // TODO this needs to be changed. As described in the TODO in overleaf, currently
    // we need to keep track of the opened commitments.
    sig->seed_storage[num_seeds_published*SEED_LENGTH_BYTES] = num_seeds_published;
    return num_seeds_published;
} /* end LESS_sign */


int LESS_verify(const pubkey_t *const PK,
                const char *const m,
                const uint64_t mlen,
                const sign_t *const sig) {

    uint8_t fixed_weight_string[T] = {0};
    uint8_t g_initial_pivot_flags [N];
    uint8_t g_permuted_pivot_flags [N];
    expand_digest_to_fixed_weight(fixed_weight_string, sig->digest);

    uint8_t published_seed_indexes[T];
    for (uint32_t i = 0; i < T; i++) {
        published_seed_indexes[i] = !!(fixed_weight_string[i]);
    }

    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    rebuild_seed_tree_leaves(seed_tree, published_seed_indexes,
                             (unsigned char *) &sig->seed_storage, sig->tree_salt);

    unsigned char *ephem_monomial_seeds = seed_tree +
                                          SEED_LENGTH_BYTES * (NUM_LEAVES_OF_SEED_TREE - 1);

    int employed_monoms = 0;

    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, PK->G_0_seed);

    uint8_t is_pivot_column[N] = {0};
    generator_mat_t tmp_full_G;
    generator_mat_t G_hat;
    monomial_action_IS_t Q_to_discard;
    normalized_IS_t V_array;
    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    for (uint32_t i = 0; i < T; i++) {
        if (fixed_weight_string[i] == 0) {
            generator_get_pivot_flags(&G0_rref, g_initial_pivot_flags);
          
            // TODO mov this out of the loop and keep `tmp_full_G` constant
            generator_rref_expand(&tmp_full_G, &G0_rref);
            monomial_t Q_to_multiply;
            monomial_mat_seed_expand_salt_rnd(&Q_to_multiply,
                                              ephem_monomial_seeds + i * SEED_LENGTH_BYTES,
                                              sig->tree_salt,
                                              i);
#if defined(LESS_REUSE_PIVOTS_VY)
            // TODO half of these operations can be optimized away
            prepare_digest_input_pivot_reuse(&V_array,
                                             &Q_to_discard,
                                             &tmp_full_G,
                                             &Q_to_multiply,
                                             g_initial_pivot_flags,
                                             VERIFY_PIVOT_REUSE_LIMIT);
#else
            prepare_digest_input(&V_array,
                                 &Q_to_discard,
                                 &tmp_full_G,
                                 &Q_to_multiply);
#endif

            const int r = cf5_nonct(&V_array);
            if (r == 0) {
                // NOTE: we just silently reject the signature, if we do not
                // have a valid CF input. This should only happen with an
                // negl. probability.
                return 0;
            }
            LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) &V_array, sizeof(normalized_IS_t));
        } else {
            expand_to_rref(&tmp_full_G, PK->SF_G[fixed_weight_string[i] - 1], g_initial_pivot_flags);
            if (!is_cf_monom_action_valid(sig->cf_monom_actions[employed_monoms])) {
                return 0;
            }

#if defined(LESS_REUSE_PIVOTS_VY)
            apply_cf_action_to_G_with_pivots(&G_hat,
                                             &tmp_full_G,
                                             sig->cf_monom_actions[employed_monoms],
                                             g_initial_pivot_flags,
                                             g_permuted_pivot_flags);
            const int ret = generator_RREF_pivot_reuse(&G_hat, is_pivot_column,
                                                       g_permuted_pivot_flags,
                                                       VERIFY_PIVOT_REUSE_LIMIT);
#else
            apply_cf_action_to_G(&G_hat, &tmp_full_G, sig->cf_monom_actions[employed_monoms]);
            const int ret = generator_RREF(&G_hat, is_pivot_column);
#endif
            if(ret != 1) {
                return 0;
            }

            // TODO not CT, not correct if more than 1 col is not a pivot column. Somehow merge with the loop just below
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
            if (r == 0) {
                // NOTE: we just silently reject the signature, if we do not
                // have a valid CF input. This should only happen with an
                // negl. probability.
                return 0;
            }
            LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) &V_array.values, sizeof(normalized_IS_t));
            employed_monoms++;
        }
    }

    uint8_t recomputed_digest[HASH_DIGEST_LENGTH] = {0};
    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) m, mlen);
    LESS_SHA3_INC_ABSORB(&state, sig->tree_salt, HASH_DIGEST_LENGTH);

    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(recomputed_digest, &state);
    return (verify(recomputed_digest, sig->digest, HASH_DIGEST_LENGTH) == 0);
} /* end LESS_verify */

#endif
