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
#include <string.h> // memcpy, memset
#include "LESS.h"
#include "canonical.h"
#include "seedtree.h"
#include "rng.h"
#include "utils.h"
#include "fips202.h"
#include "sha3.h"

void LESS_keygen(prikey_t *SK,
                 pubkey_t *PK) {
    /* generating private key from a single seed */
    randombytes(SK->compressed_sk, PRIVATE_KEY_SEED_LENGTH_BYTES);

    /* expanding it onto private seeds */
    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state, SK->compressed_sk, PRIVATE_KEY_SEED_LENGTH_BYTES);

    /* Generating public code G_0 */
    csprng_randombytes(PK->G_0_seed, SEED_LENGTH_BYTES, &sk_shake_state);

    rref_generator_mat_t G0_rref;
    generator_sample(&G0_rref, PK->G_0_seed);

    generator_mat_t tmp_full_G;
    generator_rref_expand(&tmp_full_G, &G0_rref);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    unsigned char private_monomial_seeds[NUM_KEYPAIRS - 1][PRIVATE_KEY_SEED_LENGTH_BYTES];
    for (uint32_t i = 0; i < NUM_KEYPAIRS - 1; i++) {
        csprng_randombytes(private_monomial_seeds[i],
                           PRIVATE_KEY_SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }

    /* note that the first "keypair" is just the public generator G_0, stored
     * as a seed and the identity matrix (not stored) */
    for (uint32_t i = 0; i < NUM_KEYPAIRS - 1; i++) {
        uint8_t is_pivot_column[N_pad] = {0};
        /* expand inverse monomial from seed */
        monomial_t private_Q;
        monomial_t private_Q_inv;
        monomial_sample_prikey(&private_Q_inv, private_monomial_seeds[i]);
        monomial_inv(&private_Q, &private_Q_inv);

        generator_mat_t result_G = {0};
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

/// returns the number of opened seeds in the tree.
/// \param SK[in]: secret key
/// \param m[in]: message to sign
/// \param mlen[in]: length of the message to sign in bytes
/// \param sig[out]: signature
/// \return: x: number of leaves opened by the algorithm
size_t LESS_sign(const prikey_t *SK,
                 const char *const m,
                 const uint64_t mlen,
                 sign_t *sig) {
    uint8_t g0_initial_pivot_flags [N_pad];
    uint8_t is_pivot_column[N_pad];

    /*         Private key expansion        */
    /* expand sequence of seeds for private inverse-monomial matrices */
    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state, SK->compressed_sk, PRIVATE_KEY_SEED_LENGTH_BYTES);

    /* Generating seed for public code G_0 */
    unsigned char G_0_seed[SEED_LENGTH_BYTES];
    csprng_randombytes(G_0_seed, SEED_LENGTH_BYTES, &sk_shake_state);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    unsigned char private_monomial_seeds[NUM_KEYPAIRS - 1][PRIVATE_KEY_SEED_LENGTH_BYTES];
    for (uint32_t i = 0; i < NUM_KEYPAIRS - 1; i++) {
        csprng_randombytes(private_monomial_seeds[i],
                           PRIVATE_KEY_SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }

    // generate the salt from a TRNG
    randombytes(sig->salt, HASH_DIGEST_LENGTH);

    /*         Ephemeral monomial generation        */
    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    csprng_randombytes(ephem_monomials_seed,
                       SEED_LENGTH_BYTES,
                       &sk_shake_state);

    /* create the prng for the "blinding" monomials for the canonical form computation */
    uint8_t cf_seed[SEED_LENGTH_BYTES];
    csprng_randombytes(cf_seed,
                       SEED_LENGTH_BYTES,
                       &sk_shake_state);
    SHAKE_STATE_STRUCT cf_shake_state;
    initialize_csprng(&cf_shake_state, cf_seed, SEED_LENGTH_BYTES);

    unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    BuildGGM(seed_tree, ephem_monomials_seed, sig->salt);

    unsigned char linearized_rounds_seeds[T*SEED_LENGTH_BYTES] = {0};
    seed_leaves(linearized_rounds_seeds,seed_tree);

    /*         Public G_0 expansion                  */
    rref_generator_mat_t G0_rref;
    generator_sample(&G0_rref, G_0_seed);
    generator_get_pivot_flags (&G0_rref, g0_initial_pivot_flags);
    generator_mat_t full_G0, G0;
    generator_rref_expand(&full_G0, &G0_rref);

    monomial_t mu_tilde;
    monomial_action_IS_t pi_tilde[T];
    normalized_IS_t A_i = {0};

    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    for (uint32_t i = 0; i < T; i++) {
        monomial_sample_salt(&mu_tilde,
                             linearized_rounds_seeds + i * SEED_LENGTH_BYTES,
                             sig->salt,
                             i);
        generator_monomial_mul(&G0, &full_G0, &mu_tilde);
        memset(is_pivot_column, 0, N_pad);
#if defined(LESS_REUSE_PIVOTS_SG)
        uint8_t permuted_pivot_flags[N_pad];
        for (uint32_t t = 0; t < N; t++) {
            permuted_pivot_flags[mu_tilde.permutation[t]] = g0_initial_pivot_flags[t];
        }
        if (generator_RREF_pivot_reuse(&G0, is_pivot_column, permuted_pivot_flags, SIGN_PIVOT_REUSE_LIMIT) == 0) {
            return 0;
        }
#else
        if (generator_RREF(&G0, is_pivot_column) == 0) {
            return 0;
        }
#endif
        // just copy the non IS
        uint32_t ctr = 0;
        for(uint32_t j = 0; j < N-K; j++) {
            while (is_pivot_column[ctr]) {
                ctr += 1;
            }
            /// copy column
            for (uint32_t k = 0; k < K; k++) {
                A_i.values[k][j] = G0.values[k][ctr];
            }
            ctr += 1;
        }

        POSITION_T piv_idx = 0;
        for(uint32_t col_idx = 0; col_idx < N; col_idx++) {
            POSITION_T row_idx = 0;
            for(uint32_t t = 0; t < N; t++) {
               if (mu_tilde.permutation[t] == col_idx) {
                   row_idx = t;
                   break;
               }
            }
            if(is_pivot_column[col_idx] == 1) {
               pi_tilde[i].permutation[piv_idx] = row_idx;
               piv_idx++;
            }
        }

        blind(&A_i, &cf_shake_state);
        const int t = CF(&A_i);

        if (t == 0) {
            *(linearized_rounds_seeds + i*SEED_LENGTH_BYTES) += 1;
            i -= 1;
        } else {
            // NOTE: as we increase the size of the `normalized_IS_t`
            // we need to hash the values row by row.
#if defined(USE_AVX2) || defined(USE_NEON)
            for (uint32_t sl = 0; sl < K; sl++) {
                LESS_SHA3_INC_ABSORB(&state, A_i.values[sl], K);
            }
#else
            LESS_SHA3_INC_ABSORB(&state, (uint8_t *)&A_i, sizeof(normalized_IS_t));
#endif
        }
    }

    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *)m, mlen);
    LESS_SHA3_INC_ABSORB(&state, sig->salt, HASH_DIGEST_LENGTH);

    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(sig->digest, &state);
    // (b_0, ..., b_{t-1})
    uint8_t fixed_weight_string[T] = {0};
    SampleChallenge(fixed_weight_string, sig->digest);

    uint8_t indices_to_publish[T];
    for (uint32_t i = 0; i < T; i++) {
        indices_to_publish[i] = !!(fixed_weight_string[i]);
    }

    int emitted_monoms = 0;
    memset(&sig->seed_storage, 0, SEED_TREE_MAX_PUBLISHED_BYTES);

    const uint32_t num_seeds_published =
            GGMPath(seed_tree,
                           indices_to_publish,
                           (unsigned char *) &sig->seed_storage);

    monomial_action_IS_t mono_action;
    for (uint32_t i = 0; i < T; i++) {
        monomial_t Q_to_multiply;
        if (fixed_weight_string[i] != 0) {
            const int sk_monom_seed_to_expand_idx = fixed_weight_string[i];

            monomial_sample_prikey(&Q_to_multiply,
                                   private_monomial_seeds[sk_monom_seed_to_expand_idx - 1]);
            monomial_compose_action(&mono_action, &Q_to_multiply, &pi_tilde[i]);

            CosetRep(sig->cf_monom_actions[emitted_monoms], &mono_action);
            emitted_monoms++;
        }
    }
    return num_seeds_published;
} /* end LESS_sign */

/// NOTE: non-constant time
/// \param PK[in]: public key
/// \param m[in]: message for which a signature was computed
/// \param mlen[in]: length of the message in bytes
/// \param sig[in]: signature
/// \return 0: on failure
///         1: on success
int LESS_verify(const pubkey_t *const PK,
                const char *const m,
                const uint64_t mlen,
                const sign_t *const sig) {
    uint8_t fixed_weight_string[T] = {0};
    uint8_t is_pivot_column[N_pad];
    uint8_t g0_initial_pivot_flags[N];
    uint8_t gi_initial_pivot_flags[N];
#ifdef LESS_REUSE_PIVOTS_VY
    uint8_t g0_permuted_pivot_flags[N];
#endif
    SampleChallenge(fixed_weight_string, sig->digest);

    uint8_t published_seed_indexes[T];
    for (uint32_t i = 0; i < T; i++) {
        published_seed_indexes[i] = !!(fixed_weight_string[i]);
    }

    unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] = {0};
    uint32_t rebuilding_seeds_went_fine;
    rebuilding_seeds_went_fine = RebuildGGM(seed_tree,
                                                          published_seed_indexes,
                                                          (unsigned char *) &sig->seed_storage,
                                                          sig->salt);
    if (!rebuilding_seeds_went_fine) {
        return 0;
    }

    unsigned char linearized_rounds_seeds[T*SEED_LENGTH_BYTES] = {0};
    seed_leaves(linearized_rounds_seeds,seed_tree);

    int employed_monoms = 0;

    rref_generator_mat_t G0_rref;
    generator_sample(&G0_rref, PK->G_0_seed);

    generator_mat_t G0 = {0}, G0_full = {0};
    generator_mat_t G_prime = {0};
    monomial_t mu_tilde;
    normalized_IS_t Ai = {0};
    LESS_SHA3_INC_CTX state;
    LESS_SHA3_INC_INIT(&state);

    generator_get_pivot_flags(&G0_rref, g0_initial_pivot_flags);
    generator_rref_expand(&G0_full, &G0_rref);

    for (uint32_t i = 0; i < T; i++) {
        memset(is_pivot_column, 0, N_pad);
        if (fixed_weight_string[i] == 0) {
            monomial_sample_salt(&mu_tilde,
                                 linearized_rounds_seeds + i * SEED_LENGTH_BYTES,
                                 sig->salt,
                                 i);

            generator_monomial_mul(&G_prime, &G0_full, &mu_tilde);
#if defined(LESS_REUSE_PIVOTS_VY)
            uint8_t permuted_pivot_flags[N_pad] = {0};

            for (uint32_t t = 0; t < N; t++) {
                permuted_pivot_flags[mu_tilde.permutation[t]] = g0_initial_pivot_flags[t];
            }
            if (generator_RREF_pivot_reuse(&G_prime,is_pivot_column, permuted_pivot_flags, VERIFY_PIVOT_REUSE_LIMIT) == 0) {
                return 0;
            }
#else
            if (generator_RREF(&G_prime, is_pivot_column) == 0) {
                return 0;
            }
#endif
        } else {
            expand_to_rref(&G0, PK->SF_G[fixed_weight_string[i] - 1], gi_initial_pivot_flags);
            if (!CheckCanonicalAction(sig->cf_monom_actions[employed_monoms])) {
                return 0;
            }

#if defined(LESS_REUSE_PIVOTS_VY)
            apply_cf_action_to_G_with_pivots(&G_prime,
                                             &G0,
                                             sig->cf_monom_actions[employed_monoms],
                                             gi_initial_pivot_flags,
                                             g0_permuted_pivot_flags);
            const int ret = generator_RREF_pivot_reuse(&G_prime, is_pivot_column,
                                                       g0_permuted_pivot_flags,
                                                       VERIFY_PIVOT_REUSE_LIMIT);
#else
            apply_cf_action_to_G(&G_prime, &G0, sig->cf_monom_actions[employed_monoms]);
            const int ret = generator_RREF(&G_prime, is_pivot_column);
#endif
            if(ret == 0) {
                return 0;
            }

            employed_monoms++;
        }
        // just copy the non IS
        uint32_t ctr = 0;
        for(uint32_t j = 0; j < N-K; j++) {
            while (is_pivot_column[ctr]) {
                ctr += 1;
            }
            /// copy column
            for (uint32_t k = 0; k < K; k++) {
                Ai.values[k][j] = G_prime.values[k][ctr];
            }
            ctr += 1;
        }
        const int r = CF(&Ai);
        if (r == 0) {
            return 0;
        }
#if defined(USE_AVX2) || defined(USE_NEON)
        for (uint32_t sl = 0; sl < K; sl++) {
            LESS_SHA3_INC_ABSORB(&state, Ai.values[sl], K);
        }
#else
        LESS_SHA3_INC_ABSORB(&state, (uint8_t *)&Ai, sizeof(normalized_IS_t));
#endif
    }

    uint8_t recomputed_digest[HASH_DIGEST_LENGTH] = {0};
    LESS_SHA3_INC_ABSORB(&state, (const uint8_t *) m, mlen);
    LESS_SHA3_INC_ABSORB(&state, sig->salt, HASH_DIGEST_LENGTH);
    /* Squeeze output */
    LESS_SHA3_INC_FINALIZE(recomputed_digest, &state);

    return (verify(recomputed_digest, sig->digest,
                   HASH_DIGEST_LENGTH) == 0);
} /* end LESS_verify */

