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

#include "LESS.h"
#include "seedtree.h"
#include "rng.h"
#include "utils.h"

void LESS_keygen(prikey_t *SK,
                 pubkey_t *PK)
{
   /* generating private key from a single seed */
   randombytes(SK->compressed_sk, SEED_LENGTH_BYTES);
   /* expanding it onto private (inverse) monomial seeds */
   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,SK->compressed_sk,SEED_LENGTH_BYTES);

   /* The first private key monomial is an ID matrix, no need for random
    * generation, hence NUM_KEYPAIRS-1 */
   unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                         SEED_LENGTH_BYTES,
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
   for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      /* expand inverse monomial from seed */
      monomial_t private_Q;
      monomial_t private_Q_inv;
      monomial_mat_seed_expand(&private_Q_inv, monomial_seeds[i]);
      monomial_mat_inv(&private_Q, &private_Q_inv);

      generator_mat_t result_G;
      generator_monomial_mul(&result_G,
                             &tmp_full_G,
                             &private_Q);
      memset(is_pivot_column,0,sizeof(is_pivot_column));
      generator_RREF(&result_G,is_pivot_column);
      /* note that the result is stored at i-1 as the first
       * public key element is just a seed */
      compress_rref(PK->SF_G[i],
                             &result_G,
                             is_pivot_column);
   }
} /* end LESS_keygen */

void LESS_sign(const prikey_t *SK,
               const char *const m,
               const uint64_t mlen,
               sig_t *sig)
{
   /*         Private key expansion        */
   /* expand sequence of seeds for private inverse-monomial matrices */
   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,SK->compressed_sk,SEED_LENGTH_BYTES);

   /* The first private key monomial is an ID matrix, no need for random
    * generation, hence NUM_KEYPAIRS-1 */
   unsigned char private_monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(private_monomial_seeds[i],
                         SEED_LENGTH_BYTES,
                         &sk_shake_state);
   }

   /*         Ephemeral monomial generation        */
   unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
   randombytes(ephem_monomials_seed, SEED_LENGTH_BYTES);
   randombytes(sig->tree_salt, HASH_DIGEST_LENGTH);


   unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};
   generate_seed_tree_from_root(seed_tree, ephem_monomials_seed, sig->tree_salt);
   unsigned char *ephem_monomial_seeds = seed_tree +
                                         SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

   /*         Public G_0 expansion                  */
   rref_generator_mat_t G0_rref;
   generator_SF_seed_expand(&G0_rref, SK->G_0_seed);
   generator_mat_t full_G0;
   generator_rref_expand(&full_G0,&G0_rref);

   monomial_t Q_tilde;
   normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
#ifdef COMPRESS_CMT_COLUMNS
   uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
#endif
   monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));

   for(int i = 0; i < T; i++) {
      monomial_mat_seed_expand(&Q_tilde,
                               ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
      prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);
#ifdef COMPRESS_CMT_COLUMNS
      compress_columns(&V_array_compressed[i*RREF_IS_COLUMNS_PACKEDBYTES], &V_array[i]);
#endif
   }

#ifdef COMPRESS_CMT_COLUMNS
   compute_digest_array(sig->digest, m, mlen, V_array_compressed, sig->tree_salt);
#else
   compute_digest(sig->digest, m, mlen, V_array, sig->tree_salt);
#endif

   uint8_t fixed_weight_string[T] = {0};
   expand_digest_to_fixed_weight(fixed_weight_string,sig->digest);

   uint8_t indices_to_publish[T];
   for (int i=0; i<T ; i++) {
      indices_to_publish[i] = !!(fixed_weight_string[i]);
   }

   int emitted_monoms=0;
   memset(&sig->seed_storage,0,SEED_TREE_MAX_PUBLISHED_BYTES);
   publish_seeds(seed_tree,
                 indices_to_publish,
                 (unsigned char *) &sig->seed_storage);


   monomial_action_IS_t mono_action;
   for(int i = 0; i < T; i++) {
      monomial_t Q_to_multiply;
      if(fixed_weight_string[i] != 0) {
         int sk_monom_seed_to_expand_idx = fixed_weight_string[i];
         monomial_mat_seed_expand(&Q_to_multiply,
                                  private_monomial_seeds[sk_monom_seed_to_expand_idx-1]);
         monomial_compose_action(&mono_action,
                                 &Q_to_multiply,
                                 &Q_bar_actions[i]);

         compress_monom_action(sig->monom_actions[emitted_monoms], &mono_action);

         emitted_monoms++;
      }
   }
   free(V_array);
   free(Q_bar_actions);
#ifdef COMPRESS_CMT_COLUMNS
   free(V_array_compressed);
#endif
} /* end LESS_sign */

int timing_safe_memcmp(const unsigned char *a,
                       const unsigned char *b,
                       unsigned int bytelen)
{
   int are_different = 0;
   for (int i =0; i< bytelen; i++) {
      are_different |= (a[i] != b[i]);
   }
   return are_different;
}

int LESS_verify(const pubkey_t *const PK,
                const char *const m,
                const uint64_t mlen,
                const sig_t *const sig)
{

   uint8_t fixed_weight_string[T] = {0};
   expand_digest_to_fixed_weight(fixed_weight_string,sig->digest);
   monomial_action_IS_t mono_action;


   uint8_t published_seed_indexes[T];
   for (int i=0; i<T ; i++) {
      published_seed_indexes[i] = !!(fixed_weight_string[i]);
   }

   unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};
   regenerate_leaves(seed_tree,published_seed_indexes,
                     (unsigned char *)&sig->seed_storage, sig->tree_salt);

   unsigned char *ephem_monomial_seeds = seed_tree +
                                         SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

   int employed_monoms = 0, employed_seeds = 0;

   rref_generator_mat_t G0_rref;
   generator_SF_seed_expand(&G0_rref, PK->G_0_seed);

   generator_mat_t tmp_full_G;
   normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
#ifdef COMPRESS_CMT_COLUMNS
   uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
#endif

   monomial_action_IS_t Q_to_discard;
   for(int i = 0; i < T; i++) {
      if (fixed_weight_string[i] == 0) {
         generator_rref_expand(&tmp_full_G,&G0_rref);
         monomial_t Q_to_multiply = {0};
         monomial_mat_seed_expand(&Q_to_multiply,
                                  ephem_monomial_seeds+i*SEED_LENGTH_BYTES);

         prepare_digest_input(&V_array[i],
                              &Q_to_discard,
                              &tmp_full_G,
                              &Q_to_multiply);

         employed_seeds++;
      } else {
         generator_mat_t G_hat = {0};
         expand_to_rref(&tmp_full_G, PK->SF_G[fixed_weight_string[i]-1]);
         expand_to_monom_action(&mono_action, sig->monom_actions[employed_monoms]);
         
         apply_action_to_G(&G_hat,
                           &tmp_full_G,
                           &mono_action);

         uint8_t is_pivot_column[N] = {0};
         generator_RREF(&G_hat,is_pivot_column);

         /* normalize non-pivot columns into the appropriate V_array slot */
         POSITION_T placed_non_pivots = 0;

         for(int col_to_norm = 0; col_to_norm < N; col_to_norm ++) {
            if (is_pivot_column[col_to_norm] == 0) {
               lex_minimize(&V_array[i],
                            placed_non_pivots,
                            &G_hat,
                            col_to_norm);
               placed_non_pivots++;
            }
         }
         lex_sort_cols(&V_array[i]);
         employed_monoms++;
      }
#ifdef COMPRESS_CMT_COLUMNS
      compress_columns(&V_array_compressed[i*RREF_IS_COLUMNS_PACKEDBYTES], &V_array[i]);
#endif
   }
   uint8_t recomputed_digest[HASH_DIGEST_LENGTH] = {0};
#ifdef COMPRESS_CMT_COLUMNS
   compute_digest_array(recomputed_digest,m,mlen,V_array_compressed, sig->tree_salt);
#else
   compute_digest(recomputed_digest,m,mlen,V_array, sig->tree_salt);
#endif

   int does_hash_match = 0;
   does_hash_match = (timing_safe_memcmp(recomputed_digest,
                                         sig->digest,
                                         HASH_DIGEST_LENGTH) == 0);
   free(V_array);
#ifdef COMPRESS_CMT_COLUMNS
   free(V_array_compressed);
#endif
   return does_hash_match;
} /* end LESS_verify */
