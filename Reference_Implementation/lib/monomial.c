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
#include "monomial_mat.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <assert.h>
#include <string.h>

#define POS_BITS BITS_TO_REPRESENT(N-1)
#define POS_MASK (((POSITION_T) 1 << POS_BITS) - 1)

///
/// \param shake_monomial_state
/// \param permutation
/// \param n <= N
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
void yt_shuffle_state(SHAKE_STATE_STRUCT *shake_monomial_state, POSITION_T permutation[N]) {
    uint32_t rand_u32[N] = {0};
    POSITION_T tmp;

    csprng_randombytes((unsigned char *) &rand_u32, sizeof(uint32_t)*N, shake_monomial_state);
    for (size_t i = 0; i < N - 1; ++i) {
        rand_u32[i] = i + rand_u32[i] % (N - i);
    }

    for (size_t i = 0; i < N - 1; ++i) {
        tmp = permutation[i];
        permutation[i] = permutation[rand_u32[i]];
        permutation[rand_u32[i]] = tmp;
    }
}

/* FY shuffle on the permutation, sampling from the global TRNG state */
void yt_shuffle(POSITION_T permutation[N]) {
    yt_shuffle_state(&platform_csprng_state, permutation);
}

/* expands a monomial matrix, given a PRNG seed and a salt (used for ephemeral
 * monomial matrices */
void monomial_mat_seed_expand_salt_rnd(monomial_t *res,
                                       const unsigned char seed[SEED_LENGTH_BYTES],
                                       const unsigned char salt[HASH_DIGEST_LENGTH],
                                       const uint16_t round_index)
{
   SHAKE_STATE_STRUCT shake_monomial_state = {0};
   const int shake_buffer_len = SEED_LENGTH_BYTES+HASH_DIGEST_LENGTH+sizeof(uint16_t);
   uint8_t shake_input_buffer[shake_buffer_len];
   memcpy(shake_input_buffer,seed,SEED_LENGTH_BYTES);  
   memcpy(shake_input_buffer+SEED_LENGTH_BYTES,
          salt,
          HASH_DIGEST_LENGTH);
   memcpy(shake_input_buffer+SEED_LENGTH_BYTES+HASH_DIGEST_LENGTH,
          &round_index,
          sizeof(uint16_t));
   
   initialize_csprng(&shake_monomial_state,shake_input_buffer,shake_buffer_len);
   fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = i;
   }

   /* FY shuffle on the permutation */
   yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */

/* expands a monomial matrix, given a double length PRNG seed (used to prevent
 * multikey attacks) */
void monomial_mat_seed_expand_prikey(monomial_t *res,
                                     const unsigned char seed[PRIVATE_KEY_SEED_LENGTH_BYTES])
{
   SHAKE_STATE_STRUCT shake_monomial_state = {0};
   initialize_csprng(&shake_monomial_state,seed,PRIVATE_KEY_SEED_LENGTH_BYTES);
   fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = i;
   }
   /* FY shuffle on the permutation */
   yt_shuffle_state(&shake_monomial_state, res->permutation);
} /* end monomial_mat_seed_expand */


void monomial_mat_seed_expand_rnd(monomial_t *res,
                                  const unsigned char seed[SEED_LENGTH_BYTES],
                                  const uint16_t round_index) {
    SHAKE_STATE_STRUCT shake_monomial_state = {0};
    const int shake_buffer_len = SEED_LENGTH_BYTES+HASH_DIGEST_LENGTH+sizeof(uint16_t);
    uint8_t shake_input_buffer[shake_buffer_len];
    memcpy(shake_input_buffer,seed,SEED_LENGTH_BYTES);
    memcpy(shake_input_buffer+SEED_LENGTH_BYTES,
           &round_index,
           sizeof(uint16_t));

    initialize_csprng(&shake_monomial_state,shake_input_buffer,shake_buffer_len);
    fq_star_rnd_state_elements(&shake_monomial_state, res->coefficients, N);
    for(uint32_t i = 0; i < N; i++) {
        res->permutation[i] = i;
    }

    /* FY shuffle on the permutation */
    yt_shuffle_state(&shake_monomial_state, res->permutation);

}

/* samples a random perm matrix */
void monomial_mat_rnd(monomial_t *res) {
   fq_star_rnd_elements(res->coefficients, N);
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = i;
   }
   /* FY shuffle on the permutation */
   yt_shuffle(res->permutation);
} /* end monomial_mat_rnd */

// samples a random monomial matrix, in which each row has
// its unique multiset spanning. ( <=> pairwise rows do not have the same values)
void monomial_mat_rnd_unique(monomial_t *res) {
    monomial_mat_rnd(res);

    res->coefficients[0] = 1;
    for(uint32_t row = 1; row < K; row++) {
        res->coefficients[row] = row;
    }

    res->coefficients[K] = 2;
    for(uint32_t row = 1; row < K; row++) {
        res->coefficients[K + row] = row;
    }
}

/// res = A*B
/// \param res[out] output monomial matrix
/// \param A[in] input monomial matrix
/// \param B[in] input monomial matrix
void monomial_mat_mul(monomial_t *res,
                      const monomial_t *const A,
                      const monomial_t *const B) {
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = B->permutation[A->permutation[i]];
      res->coefficients[i] = fq_red(
                                (FQ_DOUBLEPREC) A->coefficients[i] *
                                (FQ_DOUBLEPREC) B->coefficients[A->permutation[i]] );
   }
} /* end monomial_mat_mul */

///
/// \param res[out] = to_invert**-1
/// \param to_invert[in]
void monomial_mat_inv(monomial_t *res,
                      const monomial_t *const to_invert) {
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[to_invert->permutation[i]] = i;
      res->coefficients[to_invert->permutation[i]] = fq_inv(
               to_invert->coefficients[i]);
   }
} /* end monomial_mat_inv */

/* yields the identity matrix */
void monomial_mat_id(monomial_t *res) {
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = i;
      res->coefficients[i] = 1;
   }
} /* end monomial_mat_id */

/* pretty_print for monomial matrices */
void monomial_mat_pretty_print(const monomial_t *const to_print) {
   fprintf(stderr,"perm = [");
   for(uint32_t i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->permutation[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
   fprintf(stderr,"coeffs = [");
   for(uint32_t i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->coefficients[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->coefficients[N-1]);
} /* end monomial_mat_pretty_print */

void monomial_mat_pretty_print_name(char *name, const monomial_t *to_print)
{
   fprintf(stderr,"%s = [",name);
   for(uint32_t i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->permutation[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
   fprintf(stderr,"coeffs = [");
   for(uint32_t i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->coefficients[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->coefficients[N-1]);
} /* end monomial_mat_pretty_print_name */

void monomial_mat_print_exp_name(char *name,const monomial_t *to_print)
{
   FQ_ELEM mu[N][N]= {{0}};

   for(uint32_t i = 0; i < N; i++) {
      mu[to_print->permutation[i]][i] = to_print->coefficients[i];
   }

   fprintf(stderr,"%s = Mon([",name);
   for(uint32_t i = 0; i < N-1 ; i++ ) {
      fprintf(stderr,"[");
      for(uint32_t j = 0; j < N-1; j++) {
         fprintf(stderr,"%u, ",mu[i][j]);
      }
      fprintf(stderr,"%u ],\n",mu[i][N-1]);
   }
   fprintf(stderr,"[");
   for(uint32_t j = 0; j < N-1; j++) {
      fprintf(stderr,"%u, ",mu[N-1][j]);
   }
   fprintf(stderr,"%u ] ])\n",mu[N-1][N-1]);
} /* end monomial_mat_print_exp_name */


/* composes a compactly stored action of a monomial on an IS with a regular
 * monomial.
 * NOTE: Only the permutation is computed, as this is the only thing we need
 * since the adaption of canonical forms.
 */
void monomial_compose_action(monomial_action_IS_t* out, 
                             const monomial_t * Q_in, 
                             const monomial_action_IS_t * in){
   /* to compose with monomial_action_IS_t, reverse the convention
    * for Q storage: store in permutation[i] the idx of the source column landing
    * as the i-th after the GQ product, and in coefficients[i] the coefficient 
    * by which the column is multiplied upon landing */
   monomial_t reverse_Q;
   for(uint32_t i = 0; i < N; i++){
      reverse_Q.permutation[Q_in->permutation[i]] = i;
   }
   /* compose actions out = Q_in*in */
   for(uint32_t i = 0; i < K; i++){
      out->permutation[i] = reverse_Q.permutation[in->permutation[i]];
   }
}

/// type5 compression
/// \param compressed
/// \param mono
void cf_compress_monomial_IS_action(uint8_t *compressed,
                                    const monomial_action_IS_t *mono) {
    memset(compressed, 0, N8);
    for (uint32_t i = 0; i < K; i++) {
        const uint32_t limb = (mono->permutation[i])/8;
        const uint32_t pos  = (mono->permutation[i])%8;
        compressed[limb] ^= 1u << pos;
    }
}

/// \param mono
/// \param compressed
void cf_expand_to_monom_action(monomial_action_IS_t *mono,
                               const uint8_t *compressed) {
    for (uint32_t i = 0; i < K; i++) {
        mono->coefficients[i] = 1;
    }
    memset(mono->permutation, 0, K*sizeof(POSITION_T));

    uint32_t ctr = 0;
    for (uint32_t i = 0; i < N8; i++) {
        uint8_t tmp = compressed[i];
        while (tmp) {
            const uint32_t pos = __builtin_ctz(tmp);
            tmp ^= 1u << pos;

            mono->permutation[ctr++] = i*8 + pos;
        }
    }

    assert(ctr == K);
}

/// checks if the given (N+7/8) bytes are a valid
/// canonical form action.
/// \param mono
/// \return true: if the weight is K
///         false: if the weight is not K
int is_cf_monom_action_valid(const uint8_t* const mono) {
    uint32_t w = 0;
    for (uint32_t i = 0; i < N8; i++) {
        w += __builtin_popcount(mono[i]);
    }

    return w == K;
}
