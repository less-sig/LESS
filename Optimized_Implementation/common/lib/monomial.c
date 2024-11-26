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


/**
 * Taken from HQC.
 * @brief Constant-time comparison of two integers v1 and v2
 * FY shuffle on the permutation, sampling from the provided PRNG state shake_monomial_state */
static inline
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
static inline
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


/* samples a random perm matrix */
void monomial_mat_rnd(monomial_t *res) {
   fq_star_rnd_elements(res->coefficients, N);
   for(uint32_t i = 0; i < N; i++) {
      res->permutation[i] = i;
   }
   /* FY shuffle on the permutation */
   yt_shuffle(res->permutation);
} /* end monomial_mat_rnd */

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
   for(int i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->permutation[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
   fprintf(stderr,"coeffs = [");
   for(int i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->coefficients[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->coefficients[N-1]);
} /* end monomial_mat_pretty_print */

void monomial_mat_pretty_print_name(char *name, const monomial_t *to_print)
{
   fprintf(stderr,"%s = [",name);
   for(int i = 0; i < N-1; i++) {
      fprintf(stderr,"%03u, ",to_print->permutation[i]);
   }
   fprintf(stderr,"%03u ]\n",to_print->permutation[N-1]);
   fprintf(stderr,"coeffs = [");
   for(int i = 0; i < N-1; i++) {
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
 * monomial */
void monomial_compose_action(monomial_action_IS_t* out, 
                             const monomial_t * Q_in, 
                             const monomial_action_IS_t * in){
   /* to compose with monomial_action_IS_t, reverse the convention
    * for Q storage: store in permutation[i] the idx of the source column landing
    * as the i-th after the GQ product, and in coefficients[i] the coefficient 
    * by which the column is multiplied upon landing */
   monomial_t reverse_Q;
   for(int i=0;i<N;i++){
      reverse_Q.permutation[Q_in->permutation[i]] = i;
      reverse_Q.coefficients[Q_in->permutation[i]] = Q_in->coefficients[i];
   }
   /* compose actions out = Q_in*in */
   for(int i=0;i<K;i++){
      out->permutation[i] = reverse_Q.permutation[in->permutation[i]];
      out->coefficients[i] = fq_red(in->coefficients[i] * 
                             (FQ_DOUBLEPREC) reverse_Q.coefficients[in->permutation[i]]);
   }   
}

/* Compress MonomialAction object to byte array */
void compress_monom_action(uint8_t *compressed,
                            const monomial_action_IS_t *mono) {
    int encode_state = 0;
    int compress_idx = 0;

    // Compress Permutation
    for (int idx = 0; idx < K; idx++) {
#if defined(CATEGORY_1)
        compressed[compress_idx] = mono->permutation[idx];
        compress_idx++;
#elif defined(CATEGORY_3)
        switch(encode_state) {
        case 0:
            compressed[compress_idx] = mono->permutation[idx];
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 8;
            break;
        case 1:
            compressed[compress_idx] |= mono->permutation[idx] << 1;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 7;
            break;
        case 2:  
            compressed[compress_idx] |= mono->permutation[idx] << 2;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 6;               
            break;
        case 3:  
            compressed[compress_idx] |= mono->permutation[idx] << 3;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 5;                   
            break;
        case 4:
            compressed[compress_idx] |= mono->permutation[idx] << 4;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 4;                          
            break;
        case 5:
            compressed[compress_idx] |= mono->permutation[idx] << 5;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 3;                   
            break;
        case 6:
            compressed[compress_idx] |= mono->permutation[idx] << 6;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 2;                      
            break;
        case 7:
            compressed[compress_idx] |= mono->permutation[idx] << 7;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 1; 
            compress_idx++;                  
            break;
        }

        if (encode_state != 7) {
            encode_state++;
        } else {
            encode_state = 0;
        }
#else
        switch(encode_state) {
        case 0:
            compressed[compress_idx] = mono->permutation[idx];
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 8;
            break;
        case 1:
            compressed[compress_idx] |= mono->permutation[idx] << 2;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 6;
            break;
        case 2:  
            compressed[compress_idx] |= mono->permutation[idx] << 4;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 4;               
            break;
        case 3:  
            compressed[compress_idx] |= mono->permutation[idx] << 6;
            compress_idx++;
            compressed[compress_idx] = mono->permutation[idx] >> 2;  
            compress_idx++;                 
            break;
        }

        if (encode_state != 3) {
            encode_state++;
        } else {
            encode_state = 0;
        }
#endif
    }

#if defined(CATEGORY_1) || defined(CATEGORY_3)
    encode_state = 0;
#else
    encode_state = 4;
#endif

    // Compress Coefficients
    for (int idx = 0; idx < K; idx++) {

        switch(encode_state) {
        case 0:
            compressed[compress_idx] = mono->coefficients[idx];
            break;
        case 1:
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 7);
            compress_idx++;
            compressed[compress_idx] = (mono->coefficients[idx] >> 1);
            break;
        case 2:   
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 6);
            compress_idx++;
            compressed[compress_idx] = (mono->coefficients[idx] >> 2);                 
            break;
        case 3:    
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 5);
            compress_idx++;
            compressed[compress_idx] = (mono->coefficients[idx] >> 3);                
            break;
        case 4:             
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 4);
            compress_idx++;
            compressed[compress_idx] = (mono->coefficients[idx] >> 4);       
            break;
        case 5:    
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 3);
            compress_idx++;
            compressed[compress_idx] = (mono->coefficients[idx] >> 5);                
            break;
        case 6:       
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 2);
            compress_idx++; 
            compressed[compress_idx] = (mono->coefficients[idx] >> 6);                   
            break;
        case 7:       
            compressed[compress_idx] = compressed[compress_idx] | (mono->coefficients[idx] << 1);
            compress_idx++;                 
            break;
        }

        if (encode_state != 7) {
            encode_state++;
        } else {
            encode_state = 0;
        }
    }
}



/* Decompress byte array to MonomialAction object */
void expand_to_monom_action(monomial_action_IS_t *mono,
                            const uint8_t *compressed) {
    uint32_t decode_state = 0;
    uint32_t compress_idx = 0;

    // Decompress Permutation
    for (uint32_t idx = 0; idx < K; idx++) {
#if defined(CATEGORY_1)
        mono->permutation[idx] = compressed[compress_idx];
        compress_idx++;
#elif defined(CATEGORY_3)
        switch(decode_state) {
        case 0:
            mono->permutation[idx] = (compressed[compress_idx] | (compressed[compress_idx+1] << 8)) & MASK_N;
            compress_idx++;
            break;
        case 1:
            mono->permutation[idx] = ((compressed[compress_idx] >> 1) | (compressed[compress_idx+1] << 7)) & MASK_N;
            compress_idx++;
            break;
        case 2:
            mono->permutation[idx] = ((compressed[compress_idx] >> 2) | (compressed[compress_idx+1] << 6)) & MASK_N;
            compress_idx++;                 
            break;
        case 3:
            mono->permutation[idx] = ((compressed[compress_idx] >> 3) | (compressed[compress_idx+1] << 5)) & MASK_N;
            compress_idx++;                 
            break;
        case 4:
            mono->permutation[idx] = ((compressed[compress_idx] >> 4) | (compressed[compress_idx+1] << 4)) & MASK_N;
            compress_idx++;                  
            break;
        case 5:
            mono->permutation[idx] = ((compressed[compress_idx] >> 5) | (compressed[compress_idx+1] << 3)) & MASK_N;
            compress_idx++;               
            break;
        case 6:
            mono->permutation[idx] = ((compressed[compress_idx] >> 6) | (compressed[compress_idx+1] << 2)) & MASK_N;
            compress_idx++;                       
            break;
        case 7:
            mono->permutation[idx] = ((compressed[compress_idx] >> 7) | (compressed[compress_idx+1] << 1)) & MASK_N;
            compress_idx += 2;                   
            break;
        }
        
        if (decode_state != 7) {
            decode_state++;
        } else {
            decode_state = 0;
        }
#else
        switch(decode_state) {
        case 0:
            mono->permutation[idx] = (compressed[compress_idx] | (compressed[compress_idx+1u] << 8)) & MASK_N;
            compress_idx++;
            break;
        case 1:
            mono->permutation[idx] = ((compressed[compress_idx] >> 2u) | (compressed[compress_idx+1u] << 6u)) & MASK_N;
            compress_idx++;
            break;
        case 2:
            mono->permutation[idx] = ((compressed[compress_idx] >> 4u) | (compressed[compress_idx+1u] << 4u)) & MASK_N;
            compress_idx++;                 
            break;
        case 3:
            mono->permutation[idx] = ((compressed[compress_idx] >> 6u) | (compressed[compress_idx+1u] << 2u)) & MASK_N;
            compress_idx += 2;                 
            break;
        }
        
        if (decode_state != 3) {
            decode_state++;
        } else {
            decode_state = 0;
        }
#endif
    } 

#if defined(CATEGORY_1) || defined(CATEGORY_3)
    decode_state = 0;
#else
    decode_state = 4;
#endif

    // Decompress Coefficients
    for (uint32_t idx = 0; idx < K; idx++) {
        switch(decode_state) {
        case 0:
            mono->coefficients[idx] = compressed[compress_idx] & MASK_Q;
            break;
        case 1:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 7) | (compressed[compress_idx+1] << 1)) & MASK_Q;
            compress_idx++;
            break;
        case 2:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 6) | (compressed[compress_idx+1] << 2)) & MASK_Q;
            compress_idx++;                 
            break;
        case 3:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 5) | (compressed[compress_idx+1] << 3)) & MASK_Q;
            compress_idx++;                 
            break;
        case 4:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 4) | (compressed[compress_idx+1] << 4)) & MASK_Q;
            compress_idx++;                  
            break;
        case 5:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 3) | (compressed[compress_idx+1] << 5)) & MASK_Q;
            compress_idx++;               
            break;
        case 6:
            mono->coefficients[idx] = ((compressed[compress_idx] >> 2) | (compressed[compress_idx+1] << 6)) & MASK_Q;
            compress_idx++;                        
            break;
        case 7:
            mono->coefficients[idx] = (compressed[compress_idx] >> 1) & MASK_Q;
            compress_idx++;                   
            break;
        }
        
        if (decode_state != 7) {
            decode_state++;
        } else {
            decode_state = 0;
        }
    }

}

int is_monom_action_valid(const monomial_action_IS_t * const mono){
    for(int i = 0; i < K; i++){
        if ((mono->coefficients[i] <= 0) || (mono->coefficients[i] >= Q)){
            return 0;
        }
    }
    uint8_t flags[N] = {0};
    for(int i = 0; i < K; i++){
        if ((mono->permutation[i] < 0) || (mono->permutation[i] >= N)){
            return 0;
        }
        if (flags[mono->permutation[i]] == 1 ){
            return 0;
        }        
        flags[mono->permutation[i]] = 1;
    }
    return 1;
}


////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////


void permutation_apply_col(normalized_IS_t *G, permutation_t *P) {
    for (uint32_t i = 0; i < (N-K); i++) {
        column_swap(G, i, P->permutation[i]);
    }
}

///
void permutation_apply_row(permutation_t *P, normalized_IS_t *G) {
    for (uint32_t i = 0; i < K; i++) {
        row_swap(G, i, P->permutation[i]);
    }
}

///
void permutation_swap(permutation_t *P,
                      const uint32_t i,
                      const uint32_t j) {
    ASSERT(i < K);
    ASSERT(i < N);
    POSITION_T tmp = P->permutation[i];
    P->permutation[i] = P->permutation[j];
    P->permutation[j] = tmp;
}

///
void permutation_cswap(permutation_t *P,
                       const uint32_t i,
                       const uint32_t j,
                       const uintptr_t mask) {
    ASSERT(i < K);
    ASSERT(i < N);
    MASKED_SWAP(P->permutation[i], P->permutation[j], mask);
}

///
void permutation_mat_id(permutation_t *P) {
    for (uint32_t i = 0; i < N; ++i) {
        P->permutation[i] = i;
    }
}

///
void permutation_mat_rng(permutation_t *P) {
    permutation_mat_id(P);
    yt_shuffle(P->permutation);
}

///
void permutation_mat_id_v2(permutation_t *P, const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        P->permutation[i] = i;
    }
}

/// TODO remove
void permutation_mat_rng_v2(permutation_t *P, const uint32_t max) {
    permutation_mat_id_v2(P, max);
    // yt_shuffle_v2(P->permutation, max);
}

///
void permutation_pretty_print(const permutation_t *const P) {
    fprintf(stderr,"perm = [");
    for(uint32_t i = 0; i < N-1; i++) {
        fprintf(stderr,"%03u, ", P->permutation[i]);
    }

    fprintf(stderr,"%03u ]\n", P->permutation[N-1]);
}

////////////////////////////////////////////////////////////////////////
///                             Diagonal                             ///
////////////////////////////////////////////////////////////////////////

void diagonal_apply_col(normalized_IS_t *G,
                        diagonal_t *P) {
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < (N-K); j++) {
            G->values[i][j] = fq_mul(G->values[i][j], P->coefficients[j]);
        }
    }
}

///
void diagonal_apply_row(diagonal_t *P,
                        normalized_IS_t *G) {
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < (N-K); j++) {
            G->values[i][j] = fq_mul(G->values[i][j], P->coefficients[i]);
        }
    }
}

///
void diagonal_mat_zero(diagonal_t *D) {
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = 0;
    }
}

///
void diagonal_mat_id(diagonal_t *D) {
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = 1;
    }
}

///
void diagonal_mat_rnd(diagonal_t *D) {
    csprng_randombytes((unsigned char *) &D->coefficients, sizeof(FQ_ELEM)*N, &platform_csprng_state);
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = fq_red(D->coefficients[i]);
    }
}

///
void diagonal_mat_id_v2(diagonal_t *D,
                        const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        D->coefficients[i] = 1;
    }
}

///
void diagonal_mat_rnd_v2(diagonal_t *D,
                         const uint32_t max) {
    csprng_randombytes((unsigned char *) &D->coefficients, sizeof(FQ_ELEM)*max, &platform_csprng_state);
    for (uint32_t i = 0; i < max; ++i) {
        D->coefficients[i] = fq_red(D->coefficients[i]);
    }
}

///
void diagonal_pretty_print(const diagonal_t *const P) {
    fprintf(stderr,"diag = [");
    for(uint32_t i = 0; i < N-1; i++) {
        fprintf(stderr,"%03u, ", P->coefficients[i]);
    }

    fprintf(stderr,"%03u ]\n", P->coefficients[N-1]);
}
