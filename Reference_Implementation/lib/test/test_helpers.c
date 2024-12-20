#include <stdint.h>
#include <stdio.h>
#include "fq_arith.h"
#include "LESS.h"
#include "codes.h"
#include "monomial_mat.h"
#include "utils.h"
#include "rng.h"

// wrapper struct around the set S_n
typedef struct {
   /* considering the product GQ, permutation[...] stores into the cell with
    * index 0, the position of the DESTINATION of column 0 in G after the
    * computation of GQ.
    */
   POSITION_T permutation[N];
} permutation_t;

// wrapper struct around D_n
typedef struct {
   /* coefficients listed in order of appearance column-wise */
   FQ_ELEM coefficients[N];
} diagonal_t;


////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////

void permutation_swap(permutation_t *P, uint32_t i, uint32_t j);
void permutation_cswap(permutation_t *P, uint32_t i, uint32_t j, uintptr_t mask);
void permutation_mat_id(permutation_t *P);
void permutation_mat_rng(permutation_t *P);
void permutation_mat_id_v2(permutation_t *P, const uint32_t max);
void permutation_mat_rng_v2(permutation_t *P, const uint32_t max);
void permutation_pretty_print(const permutation_t *P);


////////////////////////////////////////////////////////////////////////
///                             Diagonal                             ///
////////////////////////////////////////////////////////////////////////
void diagonal_mat_zero(diagonal_t *D);
void diagonal_mat_id(diagonal_t *D);
void diagonal_mat_rnd(diagonal_t *D);
void diagonal_mat_id_v2(diagonal_t *D, uint32_t max);
void diagonal_mat_rnd_v2(diagonal_t *D, uint32_t max);
void diagonal_pretty_print(const diagonal_t *const D);

// defined in monomial.c
void permutation_apply_row(const permutation_t *P, normalized_IS_t *G);
void permutation_apply_col(normalized_IS_t *G, const permutation_t *P);

void diagonal_apply_row(diagonal_t *P, normalized_IS_t *G);
void diagonal_apply_col(normalized_IS_t *G, diagonal_t *P);


/* FY shuffle on the permutation, sampling from the provided PRNG state shake_monomial_state */
static inline
void yt_shuffle_state_v2(SHAKE_STATE_STRUCT *shake_monomial_state,
                         POSITION_T permutation[N],
                         const uint32_t max) {
    uint32_t rand_u32[N] = {0};
    POSITION_T tmp;

    csprng_randombytes((unsigned char *) &rand_u32, sizeof(uint32_t)*N, shake_monomial_state);
    for (size_t i = 0; i < max - 1; ++i) {
        rand_u32[i] = i + rand_u32[i] % (max - i);
    }

    for (size_t i = 0; i < max - 1; ++i) {
        tmp = permutation[i];
        permutation[i] = permutation[rand_u32[i]];
        permutation[rand_u32[i]] = tmp;
    }
}

/* FY shuffle on the permutation, sampling from the global TRNG state */
static inline
void yt_shuffle_v2(POSITION_T permutation[N], const uint32_t max) {
    yt_shuffle_state_v2(&platform_csprng_state, permutation, max);
}

////////////////////////////////////////////////////////////////////////
///                         Normalized                               ///
////////////////////////////////////////////////////////////////////////

/// NOTE: only for testing
void normalized_rng(normalized_IS_t *V) {
    randombytes((uint8_t *)V->values, K*(N-K));
    const uint8_t mask = 0x7F;
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < K; ++j) {
            V->values[i][j] &= mask;
        }
    }
}

void generator_to_normalized(normalized_IS_t *V,
                             const generator_mat_t *const G){
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < N - K; ++j) {
            V->values[i][j] = G->values[i][K+j];
        }
    }
}




////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////
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


////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////

///
void permutation_apply_col(normalized_IS_t *G,
                          const permutation_t *P) {
    for (uint32_t i = 0; i < (N-K); i++) {
        column_swap(G, i, P->permutation[i]);
    }
}

///
void permutation_apply_row(const permutation_t *P,
                           normalized_IS_t *G) {
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
/// @param P
/// @param i
/// @param j
/// @param mask
void permutation_cswap(permutation_t *P,
                       const uint32_t i,
                       const uint32_t j,
                       const uintptr_t mask) {
    ASSERT(i < K);
    ASSERT(i < N);
    MASKED_SWAP(P->permutation[i], P->permutation[j], mask);
}

///
/// @param P
void permutation_mat_id(permutation_t *P) {
    for (uint32_t i = 0; i < N; ++i) {
        P->permutation[i] = i;
    }
}

///
/// @param P
void permutation_mat_rng(permutation_t *P) {
    permutation_mat_id(P);
    yt_shuffle(P->permutation);
}

///
void permutation_mat_id_v2(permutation_t *P,
                           const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        P->permutation[i] = i;
    }
    for (uint32_t i = max; i < N; ++i) {
        P->permutation[i] = 0;
    }
}

///
void permutation_mat_rng_v2(permutation_t *P,
                            const uint32_t max) {
    permutation_mat_id_v2(P, max);
    yt_shuffle_v2(P->permutation, max);

    for (uint32_t i = max; i < N; ++i) {
        P->permutation[i] = 0;
    }
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
        while(D->coefficients[i] == 0) {
            D->coefficients[i] = fq_red(D->coefficients[i]+1);
        }
    }
}

///
void diagonal_mat_id_v2(diagonal_t *D,
                        const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        D->coefficients[i] = 1;
    }
    for (uint32_t i = max; i < N; ++i) {
        D->coefficients[i] = 0;
    }
}

///
void diagonal_mat_rnd_v2(diagonal_t *D,
                         const uint32_t max) {
    csprng_randombytes((unsigned char *) &D->coefficients, sizeof(FQ_ELEM)*max, &platform_csprng_state);
    for (uint32_t i = 0; i < max; ++i) {
        D->coefficients[i] = fq_red(D->coefficients[i]);
        while(D->coefficients[i] == 0) {
            D->coefficients[i] = fq_red(D->coefficients[i]+1);
        }
    }

    for (uint32_t i = max; i < N; ++i) {
        D->coefficients[i] = 0;
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
