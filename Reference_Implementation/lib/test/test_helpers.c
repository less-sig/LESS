#include <stdint.h>
#include <stdio.h>
#include <string.h>


#include "sort.h"
#include "transpose.h"
#include "fq_arith.h"
#include "codes.h"
#include "monomial_mat.h"
#include "utils.h"
#include "rng.h"
#include "test_helpers.h"

void column_swap(normalized_IS_t *V,
                 const POSITION_T col1,
                 const POSITION_T col2){
   for(uint32_t i = 0; i<K;i++ ){
      POSITION_T tmp;
      tmp = V->values[i][col2];
      V->values[i][col2] = V->values[i][col1];
      V->values[i][col1] = tmp;
   }
}

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
///                          Generator                               ///
////////////////////////////////////////////////////////////////////////

/* samples a random generator matrix */
void generator_rnd(generator_mat_t *res) {
   for(uint32_t i = 0; i < K; i++) {
      rand_range_q_elements(res->values[i], N);
   }
} /* end generator_rnd */

// generate a random matrix with full rank, where the first k columns are systemized
void generator_sf(generator_mat_t *res) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < K; ++j) {
            res->values[i][j] = i == j;
        }

        rand_range_q_elements(res->values[i] + K, N-K);
    }
}

/// \param G[in]:
void generator_pretty_print(const generator_mat_t *const G) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < N-1; ++j) {
            printf("%3d,", G->values[i][j]);
        }
        printf("%3d\n", G->values[i][N-1]);
    }

    printf("\n");
}

/* pretty_print for full generator matrices */
void generator_pretty_print_name(char *name,
                                 const generator_mat_t *const G) {
    fprintf(stderr,"%s = M([",name);
    for(uint32_t i = 0; i < K-1 ; i++ ) {
        fprintf(stderr,"[");
        for(uint32_t j = 0; j < N-1; j++) {
            fprintf(stderr,"%u, ",G->values[i][j]);
        }
        fprintf(stderr,"%u ],\n",G->values[i][N-1]);
    }
    fprintf(stderr,"[");
    for(uint32_t j = 0; j < N-1; j++) {
        fprintf(stderr,"%u, ",G->values[K-1][j]);
    }
    fprintf(stderr,"%u ] ])\n",G->values[K-1][N-1]);
} /* end generator_pretty_print_name */

/* pretty_print for generator matrices in row-reduced echelon form*/
void generator_rref_pretty_print_name(char *name,
                                      const rref_generator_mat_t *const G)
{
   fprintf(stderr,"%s =\n[",name);
   for(uint32_t i = 0; i < K-1 ; i++ ) {
      fprintf(stderr,"[");
      for(uint32_t j = 0; j < (N-K)-1; j++) {
         fprintf(stderr,"%u, ",G->values[i][j]);
      }
      fprintf(stderr,"%u ],\n",G->values[i][(N-K)-1]);
   }
   fprintf(stderr,"[");
   for(uint32_t j = 0; j < (N-K)-1; j++) {
      fprintf(stderr,"%u, ",G->values[K-1][j]);
   }
   fprintf(stderr,"%u ] ]\n",G->values[K-1][(N-K)-1]);
   fprintf(stderr,"column_pos = \n [ ");
   for(uint32_t x=0; x < K ; x++) {
      fprintf(stderr," %d ",G->column_pos[x]);
   }
   fprintf(stderr,"]\n");

} /* end generator_rref_pretty_print_name */

////////////////////////////////////////////////////////////////////////
///                         Normalized                               ///
////////////////////////////////////////////////////////////////////////

/// NOTE: only for testing
void normalized_rng(normalized_IS_t *V) {
    randombytes((uint8_t *)V->values, K*(N-K));
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < K; ++j) {
            V->values[i][j] %= Q;
        }
    }
}

/// generates a K \times (N-K) identity matrix.
/// \param V[in/out]
void normalized_ind(normalized_IS_t *V) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < N-K; ++j) {
            V->values[i][j] = i == j;
        }
    }
}

/// generates a
/// \param V [in/out]
void normalized_sf(normalized_IS_t *V) {
    normalized_ind(V);

    unsigned char x;
    for (uint32_t b = 0; b < 32; b++) {
        for (uint32_t i = 0; i < K; ++i) {
            for (uint32_t j = 0; j < K; ++j) {
                if (j == i) { continue; }

                randombytes(&x, 1);
                if (x & 1) {
                    for (uint32_t k = 0; k < N - K; ++k) {
                        if ((b&1) == 0) V->values[j][k] = fq_add(V->values[j][k], V->values[i][k]);
                        else V->values[K-1-j][k] = fq_add(V->values[K-1-j][k], V->values[K-1-i][k]);
                    }
                }
            }
        }
    }
}

/// \param V
/// \param G
void generator_to_normalized(normalized_IS_t *V,
                             const generator_mat_t *const G){
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < N - K; ++j) {
            V->values[i][j] = G->values[i][K+j];
        }
    }
}

/// \param G
void normalized_pretty_print(const normalized_IS_t *const G) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < (N-K-1); ++j) {
            printf("%3d,", G->values[i][j]);
        }
        printf("%3d\n", G->values[i][N-K-1]);
    }

    printf("\n");
}

/// \param values
void normalized_pretty_print_v(const FQ_ELEM values[K][N-K]) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < (N-K-1); ++j) {
            printf("%3d,", values[i][j]);
        }
        printf("%3d\n", values[i][N-K-1]);
    }

    printf("\n");
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

////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////

/// \param G[in/out];
/// \param P[in]:
void permutation_apply_col(normalized_IS_t *G,
                          const permutation_t *P) {
    for (uint32_t i = 0; i < (N-K); i++) {
        column_swap(G, i, P->permutation[i]);
    }
}

/// \param P[in]:
/// \param G[in/out];
void permutation_apply_row(const permutation_t *P,
                           normalized_IS_t *G) {
    for (uint32_t i = 0; i < K; i++) {
        normalized_row_swap(G, i, P->permutation[i]);
    }
}

/// \param P
/// \param i
/// \param j
void permutation_swap(permutation_t *P,
                      const uint32_t i,
                      const uint32_t j) {
    POSITION_T tmp = P->permutation[i];
    P->permutation[i] = P->permutation[j];
    P->permutation[j] = tmp;
}

/// \param P
void permutation_mat_id(permutation_t *P) {
    for (uint32_t i = 0; i < N; ++i) {
        P->permutation[i] = i;
    }
}

/// \param P
void permutation_mat_rng(permutation_t *P) {
    permutation_mat_id(P);
    yt_shuffle(P->permutation);
}

/// \param P
/// \param max
void permutation_mat_id_v2(permutation_t *P,
                           const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        P->permutation[i] = i;
    }
    for (uint32_t i = max; i < N; ++i) {
        P->permutation[i] = 0;
    }
}

/// \param P
/// \param max
void permutation_mat_rng_v2(permutation_t *P,
                            const uint32_t max) {
    permutation_mat_id_v2(P, max);
    yt_shuffle_v2(P->permutation, max);

    for (uint32_t i = max; i < N; ++i) {
        P->permutation[i] = 0;
    }
}

/// \param P
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

/// \param G[in/out]:
/// \param P[in]:
void diagonal_apply_col(normalized_IS_t *G,
                        const diagonal_t *P) {
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < (N-K); j++) {
            G->values[i][j] = fq_mul(G->values[i][j], P->coefficients[j]);
        }
    }
}

/// \param P[in]:
/// \param G[in/out]:
void diagonal_apply_row(const diagonal_t *P,
                        normalized_IS_t *G) {
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < (N-K); j++) {
            G->values[i][j] = fq_mul(G->values[i][j], P->coefficients[i]);
        }
    }
}

/// \param D[in]:
void diagonal_mat_zero(diagonal_t *D) {
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = 0;
    }
}

/// \param D[in]:
void diagonal_mat_id(diagonal_t *D) {
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = 1;
    }
}

/// \param D[in]:
void diagonal_mat_rnd(diagonal_t *D) {
    csprng_randombytes((unsigned char *) &D->coefficients, sizeof(FQ_ELEM)*N, &platform_csprng_state);
    for (uint32_t i = 0; i < N; ++i) {
        D->coefficients[i] = fq_red(D->coefficients[i]);
        while(D->coefficients[i] == 0) {
            D->coefficients[i] = fq_red(D->coefficients[i]+1);
        }
    }
}

/// \param D[in]:
/// \param max[in]:
void diagonal_mat_id_v2(diagonal_t *D,
                        const uint32_t max) {
    for (uint32_t i = 0; i < max; ++i) {
        D->coefficients[i] = 1;
    }
    for (uint32_t i = max; i < N; ++i) {
        D->coefficients[i] = 0;
    }
}

/// \param D[in]:
/// \param max[in]:
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

/// \param D[in]:
void diagonal_pretty_print(const diagonal_t *const D) {
    fprintf(stderr,"diag = [");
    for(uint32_t i = 0; i < N-1; i++) {
        fprintf(stderr,"%03u, ", D->coefficients[i]);
    }

    fprintf(stderr,"%03u ]\n", D->coefficients[N-1]);
}

////////////////////////////////////////////////////////////////////////
///                          Sorting                                 ///
////////////////////////////////////////////////////////////////////////

/// taken from djbsort
/// \param a[in/out] first input
/// \param b[in/out] second input
/// \returns nothing but: a = min(a, b),
///						  b = max(a, b)
#define int8_MINMAX(a,b)                \
do {                                    \
    int8_t ab = (int8_t)b ^ (int8_t)a;  \
    int8_t c  = (int8_t)b - (int8_t)a;  \
    c ^= ab & (c ^ (int8_t)b);          \
    c >>= 7;                            \
    c &= ab;                            \
    a ^= (FQ_ELEM)c;                    \
    b ^= (FQ_ELEM)c;                    \
} while(0)

/// NOTE: taken from djbsort
/// \param x[in/out] input array
/// \param n[in] length
void bitonic_sort_i8(FQ_ELEM *x,
                     const long long n) {
    long long p, q, r, i;
    if (n < 2) {
        return;
    }
    const long long top = 1ul << (32 - __builtin_clz(K/2));

    for (p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                int8_MINMAX(x[i], x[i + p]);
            }
        }

        i = 0;
        for (q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    FQ_ELEM a = x[i + p];
                    for (r = q; r > p; r >>= 1) {
                        int8_MINMAX(a, x[i + r]);
                    }

                    x[i + p] = a;
                }
            }
        }
    }
}

/// NOTE: helper function for type3 canonical form.
/// NOTE: specially made for the bitonic sort, which
///		operates on pointers.
/// \input rows[in/out]: K x (N-K) matrix
/// \input row1[in]: first row to compare
/// \input row2[in]: secnd row to compare
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
int compare_rows_bitonic_sort(FQ_ELEM **rows,
							  const uint32_t row1,
							  const uint32_t row2) {
    uint32_t i = 0;
    while((i < (N-K-1)) && (rows[row1][i] == rows[row2][i])) {
        i += 1;
    }
    return (((int)(rows[row1][i])) - ((int)(rows[row2][i])));
}

/// \param ptr[in/out]:
/// \param P[in/out]: a permutation to keep track of the sorting
/// \param start[in]: inclusive
/// \param end[in]: inclusive
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_quick_sort_recursive_internal(FQ_ELEM* ptr[K],
                                     uint32_t P[K],
                                     const uint32_t start,
                                     const uint32_t end) {
    if(start < end){
        const uint32_t p = 0; //(uint32_t)row_quick_sort_internal_hoare_partition(ptr, P, start, end);
    	if (p == -1u) { return 0; }
        row_quick_sort_recursive_internal(ptr, P, start, p);
        row_quick_sort_recursive_internal(ptr, P, p + 1u, end);
    }

	return 1;
}

/// NOTE: only operates on ptrs
/// NOTE: not constant time
/// \param G[in/out]: generator matrix to sort
/// \param n[in] number of elements to sort
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort_recursive(normalized_IS_t *G,
                             const uint32_t n) {
    // first sort each row into a tmp buffer
    FQ_ELEM tmp[K][Q];
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < n; ++i) {
        sort(tmp[i], G->values[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
    }

    const int ret = row_quick_sort_recursive_internal(ptr, P, 0,  n - 1u);
    if (ret == 0) { return 0; }

    // apply the permutation
    for (uint32_t t = 0; t < n; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        normalized_row_swap(G, t, ind);
    }

    return 1;
}

/// lexicographic comparison
/// \return G1[col1] <=> G2[col2]:
///         -1: G1[col1] > G2[col2]
///          0: G1[col1] == G2[col2]
///          1: G1[col1] < G2[col2]
int lex_compare_column(const generator_mat_t *G1,
					   const generator_mat_t *G2,
                       const POSITION_T col1,
                       const POSITION_T col2) {
   uint32_t i=0;
   while((i < K) &&
         (G1->values[i][col1]-G2->values[i][col2] == 0)) {
       i++;
   }

   if (i >= K) return 0;

   if (G1->values[i][col1]-G2->values[i][col2] > 0){
      return -1;
   }

   return 1;
}

/// lexicographic comparison
/// \return G1[col1] <=> G1[col2]:
///           1: G1[col1] >  G1[col2]
///           0: G1[col1] == G1[col2]
///          -1: G1[col1] <  G1[col2]
int lex_compare_col(const normalized_IS_t *G1,
                    const POSITION_T col1,
                    const POSITION_T col2) {
   uint32_t i=0;
   while((i < (K-1)) &&
         (G1->values[i][col1]-G1->values[i][col2] == 0)) {
       i++;
   }
   return G1->values[i][col1]-G1->values[i][col2];
}

/* lexicographic comparison of a column with the pivot
 * returns 1 if the pivot is greater, -1 if it is smaller,
 * 0 if it matches */
int lex_compare_with_pivot(normalized_IS_t *V,
                           const POSITION_T col_idx,
                           const FQ_ELEM pivot[K]){
   uint32_t i=0;
   while(i<K && V->values[i][col_idx]-pivot[i] == 0){
       i++;
   }
   if (i==K) return 0;
   if (V->values[i][col_idx]-pivot[i] > 0){
      return -1;
   }
   return 1;
}

/// NOTE: only used in `col_lex_quicksort`
int Hoare_partition(normalized_IS_t *V,
                    const POSITION_T col_l,
                    const POSITION_T col_h){
    FQ_ELEM pivot_col[K];
    for(uint32_t i = 0; i < K; i++){
       pivot_col[i] = V->values[i][col_l];
    }
    // NOTE double comparison, but thats no so important
    // as this function is only for testing

    POSITION_T i = col_l, j = col_h+1;
    do {
        j--;
    } while(lex_compare_with_pivot(V,j,pivot_col) == -1);
    if(i >= j){
        return j;
    }

    column_swap(V,i,j);

    while(1){
        do {
            i++;
        } while(lex_compare_with_pivot(V,i,pivot_col) == 1);
        do {
            j--;
        } while(lex_compare_with_pivot(V,j,pivot_col) == -1);
        if(i >= j){
            return j;
        }

        column_swap(V,i,j);
    }
}

/* In-place quicksort */
void col_lex_quicksort(normalized_IS_t *V,
                       int start,
                       int end){
    if(start < end){
        int p = Hoare_partition(V,start,end);
        col_lex_quicksort(V,start,p);
        col_lex_quicksort(V,p+1,end);
    }
}

/* Sorts the columns of V in lexicographic order */
void lex_sort_cols(normalized_IS_t *V){
   col_lex_quicksort(V,0,(N-K)-1);
}

////////////////////////////////////////////////////////////////////////
///                         Canonical Form                           ///
////////////////////////////////////////////////////////////////////////


/// NOTE: constant time impl
/// computes the result inplace
/// first sort the rows, then the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3_ct(normalized_IS_t *G) {
    // NOTE: not working const int ret = row_bitonic_sort(G);
    // but this is just for testing and will be removed.
    const int ret = SortRows(G, K, NULL);
    SortCols(G, K);
    // col_lex_quicksort(G, 0, N-K-1);
    return ret;
}



/// NOTE: constant time
/// NOTE: computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4_ct(normalized_IS_t *G) {
	for (uint32_t row = 0; row < K; row++) {
		if (row_all_same(G->values[row])) { continue; }

        // if we cant find a power
		FQ_ELEM s = row_acc(G->values[row]);
		FQ_ELEM sp = row_acc_inv(G->values[row]);

		if (s != 0) {
			s = fq_inv_non_ct(s);
		} else {
			s = sp;
			if (s == 0) {
				return 0;
			}
		}

		row_mul(G->values[row], s);
	}

	return compute_canonical_form_type3_ct(G);
}
/// NOTE: constant time
/// NOTE: computes the result inplace
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_ct(normalized_IS_t *G) {
    normalized_IS_t Aj, smallest;
    int touched = 0;

    // init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, sizeof(normalized_IS_t));
	memset(&Aj.values, 0, sizeof(normalized_IS_t));

    FQ_ELEM row_inv_data[N_K_pad] = {0};
    for (uint32_t row = 0; row < K; row++) {
        if (row_contains_zero(G->values[row])) { continue; }

        row_inv2(row_inv_data, G->values[row]);
        for (uint32_t row2 = 0; row2 < K; row2++) {
            row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
        }

        const int ret = compute_canonical_form_type4_ct(&Aj);
        if ((ret == 1) && (compare_matrices(&Aj, &smallest, K) < 0)) {
            touched = 1;
            normalized_copy(&smallest, &Aj);
        }
    }

    if (!touched) { return 0; }

    normalized_copy(G, &smallest);
    return 1;
}

/// NOTE: only for "fastest" CF computation
/// NOTE: non-constant time
/// implements a total order on matrices
/// we simply compare the columns lexicographically
/// \param V1[in]
/// \param V2[in]
/// \param z[in]: number of rows within both matrices
/// \return -x if V2 > V1
///			 0 if V2 == V1
///			+x
int compare_matrices_fastest(const normalized_IS_t *__restrict__ V1,
                             const normalized_IS_t *__restrict__ V2,
                             const uint32_t z) {
	for (uint32_t row = 0; row < z; row++) {
        uint32_t i=0;
        while((i < N-K) &&
             ((V1->values[row][i] - V2->values[row][i]) == 0)) {
            i++;
        }

        if (i >= K) { continue; }

        return (int)(V2->values[row][i]) - (int)(V1->values[row][i]);
	}

	// if we are here the two matrices are equal
	return 0;
}

/// NOTE: non-constant time
/// Applies the row scaling to only z rows
/// \param G[in/out] sub matrix with only z rows
/// \param z[in] number of rows in G
/// \return 0: if a row cannot be scaled any further
///         1: on success
int compute_canonical_form_type4_sub_v2(normalized_IS_t *G,
                                        const uint32_t z) {
    for (uint32_t i = 0; i < z; i++) {
		FQ_ELEM s = row_acc(G->values[i]);

		if (s != 0) {
			s = fq_inv_non_ct(s);
		} else {
			s = row_acc_inv(G->values[i]);
			if (s == 0) { return 0; }
		}

		row_mul(G->values[i], s);
    }

    return 1;
}


/****************************** Pretty Printers ******************************/
#include <stdio.h>

/* pretty-prints a seed */
void pseed(unsigned char seed[SEED_LENGTH_BYTES]){
     fprintf(stderr,"-");
   for (int i = 0 ; i < SEED_LENGTH_BYTES; i++){
     fprintf(stderr,"%02X", seed[i]);
   }
     fprintf(stderr,"- ");
}

/* pretty-prints a salt */
void psalt(unsigned char salt[SALT_LENGTH_BYTES]){
     fprintf(stderr,"-");
   for (int i = 0 ; i < SALT_LENGTH_BYTES; i++){
     fprintf(stderr,"%02X", salt[i]);
   }
     fprintf(stderr,"- ");
}

/* pretty-prints a seed tree */
void ptree(unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES]){
   const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;
   int node_idx =0;
   fprintf(stderr,"Tree dump\n");
   int ancestors = 0;
   for (int level = 0; level < LOG2(T)+1; level++){
      fprintf(stderr," * Level %d \n", level);
      for (int idx_in_level = 0; idx_in_level < npl[level]; idx_in_level++ ) {
          node_idx = ancestors + idx_in_level ;
          fprintf(stderr," [%d] ",node_idx);
          pseed(seed_tree+node_idx*SEED_LENGTH_BYTES);
      }
      ancestors += npl[level];
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");
}

/* pretty-prints a seed tree pairing each node with the corresponding flag
 * in the stencil flag tree */
void p_fs_tree(unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
               unsigned char flag_tree[NUM_NODES_SEED_TREE]){
   const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;
   int node_idx =0;
   fprintf(stderr,"Tree dump\n");
   int ancestors = 0;
   for (int level = 0; level < LOG2(T)+1; level++){
      fprintf(stderr," * Level %d \n", level);
      for (int idx_in_level = 0; idx_in_level < npl[level]; idx_in_level++ ) {
          node_idx = ancestors + idx_in_level ;
          fprintf(stderr," [%d:f=%d] ",node_idx, flag_tree[node_idx]);
          pseed(seed_tree+node_idx*SEED_LENGTH_BYTES);
      }
      ancestors += npl[level];
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");
}


/************************ End of Pretty Printers ******************************/
