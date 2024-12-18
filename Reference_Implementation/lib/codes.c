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

#include <string.h>
#include <stdio.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"

/// TODO remove and replace with
//          void row_mul2
// computes G[row] = a*G[row]
void scale_row(generator_mat_t *G,
               const uint32_t row,
               const FQ_ELEM a) {
	for (uint32_t col = 0; col < N; col++) {
		G->values[row][col] = fq_mul(G->values[row][col], a);
	}
}

/* Calculate pivot flag array */
void generator_get_pivot_flags (const rref_generator_mat_t *const G,
                                uint8_t pivot_flag [N]) {
    for (uint32_t i = 0; i < N; i = i + 1) {
        pivot_flag[i] = 1;
    }

    for (uint32_t i = 0; i < K; i = i + 1) {
        pivot_flag[G->column_pos[i]] = 0;
    }
}

/* right-multiplies a generator by a monomial */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *const G,
                            const monomial_t *const monom)
{
   for(uint32_t src_col_idx = 0; src_col_idx < N; src_col_idx++) {
      for(uint32_t row_idx = 0; row_idx < K; row_idx++) {
         res->values[row_idx][monom->permutation[src_col_idx]] =
            fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] *
                    (FQ_DOUBLEPREC) monom->coefficients[src_col_idx] );
      }
   }
} /* end generator_monomial_mul */


/* right-multiplies a generator by a monomial, input generator in compact form */
void rref_generator_monomial_mul(generator_mat_t *res,
                                 const generator_mat_t *G,
                                 const monomial_t *const monom)
{
   for(uint32_t src_col_idx = 0; src_col_idx < N; src_col_idx++) {
      for(uint32_t row_idx = 0; row_idx < K; row_idx++) {
         res->values[row_idx][monom->permutation[src_col_idx]] =
            fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] *
                    (FQ_DOUBLEPREC) monom->coefficients[src_col_idx] );
      }
   }
} /* end rref_generator_monomial_mul */


static inline
void swap_rows(FQ_ELEM r[N], FQ_ELEM s[N]){
   FQ_ELEM tmp;
   for(uint32_t i=0; i<N; i++) {
      tmp = r[i];
      r[i] = s[i];
      s[i] = tmp;
   }
} /* end swap_rows */

/// 
/// @param G 
/// @param is_pivot_column 
/// @return 
int generator_RREF(generator_mat_t *G,
                   uint8_t is_pivot_column[N]) {
   for(uint32_t row_to_reduce = 0; row_to_reduce < K; row_to_reduce++) {
      uint32_t pivot_row = row_to_reduce;
      /*start by searching the pivot in the col = row*/
      uint32_t pivot_column = row_to_reduce;
      while( (pivot_column < N) &&
             (G->values[pivot_row][pivot_column] == 0) ) {

         while ( (pivot_row < K) &&
                 (G->values[pivot_row][pivot_column] == 0) ) {
            pivot_row++;
         }

         if(pivot_row >= K) { /*entire column tail swept*/
            pivot_column++; /* move to next col */
            pivot_row = row_to_reduce; /*starting from row to red */
         }
      }

      if (pivot_column >= N) {
         return 0; /* no pivot candidates left, report failure */
      }
      is_pivot_column[pivot_column] = 1; /* pivot found, mark the column*/


      /* if we found the pivot on a row which has an index > pivot_column
       * we need to swap the rows */
      if (row_to_reduce != pivot_row) {
         swap_rows(G->values[row_to_reduce],G->values[pivot_row]);
      }
      pivot_row = row_to_reduce; /* row with pivot now in place */


      /* Compute rescaling factor */
      FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);

      /* rescale pivot row to have pivot = 1. Values at the left of the pivot
       * are already set to zero by previous iterations */
      for(uint32_t i = pivot_column; i < N; i++) {
         G->values[pivot_row][i] = fq_red( (FQ_DOUBLEPREC) scaling_factor *
                                           (FQ_DOUBLEPREC) (G->values[pivot_row][i]) );
      }

      /* Subtract the now placed and reduced pivot rows, from the others,
       * after rescaling it */
      for(uint32_t row_idx = 0; row_idx < K; row_idx++) {
         if (row_idx != pivot_row) {
            FQ_DOUBLEPREC multiplier = G->values[row_idx][pivot_column];
            /* all elements before the pivot in the pivot row are null, no need to
             * subtract them from other rows. */
            for(uint32_t col_idx = 0; col_idx < N; col_idx++) {
               FQ_DOUBLEPREC tmp;
               tmp = fq_red( (FQ_DOUBLEPREC) multiplier *
                             (FQ_DOUBLEPREC) G->values[pivot_row][col_idx] );

               tmp = (FQ_DOUBLEPREC) Q + (FQ_DOUBLEPREC) G->values[row_idx][col_idx] - tmp;
               tmp = fq_red(tmp);

               G->values[row_idx][col_idx] = tmp;
            }
         }
      }
   }

   return 1;
} /* end generator_RREF */

int generator_RREF_pivot_reuse(generator_mat_t *G,
                   uint8_t is_pivot_column[N],
                   uint8_t was_pivot_column[N],
                   const int pvt_reuse_limit)
{
   int pvt_reuse_cnt;
   int row_red_pvt_skip_cnt;
   pvt_reuse_cnt = 0;

    // row swap pre-process - swap previous pivot elements to corresponding row to reduce likelihood of corruption
   int pivot_el_row;

   if (pvt_reuse_limit != 0) {
      for(int preproc_col = K-1; preproc_col >= 0; preproc_col--) {
           if (was_pivot_column[preproc_col] == 1) {
               // find pivot row
               pivot_el_row = -1;
               for (int row = 0; row < K; row = row + 1) {
                   if (G->values[row][preproc_col] != 0) {
                       pivot_el_row = row;
                   }
               }
               swap_rows(G->values[preproc_col],G->values[pivot_el_row]);
           }
      }
   }

   for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++) {
      int pivot_row = row_to_reduce;
      /*start by searching the pivot in the col = row*/
      int pivot_column = row_to_reduce;
      while( (pivot_column < N) &&
             (G->values[pivot_row][pivot_column] == 0) ) {
         while ( (pivot_row < K) &&
                 (G->values[pivot_row][pivot_column] == 0) ) {
            pivot_row++;
         }
         if(pivot_row >= K) { /*entire column tail swept*/
            pivot_column++; /* move to next col */
            pivot_row = row_to_reduce; /*starting from row to red */
         }
      }
      if ( pivot_column >=N ) {
         return 0; /* no pivot candidates left, report failure */
      }
      is_pivot_column[pivot_column] = 1; /* pivot found, mark the column*/

      /* if we found the pivot on a row which has an index > pivot_column
       * we need to swap the rows */
      if (row_to_reduce != pivot_row) {
         was_pivot_column[pivot_row] = 0; // pivot no longer reusable - will be corrupted during reduce row
         swap_rows(G->values[row_to_reduce],G->values[pivot_row]);
      }
      pivot_row = row_to_reduce; /* row with pivot now in place */

      /* Compute rescaling factor */
      FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);

      /* rescale pivot row to have pivot = 1. Values at the left of the pivot
       * are already set to zero by previous iterations */
      for(int i = pivot_column; i < N; i++) {
         G->values[pivot_row][i] = fq_red( (FQ_DOUBLEPREC) scaling_factor *
                                           (FQ_DOUBLEPREC) (G->values[pivot_row][i]) );
      }

      if (was_pivot_column[pivot_column] == 0 || (pvt_reuse_cnt >= pvt_reuse_limit) || (pivot_column >= K)) { // Skip row-reduce on previous pivots
      /* Subtract the now placed and reduced pivot rows, from the others,
       * after rescaling it */
          for(int row_idx = 0; row_idx < K; row_idx++) {
             if (row_idx != pivot_row) { 
                FQ_DOUBLEPREC multiplier = G->values[row_idx][pivot_column];
                /* all elements before the pivot in the pivot row are null, no need to
                 * subtract them from other rows. */
                row_red_pvt_skip_cnt = 0;
                for(int col_idx = 0; col_idx < N; col_idx++) {
                    if (!(col_idx < K && was_pivot_column[col_idx]) || (row_red_pvt_skip_cnt >= pvt_reuse_limit)) { // skip row reduce of pivots we will reuse
                       FQ_DOUBLEPREC tmp;
                       tmp = fq_red( (FQ_DOUBLEPREC) multiplier *
                                     (FQ_DOUBLEPREC) G->values[pivot_row][col_idx] );

                       tmp = (FQ_DOUBLEPREC) Q + (FQ_DOUBLEPREC) G->values[row_idx][col_idx] - tmp;
                       tmp = fq_red(tmp);

                       G->values[row_idx][col_idx] = tmp;
                   } else {
                     row_red_pvt_skip_cnt++;
                   }
                }
             }
          }
      } else {
         pvt_reuse_cnt++;
      }
   }


   return 1;
} /* end generator_RREF_pivot_reuse */

/* constant time lexical minimization of a column of G, writes the result
 * as the desired column of a normalized_IS_t matrix */
void lex_minimize(normalized_IS_t *V,
                  POSITION_T dst_col_idx,
                  const generator_mat_t *const G,
                  const POSITION_T col_idx)
{
   POSITION_T first_nonzero_idx = 0;
   FQ_ELEM first_nonzero_val = 0;
   int nonzero_idx_found = 0;
   for (POSITION_T i = 0; i < K; i++) {
      first_nonzero_idx = ( i *
                            ( (G->values[i][col_idx] != 0) && !(nonzero_idx_found) ) ) +
                          ( first_nonzero_idx *
                            ( (G->values[i][col_idx] == 0) || nonzero_idx_found ) );
      first_nonzero_val = ( G->values[i][col_idx] *
                            ( (G->values[i][col_idx] != 0) && !(nonzero_idx_found) ) )+
                          ( first_nonzero_val *
                            ( (G->values[i][col_idx] == 0) || nonzero_idx_found ) );
      nonzero_idx_found = nonzero_idx_found || (G->values[i][col_idx] != 0);
   }
   FQ_ELEM inv_first = fq_inv(first_nonzero_val);
   for (uint32_t i = 0; i < K; i++) {
      V->values[i][dst_col_idx] = fq_red(G->values[i][col_idx] *
                                         (FQ_DOUBLEPREC) inv_first);
   }
} /* end lex_minimize */


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

/// \param V
/// \param col1
/// \param col2
/// \param mask
void column_cswap(normalized_IS_t *V,
                  const POSITION_T col1,
                  const POSITION_T col2,
                  const uintptr_t mask){
    for(uint32_t i = 0; i<K;i++ ){
       MASKED_SWAP(V->values[i][col1], V->values[i][col2], mask);
    }
}

///
/// \param V
/// \param row1
/// \param row2
void row_swap(normalized_IS_t *V,
              const POSITION_T row1,
              const POSITION_T row2) {
    ASSERT(row1 < K);
    ASSERT(row2 < K);
    if (row1 == row2) { return; }
    for(uint32_t i = 0; i < N-K; i++){
        POSITION_T tmp;
        tmp = V->values[row1][i];
        V->values[row1][i] = V->values[row2][i];
        V->values[row2][i] = tmp;
    }
}

/// TODO avx/neon for all of this stuff
void row_cswap(normalized_IS_t *V,
              const POSITION_T row1,
              const POSITION_T row2,
              const uintptr_t mask) {
    ASSERT(row1 < K);
    ASSERT(row2 < K);

    for(uint32_t i = 0; i < N-K;i++ ){
        MASKED_SWAP(V->values[row1][i], V->values[row2][i], mask);
    }
}

/// \param V
/// \param row1
/// \param row2
void generator_row_swap(generator_mat_t *V,
                        const POSITION_T row1,
                        const POSITION_T row2) {
   for(uint32_t i = 0; i<N;i++ ){
      POSITION_T tmp;
      tmp = V->values[row1][i];
      V->values[row1][i] = V->values[row2][i];
      V->values[row2][i] = tmp;
   }
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
                           FQ_ELEM pivot[K]){
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

///
int Hoare_partition(normalized_IS_t *V, 
                    const POSITION_T col_l,
                    const POSITION_T col_h){
    FQ_ELEM pivot_col[K];
    for(uint32_t i = 0; i < K; i++){
       pivot_col[i] = V->values[i][col_l];
    }
    // TODO double comparison

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

/// TODO doc
/// \param V
/// \param Q_bar_IS
/// \param G
/// \param Q_tilde
void prepare_digest_input(normalized_IS_t *V,
                          monomial_action_IS_t *Q_bar_IS,
                          const generator_mat_t *const G,
                          const monomial_t *const Q_tilde) {
    generator_mat_t G_dagger;
    memset(&G_dagger,0,sizeof(generator_mat_t));
    generator_monomial_mul(&G_dagger, G, Q_tilde);

    uint8_t is_pivot_column[N] = {0};
    int rref_ok = generator_RREF(&G_dagger, is_pivot_column);
    /// TODO, this is kind of bad, should be removed, and proper error handling should be applied
    ASSERT(rref_ok != 0);

    // just copy the non IS
    // TODO not CT, not correct if more than 1 col is not a pivot column.
    // TODO Somehow merge with the loop just below
    uint32_t ctr = 0, offset = K;
    for(uint32_t j = 0; j < N-K; j++) {
        if (is_pivot_column[j+K]) {
            ctr += 1;
            offset = K - ctr;
        }

        for (uint32_t i = 0; i < K; i++) {
            V->values[i][j] = G_dagger.values[i][j + offset];
        }

        offset = K;
    }

    POSITION_T piv_idx = 0;
    for(uint32_t col_idx = 0; col_idx < N; col_idx++) {
        POSITION_T row_idx;
        for(uint32_t i = 0; i < N; i++) {
           if (Q_tilde->permutation[i] == col_idx) {
              row_idx = i;
           }
        }

        if(is_pivot_column[col_idx] == 1) {
           Q_bar_IS->permutation[piv_idx] = row_idx;
           piv_idx++;
        }
    }
} /* end prepare_digest_input */

/// TODO doc
/// \param V
/// \param Q_bar_IS
/// \param G
/// \param Q_tilde
/// \param initial_pivot_flags
/// \param pvt_reuse_limit
void prepare_digest_input_pivot_reuse(normalized_IS_t *V,
                                      monomial_action_IS_t *Q_bar_IS,
                                      const generator_mat_t *const G,
                                      const monomial_t *const Q_tilde,
                                      const uint8_t initial_pivot_flags [N],
                                      const int pvt_reuse_limit) {
   uint8_t g_permuated_pivot_flags[N];
   generator_mat_t G_dagger;
    // TODO unneeded
   memset(&G_dagger,0,sizeof(generator_mat_t));
   generator_monomial_mul(&G_dagger, G, Q_tilde);

   for (uint32_t i = 0; i < N; i++) {
       g_permuated_pivot_flags[Q_tilde->permutation[i]] = initial_pivot_flags[i];
   }

   uint8_t is_pivot_column[N] = {0};
   int rref_ok = generator_RREF_pivot_reuse(&G_dagger,is_pivot_column, g_permuated_pivot_flags, pvt_reuse_limit);
    /// TODO, this is kind of bad, should be removed, and proper error handling should be applied
   ASSERT(rref_ok != 0);

    // just copy the non IS
    // TODO not CT, not correct if more than 1 col is not a pivot column.
    // TODO Somehow merge with the loop just below
    uint32_t ctr = 0, offset = K;
    for(uint32_t j = 0; j < N-K; j++) {
        if (is_pivot_column[j+K]) {
            ctr += 1;
            offset = K - ctr;
        }

        for (uint32_t i = 0; i < K; i++) {
            V->values[i][j] = G_dagger.values[i][j + offset];
        }

        offset = K;
    }


    POSITION_T piv_idx = 0;
    for(uint32_t col_idx = 0; col_idx < N; col_idx++) {
        POSITION_T row_idx;
        for(uint32_t i = 0; i < N; i++) {
            if (Q_tilde->permutation[i] == col_idx) {
                row_idx = i;
            }
        }

        if(is_pivot_column[col_idx] == 1) {
            Q_bar_IS->permutation[piv_idx] = row_idx;
            piv_idx++;
        }
    }

} /* end prepare_digest_input_pivot_reuse */

/// \param res
/// \param G
/// \param Q_IS
void apply_action_to_G(generator_mat_t* res,
                       const generator_mat_t* G,
                       const monomial_action_IS_t* Q_IS,
                       uint8_t initial_G_col_pivot[N],
                       uint8_t permutated_G_col_pivot[N]) {
    /* sweep inorder Q_IS, pick cols from G, and note unpicked cols in support
     * array */
    uint8_t is_G_col_pivot[N];
    for (int i = 0; i < N; i++) {
        is_G_col_pivot[i] = 0;
    }

    for(uint32_t dst_col_idx = 0; dst_col_idx < K; dst_col_idx++){
        POSITION_T src_col_idx;
        src_col_idx = Q_IS->permutation[dst_col_idx];
        for(uint32_t i = 0; i < K; i++){
            res->values[i][dst_col_idx] = fq_red((FQ_DOUBLEPREC) G->values[i][src_col_idx] * Q_IS->coefficients[dst_col_idx]);
        }
        permutated_G_col_pivot[dst_col_idx] = initial_G_col_pivot[src_col_idx];
        is_G_col_pivot[src_col_idx] = 1;
    }

    uint32_t dst_col_idx = K;
    for(uint32_t src_col_idx = 0; src_col_idx<N; src_col_idx++){
        if (!is_G_col_pivot[src_col_idx]){
            for(uint32_t i = 0; i < K; i++){
                res->values[i][dst_col_idx] = G->values[i][src_col_idx];
            }
            permutated_G_col_pivot[dst_col_idx] = initial_G_col_pivot[src_col_idx];
            dst_col_idx++;
        }
    }    
}

/// NOTE: not constant time
/// \param res
/// /param G
/// param c
void apply_cf_action_to_G(generator_mat_t* res,
                          const generator_mat_t *G,
                          const uint8_t *const c) {
    uint32_t l = 0, r = 0;
    for (uint32_t i = 0; i < N8; i++) {
        for (uint32_t j = 0; j < 8; j++) {
            if ((i*8 + j) >= N) { goto finish; }

            const uint8_t bit = (c[i] >> j) & 1u;
            uint32_t pos;
            if (bit) {
                pos = l;
                l += 1;
            } else {
                pos = K + r;
                r += 1;
            }

            // copy the column
            for (uint32_t k = 0; k < K; k++) {
                res->values[k][pos] = G->values[k][i*8 + j];
            }
        }
    }
finish:
    assert(l == K);
    assert(r == (N-K));
    return;
}

/// NOTE: not constant time
/// @param res
/// @param G
/// @param c
void apply_cf_action_to_G_with_pivots(generator_mat_t* res,
                                      const generator_mat_t *G,
                                      const uint8_t *const c,
                                      uint8_t initial_G_col_pivot[N],
                                      uint8_t permuted_G_col_pivot[N]) {
    uint32_t l = 0, r = 0;
    for (uint32_t i = 0; i < N8; i++) {
        for (uint32_t j = 0; j < 8; j++) {
            if ((i*8 + j) >= N) { goto finish; }

            const uint8_t bit = (c[i] >> j) & 1u;
            uint32_t pos;
            if (bit) {
                pos = l;
                l += 1;
            } else {
                pos = K + r;
                r += 1;
            }

            permuted_G_col_pivot[pos] = initial_G_col_pivot[i*8+j];

            // copy the column
            for (uint32_t k = 0; k < K; k++) {
                res->values[k][pos] = G->values[k][i*8 + j];
            }
        }
    }
    finish:
    assert(l == K);
    assert(r == (N-K));
    return;
}

/* Compresses a generator matrix in RREF storing only non-pivot columns and
 * their position */
void generator_rref_compact(rref_generator_mat_t *compact,
                            const generator_mat_t *const full,
                            const uint8_t is_pivot_column[N] )
{
   int dst_col_idx = 0;
   for (uint32_t src_col_idx = 0; src_col_idx < N; src_col_idx++) {
      if(!is_pivot_column[src_col_idx]) {
         for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
            compact->values[row_idx][dst_col_idx] = full->values[row_idx][src_col_idx];
         }
         compact->column_pos[dst_col_idx] = src_col_idx;
         dst_col_idx++;
      }
   }
} /* end generator_rref_compact */

void generator_to_normalized(normalized_IS_t *V,
                             const generator_mat_t *const G){
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < N - K; ++j) {
            V->values[i][j] = G->values[i][K+j];
        }
    }
}

/* Compresses a generator matrix in RREF into a array of bytes */
void compress_rref(uint8_t *compressed, const generator_mat_t *const full,
                   const uint8_t is_pivot_column[N]) {
    // Compress pivot flags
    for (uint32_t col_byte = 0; col_byte < N / 8; col_byte++) {
        compressed[col_byte] = is_pivot_column[8 * col_byte + 0] |
                               (is_pivot_column[8 * col_byte + 1] << 1) |
                               (is_pivot_column[8 * col_byte + 2] << 2) |
                               (is_pivot_column[8 * col_byte + 3] << 3) |
                               (is_pivot_column[8 * col_byte + 4] << 4) |
                               (is_pivot_column[8 * col_byte + 5] << 5) |
                               (is_pivot_column[8 * col_byte + 6] << 6) |
                               (is_pivot_column[8 * col_byte + 7] << 7);
    }

#if defined(CATEGORY_1) || defined(CATEGORY_5)
    // Compress last flags
    compressed[N / 8] = is_pivot_column[N - 4] | (is_pivot_column[N - 3] << 1) |
                        (is_pivot_column[N - 2] << 2) |
                        (is_pivot_column[N - 1] << 3);

    int compress_idx = N / 8 + 1;
#else
    int compress_idx = N / 8;
#endif

    // Compress non-pivot columns row-by-row
    int encode_state = 0;
    for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
        for (uint32_t col_idx = 0; col_idx < N; col_idx++) {
            if (!is_pivot_column[col_idx]) {
                switch (encode_state) {
                    case 0:
                        compressed[compress_idx] = full->values[row_idx][col_idx];
                        break;
                    case 1:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 7);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 1);
                        break;
                    case 2:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 6);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 2);
                        break;
                    case 3:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 5);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 3);
                        break;
                    case 4:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 4);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 4);
                        break;
                    case 5:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 3);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 5);
                        break;
                    case 6:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 2);
                        compress_idx++;
                        compressed[compress_idx] = (full->values[row_idx][col_idx] >> 6);
                        break;
                    case 7:
                        compressed[compress_idx] =
                                compressed[compress_idx] | (full->values[row_idx][col_idx] << 1);
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
    } /* end compress_rref */
}

/* Expands a compressed RREF generator matrix into a full one */
void expand_to_rref(generator_mat_t *full,
                    const uint8_t *compressed,
                    uint8_t is_pivot_column[N]) {
    /// TODO: this whole decompression code, can be removed. Specially since we moved to
    /// the `reusage of pivots` strategy, which makes this very inefficient.
    // Decompress pivot flags
    for (int i = 0; i < N; i++) {
        is_pivot_column[i] = 0;
    }

    for (int col_byte = 0; col_byte < N / 8; col_byte++) {
        is_pivot_column[col_byte * 8 + 0] = compressed[col_byte] & 0x1;
        is_pivot_column[col_byte * 8 + 1] = (compressed[col_byte] >> 1) & 0x1;
        is_pivot_column[col_byte * 8 + 2] = (compressed[col_byte] >> 2) & 0x1;
        is_pivot_column[col_byte * 8 + 3] = (compressed[col_byte] >> 3) & 0x1;
        is_pivot_column[col_byte * 8 + 4] = (compressed[col_byte] >> 4) & 0x1;
        is_pivot_column[col_byte * 8 + 5] = (compressed[col_byte] >> 5) & 0x1;
        is_pivot_column[col_byte * 8 + 6] = (compressed[col_byte] >> 6) & 0x1;
        is_pivot_column[col_byte * 8 + 7] = (compressed[col_byte] >> 7) & 0x1;
    }

#if defined(CATEGORY_1) || defined(CATEGORY_5)
    // Decompress last flags
    is_pivot_column[N - 4] = compressed[N / 8] & 0x1;
    is_pivot_column[N - 3] = (compressed[N / 8] >> 1) & 0x1;
    is_pivot_column[N - 2] = (compressed[N / 8] >> 2) & 0x1;
    is_pivot_column[N - 1] = (compressed[N / 8] >> 3) & 0x1;

    int compress_idx = N / 8 + 1;
#else
    int compress_idx = N / 8;
#endif

    // Decompress columns row-by-row
    int decode_state = 0;
    for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
        int pivot_idx = 0;
        for (uint32_t col_idx = 0; col_idx < N; col_idx++) {
            if (!is_pivot_column[col_idx]) {
                // Decompress non-pivot
                switch (decode_state) {
                    case 0:
                        full->values[row_idx][col_idx] = compressed[compress_idx] & MASK_Q;
                        break;
                    case 1:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 7) |
                                 (compressed[compress_idx + 1] << 1)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 2:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 6) |
                                 (compressed[compress_idx + 1] << 2)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 3:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 5) |
                                 (compressed[compress_idx + 1] << 3)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 4:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 4) |
                                 (compressed[compress_idx + 1] << 4)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 5:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 3) |
                                 (compressed[compress_idx + 1] << 5)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 6:
                        full->values[row_idx][col_idx] =
                                ((compressed[compress_idx] >> 2) |
                                 (compressed[compress_idx + 1] << 6)) &
                                MASK_Q;
                        compress_idx++;
                        break;
                    case 7:
                        full->values[row_idx][col_idx] =
                                (compressed[compress_idx] >> 1) & MASK_Q;
                        compress_idx++;
                        break;
                }

                if (decode_state != 7) {
                    decode_state++;
                } else {
                    decode_state = 0;
                }
            } else {
                // Decompress pivot
                full->values[row_idx][col_idx] = ((uint32_t)row_idx == (uint32_t)pivot_idx);
                pivot_idx++;
            }
        }
    }

} /* end expand_to_rref */


/* Expands a compressed RREF generator matrix into a full one */
void generator_rref_expand(generator_mat_t *full,
                           const rref_generator_mat_t *const compact)
{
   int placed_dense_cols = 0;
   for (uint32_t col_idx = 0; col_idx < N; col_idx++) {
      if ( (placed_dense_cols< N-K) &&
            (col_idx == compact->column_pos[placed_dense_cols])) {
         /* non-pivot column, restore one full column */
         for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
            full->values[row_idx][col_idx] = compact->values[row_idx][placed_dense_cols];
         }
         placed_dense_cols++;
      } else {
         /* regenerate the appropriate pivot column */
         for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
            full->values[row_idx][col_idx] = (row_idx == col_idx-placed_dense_cols);
         }
      }
   }
} /* end generator_rref_expand */

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
/// @param V [in/out]
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

// V1 =V2
void normalized_copy(normalized_IS_t *V1,
                     const normalized_IS_t *V2) {
    memcpy(V1->values, V2->values, sizeof(normalized_IS_t));
}

/// \param G: non IS part: G[row] *= a
/// \param row
/// \param a
void normalized_mat_scale_row(normalized_IS_t *G,
                              const uint32_t row,
                              const FQ_ELEM a) {
    for (uint32_t col = 0; col < N-K; col++) {
        G->values[row][col] = fq_mul(G->values[row][col], a);
    }
}

///
/// \param V
/// \param col
/// \return 0 if the column is zero
///         1 if the columns is non-zero
int normalized_is_zero_column(const normalized_IS_t *const V,
                              const uint32_t col) {
    for (uint32_t i = 0; i < K; i++) {
        if (V->values[i][col] > 0) {
            return 1;
        }
    }

    return 0;
}

///
/// \param V
/// \param col
/// \return 0 if every value in columns is non zerp
///         1 otherwise
int normalized_is_zero_in_column(const normalized_IS_t *const V,
                              const uint32_t col) {
    for (uint32_t i = 0; i < K; i++) {
        if (V->values[i][col] == 0) {
            return 1;
        }
    }

    return 0;
}


void generator_SF_seed_expand(rref_generator_mat_t *res,
                              const unsigned char seed[SEED_LENGTH_BYTES])
{
   SHAKE_STATE_STRUCT csprng_state;
   initialize_csprng(&csprng_state,seed,SEED_LENGTH_BYTES);
   for(uint32_t i = 0; i < K; i++) {
      rand_range_q_state_elements(&csprng_state, res->values[i], N-K);
   }
   for(uint32_t i = 0; i < N-K ; i++) {
      res->column_pos[i]=i+K;
   }


} /* end generator_seed_expand */

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
void generator_pretty_print_name(char *name, const generator_mat_t *const G)
{
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

void normalized_pretty_print(const normalized_IS_t *const G) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < (N-K-1); ++j) {
            printf("%3d,", G->values[i][j]);
        }
        printf("%3d\n", G->values[i][N-K-1]);
    }

    printf("\n");
}
void normalized_pretty_print_v(const FQ_ELEM values[K][N-K]) {
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < (N-K-1); ++j) {
            printf("%3d,", values[i][j]);
        }
        printf("%3d\n", values[i][N-K-1]);
    }

    printf("\n");
}
