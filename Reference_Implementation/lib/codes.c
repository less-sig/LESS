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

#include <string.h>
#include <stdio.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include <assert.h>

/// swap N uint8 in r and s.
/// \param r[in/out]
/// \param s[in/out]
void swap_rows(FQ_ELEM r[N],
               FQ_ELEM s[N]) {
    FQ_ELEM tmp[N];
    memcpy(tmp, r, sizeof(FQ_ELEM) * N);
    memcpy(r, s, sizeof(FQ_ELEM) * N);
    memcpy(s, tmp, sizeof(FQ_ELEM) * N);
} /* end swap_rows */

/* Calculate pivot flag array */
void generator_get_pivot_flags(const rref_generator_mat_t *const G,
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
                            const monomial_t *const monom) {
   for(uint32_t src_col_idx = 0; src_col_idx < N; src_col_idx++) {
      for(uint32_t row_idx = 0; row_idx < K; row_idx++) {
         res->values[row_idx][monom->permutation[src_col_idx]] =
            fq_mul(G->values[row_idx][src_col_idx], monom->coefficients[src_col_idx]);
      }
   }
} /* end generator_monomial_mul */

/// \param G[in/out]: generator matrix
/// \param is_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \return 0 on failure
///         1 on success

int generator_RREF(generator_mat_t *G,
                   uint8_t is_pivot_column[N]) {
   for(unsigned row_to_reduce = 0; row_to_reduce < K; row_to_reduce++) {
      unsigned pivot_row = row_to_reduce;
      /*start by searching the pivot in the col = row*/
      unsigned pivot_column = row_to_reduce;
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
      FQ_ELEM scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);

      /* rescale pivot row to have pivot = 1. Values at the left of the pivot
       * are already set to zero by previous iterations */
      for(unsigned i = pivot_column; i < N; i++) {
         G->values[pivot_row][i] = fq_mul(scaling_factor, G->values[pivot_row][i]);
      }

      /* Subtract the now placed and reduced pivot rows, from the others,
       * after rescaling it */
      for(unsigned row_idx = 0; row_idx < K; row_idx++) {
         if (row_idx != pivot_row) {
            FQ_ELEM multiplier = G->values[row_idx][pivot_column];
            /* all elements before the pivot in the pivot row are null, no need to
             * subtract them from other rows. */
            for(unsigned col_idx = 0; col_idx < N; col_idx++) {
               FQ_ELEM tmp = fq_mul(multiplier, G->values[pivot_row][col_idx]);
               G->values[row_idx][col_idx] = fq_sub(G->values[row_idx][col_idx], tmp);
            }
         }
      }
   }

   return 1;
} /* end generator_RREF */

/// \param G[in/out]: generator matrix K \times N
/// \param is_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param was_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param pvt_reuse_limit:[in]:
/// \return 0 on failure
///         1 on success
int generator_RREF_pivot_reuse(generator_mat_t *G,
                               uint8_t is_pivot_column[N],
                               uint8_t was_pivot_column[N],
                               const int pvt_reuse_limit) {
   int pvt_reuse_cnt = 0;

    // row swap pre-process - swap previous pivot elements to corresponding row to reduce likelihood of corruption
    if (pvt_reuse_limit != 0) {
        for (int preproc_col = K - 1; preproc_col >= 0; preproc_col--) {
            if (was_pivot_column[preproc_col] == 1) {
                // find pivot row
                uint32_t pivot_el_row = -1;
                for (uint32_t row = 0; row < K; row = row + 1) {
                    if (G->values[row][preproc_col] != 0) {
                        pivot_el_row = row;
                    }
                }
                swap_rows(G->values[preproc_col], G->values[pivot_el_row]);
            }
        }
    }

    for (uint32_t row_to_reduce = 0; row_to_reduce < K; row_to_reduce++) {
        uint32_t pivot_row = row_to_reduce;
        /*start by searching the pivot in the col = row*/
        uint32_t pivot_column = row_to_reduce;
        while ((pivot_column < N) && (G->values[pivot_row][pivot_column] == 0)) {
            while ((pivot_row < K) && (G->values[pivot_row][pivot_column] == 0)) {
                pivot_row++;
            }
            if (pivot_row >= K) { /*entire column tail swept*/
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
            was_pivot_column[pivot_row] = 0; // pivot no longer reusable - will be corrupted during reduce row
            swap_rows(G->values[row_to_reduce], G->values[pivot_row]);
        }
        pivot_row = row_to_reduce; /* row with pivot now in place */

        /// NOTE: this needs explenation. We can skip the reduction of the pivot row, because for
        /// the CF it doesnt matter. The only thing that is important for the CF is the number of
        /// zeros, and this doest change if we reduce a reused pivot row.
        if (((was_pivot_column[pivot_column] == 1) && (pvt_reuse_cnt < pvt_reuse_limit) && (pivot_column < K))) {
            pvt_reuse_cnt++;
            continue;
        }

        /* Compute rescaling factor */
        const FQ_ELEM scaling_factor = fq_inv(G->values[pivot_row][pivot_column]);

        /* rescale pivot row to have pivot = 1. Values at the left of the pivot
         * are already set to zero by previous iterations */
        for (uint32_t i = pivot_column; i < N; i++) {
            G->values[pivot_row][i] = fq_mul(scaling_factor, G->values[pivot_row][i]);
        }

        /* Subtract the now placed and reduced pivot rows, from the others,
         * after rescaling it */
        for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
            if (row_idx != pivot_row) {
                FQ_ELEM multiplier = G->values[row_idx][pivot_column];
                /* all elements before the pivot in the pivot row are null, no need to
                 * subtract them from other rows. */
                for (int col_idx = 0; col_idx < N; col_idx++) {
                    FQ_ELEM tmp = fq_mul(multiplier, G->values[pivot_row][col_idx]);
                    G->values[row_idx][col_idx] = fq_sub(G->values[row_idx][col_idx], tmp);
                }
            }
        }
    }

    return 1;
} /* end generator_RREF_pivot_reuse */

/// NOTE: not constant time
/// \param res[out]: G*c a generator matrix: K \times N-K
/// \param G[in]: current generator matrix: K \times N-K
/// \param c[in]: compressed cf action
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
    return;
}

/// NOTE: not constant time
/// \param res[out]: G*c a generator matrix: K \times N-K
/// \param G[in]: current generator matrix: K \times N-K
/// \param c[in]: compressed cf action
void apply_cf_action_to_G_with_pivots(generator_mat_t* res,
                                      const generator_mat_t *G,
                                      const uint8_t *const c,
                                      const uint8_t initial_G_col_pivot[N],
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
    return;
}

/* Compresses a generator matrix in RREF storing only non-pivot columns and
 * their position */
void generator_rref_compact(rref_generator_mat_t *compact,
                            const generator_mat_t *const full,
                            const uint8_t is_pivot_column[N]) {
    int dst_col_idx = 0;
    for (uint32_t src_col_idx = 0; src_col_idx < N; src_col_idx++) {
        if (!is_pivot_column[src_col_idx]) {
            for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
                compact->values[row_idx][dst_col_idx] = full->values[row_idx][src_col_idx];
            }
            compact->column_pos[dst_col_idx] = src_col_idx;
            dst_col_idx++;
        }
    }
} /* end generator_rref_compact */

/* Compresses a generator matrix in RREF into an array of bytes */
void compress_rref(uint8_t *compressed,
                   const generator_mat_t *const full,
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

#if (CATEGORY == 252) || (CATEGORY == 548)
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

#if (CATEGORY == 252) || (CATEGORY == 548)
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
                           const rref_generator_mat_t *const compact) {
    int placed_dense_cols = 0;
    for (uint32_t col_idx = 0; col_idx < N; col_idx++) {
        if ((placed_dense_cols < N - K) && (col_idx == compact->column_pos[placed_dense_cols])) {
            /* non-pivot column, restore one full column */
            for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
                full->values[row_idx][col_idx] = compact->values[row_idx][placed_dense_cols];
            }
            placed_dense_cols++;
        } else {
            /* regenerate the appropriate pivot column */
            for (uint32_t row_idx = 0; row_idx < K; row_idx++) {
                full->values[row_idx][col_idx] = (row_idx == col_idx - placed_dense_cols);

            }
        }
    }
} /* end generator_rref_expand */

// V1 =V2
void normalized_copy(normalized_IS_t *V1,
                     const normalized_IS_t *V2) {
    memcpy(V1->values, V2->values, sizeof(normalized_IS_t));
}

/// \param V[in/out]: K \times N-K matrix in which row `row1` and
///                 row `row2` are swapped
/// \param row1[in]: first row
/// \param row2[in]: second row
void normalized_row_swap(normalized_IS_t *V,
                         const POSITION_T row1,
                         const POSITION_T row2) {
    if (row1 == row2) { return; }
    for(uint32_t i = 0; i < N-K; i++){
        POSITION_T tmp = V->values[row1][i];
        V->values[row1][i] = V->values[row2][i];
        V->values[row2][i] = tmp;
    }
}

/// \param res[out]: full rank generator matrix K \times N-K
/// \param seed[int] seed for the prng
void generator_sample(rref_generator_mat_t *res, const unsigned char seed[SEED_LENGTH_BYTES]) {
    SHAKE_STATE_STRUCT csprng_state;
    initialize_csprng(&csprng_state, seed, SEED_LENGTH_BYTES);
    for (uint32_t i = 0; i < K; i++) {
        rand_range_q_state_elements(&csprng_state, res->values[i], N - K);
    }
    for (uint32_t i = 0; i < N - K; i++) {
        res->column_pos[i] = i + K;
    }
} /* end generator_seed_expand */
