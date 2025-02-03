/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (May 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Duc Tri Nguyen <dnguye69@gmu.edu>
 *
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
#include <stdint.h>

#include "codes.h"
#include "fq_arith.h"
#include "macro.h"
#include "utils.h"
#include "parameters.h"

// Select low 8-bit, skip the high 8-bit in 16 bit type
const uint8_t shuff_low_half[32] = {
        0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe,
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
        0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe,
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
};

/* Calculate pivot flag array */
void generator_get_pivot_flags (const rref_generator_mat_t *const G, uint8_t pivot_flag [N]) {
    for (int i = 0; i < N; i = i + 1) {
        pivot_flag[i] = 1;
    }

    for (int i = 0; i < K; i = i + 1) {
        pivot_flag[G->column_pos[i]] = 0;
    }
}

/// NOTE: not constant time
/// \param res
/// \param G
/// \param c
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
/// \param res
/// \param G
/// \param c
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

/* right-multiplies a generator by a monomial */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *const G,
                            const monomial_t *const monom) {

    FQ_ELEM buffer[N_pad];
    // Total SIMD registers: 9 = 2 + 3 + 4
    vec256_t shuffle, t, c8_127, c8_1;     // 2
    vec256_t g, mono, s;                   // 3
    vec256_t g_lo, g_hi, mono_lo, mono_hi; // 4
    vec128_t tmp;

    vload256(shuffle, (vec256_t *) shuff_low_half);
    vset8(c8_127, 127);
    vset8(c8_1, 1);
    for (uint32_t row_idx = 0; row_idx < K; row_idx++) {

        const FQ_ELEM *G_pointer = G->values[row_idx];
        for (uint32_t src_col_idx = 0; (src_col_idx+32) <= N_pad; src_col_idx += 32) {
            vload256(g, (vec256_t *) &G_pointer[src_col_idx]);
            vload256(mono, (vec256_t *) &(monom->coefficients[src_col_idx]));

            vget_lo(tmp, g);
            vextend8_16(g_lo, tmp);
            vget_hi(tmp, g);
            vextend8_16(g_hi, tmp);
            vget_lo(tmp, mono);
            vextend8_16(mono_lo, tmp);
            vget_hi(tmp, mono);
            vextend8_16(mono_hi, tmp);

            barrett_mul_u16(g_lo, g_lo, mono_lo, s);
            barrett_mul_u16(g_hi, g_hi, mono_hi, t);

            vshuffle8(g_lo, g_lo, shuffle);
            vshuffle8(g_hi, g_hi, shuffle);

            vpermute_4x64(g_lo, g_lo, 0xd8);
            vpermute_4x64(g_hi, g_hi, 0xd8);

            vpermute2(s, g_lo, g_hi, 0x20);

            barrett_red8(s, t, c8_127, c8_1);
            vstore256((vec256_t *) &buffer[src_col_idx], s);
        }

        for (uint32_t i = 0; i < N; i++) {
            res->values[row_idx][monom->permutation[i]] = buffer[i];
        }
    }
} /* end generator_monomial_mul */



void swap_rows(FQ_ELEM r[N_pad], 
               FQ_ELEM s[N_pad]){
    vec256_t a, b;
    for(uint32_t i=0; i<N_pad; i+=32) {
         vload256(a, (const vec256_t *)(r + i));
         vload256(b, (const vec256_t *)(s + i));
         vstore256((vec256_t *)(r + i), b);
         vstore256((vec256_t *)(s + i), a);
    }
} /* end swap_rows */

int generator_RREF(generator_mat_t *G, uint8_t is_pivot_column[N_pad]) {
    int i, j, pivc;
    uint8_t tmp, sc;

    vec256_t *gm[K] __attribute__((aligned(32)));
    vec256_t em[0x80][NW];
    vec256_t *ep[0x80];
    vec256_t c01, c7f;
    vec256_t x, t, *rp, *rg;

    vset8(c01, 0x01);
    vset8(c7f, 0x7f);

    for (i = 0; i < K; i++) {
        gm[i] = (vec256_t *) G->values[i];
    }

    for (i = 0; i < K; i++) {
        j = i;
        /*start by searching the pivot in the col = row*/
        pivc = i;

        while (pivc < N) {
            while (j < K) {
                sc = G->values[j][pivc];
                if (sc != 0) {
                    goto found;
                }

                j++;
            }
            pivc++;     /* move to next col */
            j = i;      /*starting from row to red */
        }

        if (pivc >= N) {
            /* no pivot candidates left, report failure */
            return 0;
        }

        found:
        is_pivot_column[pivc] = 1; /* pivot found, mark the column*/

        /* if we found the pivot on a row which has an index > pivot_column
         * we need to swap the rows */
        if (i != j) {
            for (uint32_t k = 0; k < NW; k++) {
                t = gm[i][k];
                gm[i][k] = gm[j][k];
                gm[j][k] = t;
            }
        }

        /* Compute rescaling factor */
        /* rescale pivot row to have pivot = 1. Values at the left of the pivot
         * are already set to zero by previous iterations */

        //	generate the em matrix
        rg = gm[i];
        memcpy(em[1], rg, LESS_WSZ * NW);

        for (j = 2; j < 127; j++) {
            for (uint32_t k = 0; k < NW; k++) {
                vadd8(x, em[j - 1][k], rg[k])
                W_RED127_(x);
                em[j][k] = x;
            }
        }

        //	shuffle the pointers into ep
        sc = ((uint8_t *) rg)[pivc];
        ep[0] = em[0];
        tmp = sc;
        for (j = 1; j < 127; j++) {
            ep[tmp] = em[j];
            tmp += sc;
            tmp = (tmp + (tmp >> 7)) & 0x7F;
        }
        ep[0x7F] = em[0];

        //	copy back the normalized one
        memcpy(rg, ep[1], LESS_WSZ * NW);

        /* Subtract the now placed and reduced pivot rows, from the others,
         * after rescaling it */
        for (j = 0; j < K; j++) {
            sc = ((uint8_t *) gm[j])[pivc] & 0x7F;
            if (j != i && (sc != 0x00)) {
                rp = ep[127 - sc];
                for (uint32_t k = 0; k < NW; k++) {
                    vadd8(x, gm[j][k], rp[k])
                    W_RED127_(x);
                    gm[j][k] = x;
                }
            }
        }

    }

    return 1;
} /* end generator_RREF */

int generator_RREF_pivot_reuse(generator_mat_t *G,
                               uint8_t is_pivot_column[N],
                               uint8_t was_pivot_column[N],
                               const int pvt_reuse_limit) {
    int i, j, pivc;
    uint8_t tmp, sc;

    vec256_t *gm[K] __attribute__((aligned(32)));
    vec256_t em[0x80][NW];
    vec256_t *ep[0x80];
    vec256_t c01, c7f;
    vec256_t x, t, *rp, *rg;

    vset8(c01, 0x01);
    vset8(c7f, 0x7f);
    int pvt_reuse_cnt = 0;

    // this loop roughly takes 2.2% of the whole function runtime ()
    if (pvt_reuse_limit != 0) {
        for (int preproc_col = K - 1; preproc_col >= 0; preproc_col--) {
            if (was_pivot_column[preproc_col] == 1) {
                // find pivot row
                uint32_t pivot_el_row = 0;
                for (uint32_t row = 0; row < K; row = row + 1) {
                    if (G->values[row][preproc_col] != 0) {
                        pivot_el_row = row;
                    }
                }

                swap_rows(G->values[preproc_col], G->values[pivot_el_row]);
            }
        }
    }

    for (i = 0; i < K; i++) {
        gm[i] = (vec256_t *) G->values[i];
    }

    for (i = 0; i < K; i++) {
        j = i;
        /*start by searching the pivot in the col = row*/
        pivc = i;

        while (pivc < N) {
            while (j < K) {
                sc = G->values[j][pivc];
                if (sc != 0) {
                    goto found;
                }

                j++;
            }
            pivc++;     /* move to next col */
            j = i;      /*starting from row to red */
        }

        if (pivc >= N) {
            return 0; /* no pivot candidates left, report failure */
        }

        found:
        is_pivot_column[pivc] = 1; /* pivot found, mark the column*/

        /* if we found the pivot on a row which has an index > pivot_column
         * we need to swap the rows */
        if (i != j) {
            was_pivot_column[j] = 0; // pivot no longer reusable - will be corrupted during reduce row
            for (uint32_t k = 0; k < NW; k++) {
                t = gm[i][k];
                gm[i][k] = gm[j][k];
                gm[j][k] = t;
            }
        }

        /// NOTE: this needs explenation. We can skip the reduction of the pivot row, because for
        /// the CF it doesnt matter. The only thing that is important for the CF is the number of
        /// zeros, and this doest change if we reduce a reused pivot row.
        if ((was_pivot_column[pivc] == 1) &&
            (pvt_reuse_cnt < pvt_reuse_limit) &&
            (pivc < K)) {
            continue;
        }

        /* Compute rescaling factor */
        /* rescale pivot row to have pivot = 1. Values at the left of the pivot
         * are already set to zero by previous iterations */

        //  generate the em matrix
        rg = gm[i];
        memcpy(em[1], rg, LESS_WSZ * NW);

        for (j = 2; j < 127; j++) {
            for (uint32_t k = 0; k < NW; k++) {
                vadd8(x, em[j - 1][k], rg[k])
                W_RED127_(x);
                em[j][k] = x;
            }
        }

        //  shuffle the pointers into ep
        sc = ((uint8_t *) rg)[pivc];
        ep[0] = em[0];
        tmp = sc;
        for (j = 1; j < 127; j++) {
            ep[tmp] = em[j];
            tmp += sc;
            tmp = (tmp + (tmp >> 7)) & 0x7F;
        }
        ep[0x7F] = em[0];

        //  copy back the normalized one
        memcpy(rg, ep[1], LESS_WSZ * NW);

        /* Subtract the now placed and reduced pivot rows, from the others,
         * after rescaling it */
        for (j = 0; j < K; j++) {
            sc = ((uint8_t *) gm[j])[pivc];
            if (sc != 0x00 && j != i) {
                rp = ep[127 - sc];
                for (uint32_t k = 0; k < NW; k++) {
                    vadd8(x, gm[j][k], rp[k])
                    W_RED127_(x);
                    gm[j][k] = x;
                }
            }
        }
    }

    return 1;
} /* end generator_RREF */

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


// V1 =V2
void normalized_copy(normalized_IS_t *V1,
                     const normalized_IS_t *V2) {
    memcpy(V1->values, V2->values, sizeof(normalized_IS_t));
}

/// \param V
/// \param row1
/// \param row2
void normalized_row_swap(normalized_IS_t *V,
              const POSITION_T row1,
              const POSITION_T row2) {
    for(uint32_t i = 0; i < N-K; i++){
        POSITION_T tmp;
        tmp = V->values[row1][i];
        V->values[row1][i] = V->values[row2][i];
        V->values[row2][i] = tmp;
    }
}

/* Expands a compressed RREF generator matrix into a full one */
void generator_rref_expand(generator_mat_t *full,
                           const rref_generator_mat_t *const compact) {
   uint32_t placed_dense_cols = 0;
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

void generator_sample(rref_generator_mat_t *res,
                      const unsigned char seed[SEED_LENGTH_BYTES]) {
   SHAKE_STATE_STRUCT csprng_state;
   initialize_csprng(&csprng_state,seed,SEED_LENGTH_BYTES);
   for(uint32_t i = 0; i < K; i++) {
      rand_range_q_state_elements(&csprng_state, res->values[i], N-K);
   }
   for(uint32_t i = 0; i < N-K ; i++) {
      res->column_pos[i]=i+K;
   }


} /* end generator_seed_expand */
