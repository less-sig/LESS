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

#pragma once

#include "parameters.h"
#include "monomial_mat.h"
#include <stdint.h>

typedef struct {  /* Generator matrix, stored explicitly */
   FQ_ELEM values[K][N];
} generator_mat_t;

/* RREF Generator mat., only values and positions of non-pivot columns stored */
typedef struct {
   FQ_ELEM values[K][N-K];     /* values of the non-pivot columns    */
   POSITION_T column_pos[N-K]; /* positions of the non-pivot columns */
} rref_generator_mat_t;

/* Set of columns not constituting the IS for an RREF matrix
 * See algorithm PrepareDigestInput in specification (V matrix)*/
typedef struct {
   FQ_ELEM values[K][N-K];   /* values of the non-pivot columns */
} normalized_IS_t;


/* multiplies a monomial matrix by a generator matrix */
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *const G,
                            const monomial_t *const monom);


/** Computes the row-reduced echelon form of the generator matrix
 *  returns 1 on success, 0 on failure, computation is done in-place
 *  Provides the positions of the pivot columns, one-hot encoded in
 *  is_pivot_column
 **/
int generator_RREF(generator_mat_t *G,
                   uint8_t is_pivot_column[N]);


void lex_minimize(normalized_IS_t *V,
                  POSITION_T dst_col_idx,
                  const generator_mat_t *const G,
                  const POSITION_T col_idx);

/* performs lexicographic sorting of the IS complement */
void lex_sort_cols(normalized_IS_t *V);

void prepare_digest_input(normalized_IS_t *V,
                          monomial_action_IS_t *Q_bar_IS,
                          const generator_mat_t *const G,
                          const monomial_t *const Q_in);

/* extracts the last N-K columns from a generator matrix, filling
 * in the compact RREF representation*/
void generator_rref_compact(rref_generator_mat_t *compact,
                            const generator_mat_t *const full,
                            const uint8_t is_pivot_column[N] );

/* Compresses a columns of an IS matrix */
void compress_columns(uint8_t *compressed,
                    const normalized_IS_t *const full);

/* Compresses a generator matrix in RREF into a array of bytes */
void compress_rref(uint8_t *compressed,
                    const generator_mat_t *const full,
                    const uint8_t is_pivot_column[N]);

/* Expands a compressed RREF generator matrix into a full one */
void expand_to_rref(generator_mat_t *full,
                        const uint8_t *compressed);

/* Takes as input a compact RREF generator matrix, i.e. a set of N-K
 * columns and their position in the RREF and normalizes the columns themselves
 * i.e., rescales them to obtain the leading nonzero term equal to 1, and
 * sorts them according to the lexicographic order.
 * The action is done in-place; the column_pos vector is not altered so that
 * the normalized, compacted form can be reexpanded (albeit not needed
 * in the current protocol)
 */
void generator_compact_rref_normalize(rref_generator_mat_t *compact);

/* Expands a compact representation of a generator matrix into full form*/
void generator_rref_expand(generator_mat_t *full,
                           const rref_generator_mat_t *const compact);

/* Computes the systematic form of the generator matrix in place
 * returns 1 on success, 0 on failure */
int generator_gausselim(generator_mat_t *G);

void apply_action_to_G(generator_mat_t* res,
                       const generator_mat_t* G,
                       const monomial_action_IS_t* Q_IS);

/* samples a random monomial matrix */
void generator_rnd(generator_mat_t *res);


/* expands a systematic form generator from a seed randomly drawing only
 * non-identity portion */
void generator_SF_seed_expand(rref_generator_mat_t *res,
                              const unsigned char seed[SEED_LENGTH_BYTES]);

void generator_pretty_print_name(char *name, const generator_mat_t *const G);
void generator_rref_pretty_print_name(char *name,
                                      const rref_generator_mat_t *const G);
