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

#pragma once

#include <stdint.h>

#include "parameters.h"
#include "monomial_mat.h"

/// Generator matrix, stored explicitly
typedef struct {  
    // NOTE: the alignment is needed for the optimized AVX{2|512} implementation
    FQ_ELEM values[K][N_pad] __attribute__((aligned(32)));
} generator_mat_t;

/// RREF Generator mat., only values and positions of non-pivot columns stored
typedef struct {
    /// values of the non-pivot columns
    FQ_ELEM values[K][N_K_pad];
    // positions of the non-pivot columns
    POSITION_T column_pos[N_K_pad];
} rref_generator_mat_t;

/// Set of columns not constituting the IS for an RREF matrix.
typedef struct {
    /// values of the non-pivot columns
    FQ_ELEM values[K_pad][N_K_pad];   
} normalized_IS_t;

/// Calculate pivot flag array
/// \param G[in]:
/// \param pivot_flag[ou]:
void generator_get_pivot_flags(const rref_generator_mat_t *const G,
                               uint8_t pivot_flag [N]);

/// multiplies a monomial matrix by a generator matrix
/// \param res[out]: pointer to an uninitialized generator matrix
/// \param G[in]: full (K \times N) generator matrix
/// \param monom[in]: (random) monomial matrix
void generator_monomial_mul(generator_mat_t *res,
                            const generator_mat_t *const G,
                            const monomial_t *const monom);


/// Computes the row-reduced echelon form of the generator matrix
/// returns 1 on success, 0 on failure, computation is done in-place
/// Provides the positions of the pivot columns, one-hot encoded in
/// is_pivot_column
/// \param G[in/out] generator matrix (K \times N)
/// \param is_pivot_column[out]: position of the pivot columns, indicated 
///     by a bit flag.
int generator_RREF(generator_mat_t *G,
                   uint8_t is_pivot_column[N_pad]);

/// \param G[in/out]: generator matrix K \times N
/// \param is_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param was_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param pvt_reuse_limit:[in]:
/// \return 0 on failure
///         1 on success
int generator_RREF_pivot_reuse_ct(generator_mat_t *G,
                                  uint8_t is_pivot_column[N],
                                  uint8_t was_pivot_column[N],
                                  const int pvt_reuse_limit);

/// NOTE: not constant time
/// \param G[in/out]: generator matrix K \times N
/// \param is_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param was_pivot_column[out]: N bytes, set to 1 if this column
///                 is a pivot column
/// \param pvt_reuse_limit:[in]: max number of pivots to reuse
/// \return 0 on failure
///         1 on success
int generator_RREF_pivot_reuse(generator_mat_t *G,
                               uint8_t is_pivot_column[N],
                               uint8_t was_pivot_column[N],
                               const int pvt_reuse_limit);

/// Compresses a generator matrix in RREF into a array of bytes
/// \param compressed[out] byte array of length RREF_MAT_PACKEDBYTES
/// \param full[in]: full generator matrix (K \times N)
/// \param is_pivot_column[in]: array of length N in which K fields are 1, the 
///     rest must be zero. Indicating the positions of the pivot columns.
void compress_rref(uint8_t *compressed,
                   const generator_mat_t *const full,
                   const uint8_t is_pivot_column[N]);

/// Expands a compressed RREF generator matrix into a full one
/// \param full[out]: output full matrix (K \times N)
/// \param compressed[in]: bytestream containing the compressed maitrx
/// \param is_pivot_column[out]: N bytes will be initialized with zeros. And 
///     only 1 will be written at the column position which is a pivot column.
void expand_to_rref(generator_mat_t *full,
                    const uint8_t *compressed,
                    uint8_t is_pivot_column[N]);

/// Expands a compressed RREF generator matrix into a full one
/// \param full[out]: output generator matrix (K \times N) 
/// \param compact[out]: input compressed generator matrix (K \times N-K) 
void generator_rref_expand(generator_mat_t *full,
                           const rref_generator_mat_t *const compact);

/// expands a systematic form generator from a seed randomly drawing only
/// non-identity portion
/// \param res[out]: full rank generator matrix K \times N-K
/// \param seed[int] seed for the prng
void generator_sample(rref_generator_mat_t *res,
                      const unsigned char seed[SEED_LENGTH_BYTES]);

/// \param res[out]: G*c a generator matrix: K \times N-K
/// \param G[in]: current generator matrix: K \times N-K
/// \param c[in]: compressed cf action
/// \param initial_G_col_pivot[in]: input IS
/// \param permuted_G_col_pivot[out]: output IS, to keep track of the pivot cols
void apply_cf_action_to_G_with_pivots(generator_mat_t* res,
                                      const generator_mat_t *G,
                                      const uint8_t *const c,
                                      const uint8_t initial_G_col_pivot[N],
                                      uint8_t permuted_G_col_pivot[N]);

/// V1 = V2 
/// \param V1[out]: pointer to generator matrix (non IS part)
/// \param V2[in]: pointer to generator matrix (non IS part)
void normalized_copy(normalized_IS_t *V1,
                     const normalized_IS_t *V2);

/// \param V[in/out]: K \times N-K matrix in which row `row1` and
///     row `row2` are swapped
/// \param row1[in]: first row
/// \param row2[in]: second row
void normalized_row_swap(normalized_IS_t *V,
                         const POSITION_T row1,
                         const POSITION_T row2);

/// right-multiplies a generator by a monomial: res = G*monom
/// \param res[out] pointer to an uninitialized generator matrix (non IS part)
/// \param G[in]: pointer to an initialized generator matrix (non IS part)
/// \param monom[in]: pointer to an initialized monomial matrix
void normalized_monomial_right(normalized_IS_t *res,
                               const normalized_IS_t *const G,
                               const monomial_t *const monom);
