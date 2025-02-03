/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include "sort.h"


/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// first sort the rows, then the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G,
                                 const uint8_t *L) {
    if (SortRows(G, K, L) == 0) {
	    return 0;
    }
    SortCols(G, K_pad);
    return 1;
}

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4(normalized_IS_t *G,
                                 const uint8_t *L) {
	for (uint32_t row = 0; row < K; row++) {
		if (row_all_same(G->values[row])) { continue; }
		FQ_ELEM s = row_acc(G->values[row]);

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = row_acc_inv(G->values[row]);
			if (s == 0) {
				return 0;
			}
		}

		row_mul(G->values[row], s);
	}

	return compute_canonical_form_type3(G, L);
}

/// NOTE: non-constant time
/// \param G[in/out]: sub matrix with only z rows.
///     Nothing really is returned in the matrix, but it still
///     gets clobbered.
/// \param z[in]: number of rows in G
/// \param M[in]: the currently shortest multiset
/// \return 0: if no multiset was found < `M`
///         1: if one of the  z rows is < `M`
int compute_canonical_form_type4_sub(normalized_IS_t *G,
                                     const uint32_t z,
                                     const FQ_ELEM *M) {
    FQ_ELEM tmp[Q_pad] __attribute__((aligned(32))) = {0};
    for (uint32_t i = 0; i < z; i++) {
		FQ_ELEM s = row_acc(G->values[i]);

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = row_acc_inv(G->values[i]);
			if (s == 0) { return 0; }
		}

		row_mul(G->values[i], s);
        sort(tmp, G->values[i], N-K);
        if (compare_rows(tmp, M) < 0) {
            return 1;
        }
    }

    return 0;
}

/// NOTE: non-constant time
/// implements a total order on matrices
/// we simply compare the columns lexicographically
/// \param V1[in]: first matrix
/// \param V2[in]: second matrix
/// \param z[in]: number of rows within both matrices
/// \return -x if V2 > V1
///			 0 if V2 == V1
///			+x
int compare_matrices(const normalized_IS_t *__restrict__ V1,
                     const normalized_IS_t *__restrict__ V2,
                     const uint32_t z) {
	for (uint32_t row = 0; row < z; row++) {
        uint32_t i=0;
        while((i < N-K) &&
             ((V1->values[row][i] - V2->values[row][i]) == 0)) {
            i++;
        }

        if (i >= K) { continue; }

        return (int)(V1->values[row][i]) - (int)(V2->values[row][i]);
	}
	
	// if we are here the two matrices are equal
	return 0;
}

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5(normalized_IS_t *G) {
	normalized_IS_t A __attribute__((aligned(32))) = {0}, M;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&M.values, Q-1, sizeof(normalized_IS_t));

	static FQ_ELEM row_inv_data[N_K_pad] = {0};
	for (uint32_t row = 0; row < K; row++) {
        if (row_contains_zero(G->values[row])) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < K; row2++) {
			row_mul3(A.values[row2], G->values[row2], row_inv_data);
		}

		const int ret = compute_canonical_form_type4(&A, NULL);
		if ((ret == 1) && (compare_matrices(&A, &M, K) < 0)) {
			touched = 1;
			normalized_copy(&M, &A);
		}
	}

    if (!touched) { return 0; }
	normalized_copy(G, &M);
	return 1;
}

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_popcnt(normalized_IS_t *G) {
	normalized_IS_t M;
    int touched = 0;

	// init the output matrix to some `invalid` data
    // honestly, will be there ever a case, where more than 2 rows are compared?
	memset(&M, Q-1, sizeof(normalized_IS_t));

    /// track the rows with the most zeros.
	uint32_t J[K];
    uint32_t z = 0;

    // count zeros in each row
    uint32_t max_zeros = 0;
	uint8_t Z[K] = {0};
	for (uint32_t row = 0; row < K; row++) {
        const uint32_t num_zeros = row_count_zero(G->values[row]);

        if (num_zeros > 0) {
            Z[row] = 1;
        }

        if (num_zeros > max_zeros) {
            z = 1;
            J[0] = row;
            max_zeros = num_zeros;
        	continue;
        }

        if (num_zeros == max_zeros) {
            J[z++] = row;
        }
    }

    /// NOTE: fallback solution if everything falls apart
    if (z == (N-K)) {
	    return compute_canonical_form_type5(G);
    }

    static normalized_IS_t B __attribute__((aligned(32))) = {0};
	FQ_ELEM row_inv_data[N_K_pad] = {0};

	/// NOTE: this is already "sorted"
	FQ_ELEM L[Q_pad] __attribute__((aligned(32))) = {0};
	for (uint32_t row = 0; row < K; row++) {
        if (Z[row]) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < z; row2++) {
			row_mul3(B.values[row2], G->values[J[row2]], row_inv_data);
		}

        if (compute_canonical_form_type4_sub(&B, z, L)) {
            for (uint32_t row2 = 0; row2 < K; row2++) {
                row_mul3(B.values[row2], G->values[row2], row_inv_data);
            }

		    const int ret = compute_canonical_form_type4(&B, L);
		    if ((ret == 1) && (compare_matrices(&B, &M, K) < 0)) {
                sort(L, B.values[0], N-K);
		    	touched = 1;
		    	normalized_copy(&M, &B);
		    }
        }
	}

	normalized_copy(G, &M);
    return touched;
}

/// samples to random monomial matrices (A, B) and comptes A*G*B
/// \param G[in/out] non IS part of a generator matrix
/// \param prng[in/out]:
void blind(normalized_IS_t *G,
           SHAKE_STATE_STRUCT *prng) {
    /// NOTE: the type `monomial` allocates a N elements, which is a little
    /// bit a overkill, as we only need K or N-K.
    monomial_t left, right;

    // NOTE: init with zero to please valgrind, as the tailing (K_pad - K) rows
    // are never used/initialized.
    static normalized_IS_t B = {0};

    for(uint32_t i = 0; i < N-K; i++) {
        left.permutation[i] = i;
        right.permutation[i] = i;
    }

    // We compute the following matrix multiplication G = left * G * right
    // where `left` and `right` are randomly sampled monomials
    fq_star_rnd_state_elements(prng, right.coefficients, N-K);
    yt_shuffle_state_limit(prng, right.permutation, N-K);

    fq_star_rnd_state_elements(prng, left.coefficients, N-K);
    yt_shuffle_state_limit(prng, left.permutation, N-K);

    // apply the right multiplication
    for (uint32_t i = 0; i < K; i++) {
        const FQ_ELEM a = right.coefficients[i];
        const POSITION_T pos = right.permutation[i];

        /// NOTE: thats quite a bottleneck.
        for (uint32_t j = 0; j < N-K; j++) {
            B.values[j][i] = fq_mul(a,  G->values[j][pos]);
        }
    }

    // apply the left multiplication
    for (uint32_t i = 0; i < K; i++) {
        const FQ_ELEM a = left.coefficients[i];
        const POSITION_T pos = left.permutation[i];

		row_mul2(G->values[i], B.values[pos], a);
    }
}

/// NOTE: non-constant time implementation
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int CF(normalized_IS_t *G) {
    return compute_canonical_form_type5_popcnt(G);
}
