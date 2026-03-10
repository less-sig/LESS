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
#include <string.h>

#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include "sort.h"

/// Implementation of the `SortCF` algorithm from the specification.
/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// NOTE: early exist from the `SortRows` function if `L` is already shorter.
/// NOTE: first sort the rows, then the columns
/// \param G[in/out]: pointer to the non IS-part of a generator matrix.
/// The rows and then the columns are sorted inplace.
/// \param L[in]: pointer the currently shortest row.
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
static
int SortCF(normalized_IS_t *__restrict__ G,
           const uint8_t *__restrict__ const L) {
    if (SortRows(G, K, L) == 0) {
	    return 0;
    }
    SortCols(G, K_pad);
    return 1;
} /* end SortCF */

/// Implementation of the `ScaleCF` algorithm from the specification.
/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// NOTE: does not early exit. Only `SortCF` can do that.
/// \param G[in/out]: pointer to the non IS-part of a generator matrix.
/// \param L[in]: pointer the currently shortest row, in histogram form.
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
static
int ScaleCF(normalized_IS_t *__restrict__ G,
            const uint8_t *__restrict__ const L) {
	for (uint32_t row = 0; row < K; row++) {
		if (row_all_same(G->values[row])) { continue; }
		FQ_ELEM s = fq_red(row_acc(G->values[row]));

		if (s != 0) {
			s = fq_inv_non_ct(s);
		} else {
			s = row_acc_inv(G->values[row]);
			if (s == 0) {
				return 0;
			}
		}

		row_mul(G->values[row], s);
	}

	return SortCF(G, L);
} /* end ScaleCF */

/// NOTE: non-constant time
/// This is the `ScaleCFSubPreprocess` function for `ImprovedCFBase`. 
/// The idea is to enumerate the rows of the small z \times N-K matrix G to 
/// find a shorter row than `M`. NOTE: only the fact that a shorter row was 
/// found, is returned. Not the shortest row itself.
/// \param G[in/out]: sub matrix with only z rows.
///     Nothing really is returned in the matrix, but it still gets clobbered.
///     The rows are not in histogram from.
/// \param z[in]: number of rows in G
/// \param M[in]: the currently shortest multiset, in histogram form
/// \return 0: if no multiset was found < `M`
///         1: if one of the  z rows is < `M`
static
int ScaleCFSubPreprocessBase(normalized_IS_t *__restrict__ G,
                             const uint32_t z,
                             const FQ_ELEM *__restrict__ M) {
    FQ_ELEM tmp[Q_pad] __attribute__((aligned(32))) = {0};
    for (uint32_t i = 0; i < z; i++) {
		FQ_ELEM s = fq_red(row_acc(G->values[i]));

		if (s != 0) {
			s = fq_inv_non_ct(s);
		} else {
			s = row_acc_inv(G->values[i]);
			if (s == 0) { return 0; }
		}

		row_mul(G->values[i], s);
        // compute the histogram
        sort(tmp, G->values[i], N-K);

        if (compare_rows(tmp, M) < 0) {
            return 1;
        }
    }

    return 0;
} /* ScaleCFSubPreprocessBase */

/// NOTE: non-constant time
/// This is the `ScaleCFSubPreprocess` function for `ImprovedCF`.
/// The idea is to enumerate the rows of the small z \times N-K matrix G to 
/// find a shorter row than `M`. NOTE: this function returns the shortest row
/// \param G[in/out]: sub matrix with only z rows.
/// \param z[in]: number of rows in G
/// \param M[in]: the currently shortest multiset, will get overwritten with 
///     a shorter row, if found.
/// \return 0: if no multiset was found < `M`
///         1: if one of the  z rows is < `M`
static
int ScaleCFSubPreprocess(normalized_IS_t *__restrict__ G,
                         const uint32_t z,
                         FQ_ELEM *__restrict__ M) {
    /// NOTE: aligment is needed for the optimized avx{2|512} impelementations.
    FQ_ELEM tmp[Q_pad] __attribute__((aligned(32))) = {0};
    int ret = 0;

    for (uint32_t i = 0; i < z; i++) {
        FQ_ELEM s = fq_red(row_acc(G->values[i]));

        if (s != 0) {
            s = fq_inv_non_ct(s);
        } else {
            s = row_acc_inv(G->values[i]);
            if (s == 0) { continue; }
        }

        row_mul(G->values[i], s);
        // compute the histogram
        sort(tmp, G->values[i], N-K);

        if (compare_rows(tmp, M) < 0) {
            ret = 1;
            // copy new smallest row
            for (uint32_t j = 0; j < Q; j++) {
                M[j] = tmp[j];
            }
        }
    }

    return ret;
} /* end ScaleCFSubPreprocess */

/// NOTE: non-constant time
/// implements a total order on matrices
/// we simply compare the rows lexicographically
/// \param V1[in]: first matrix
/// \param V2[in]: second matrix
/// \param z[in]: number of rows within both matrices
/// \return -x if V2 > V1
///			 0 if V2 == V1
///			+x if V2 < V1
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
} /* compare_matrices */

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// Original slowest canonical form function. It's the backup implementation
/// which is called by all other implementation if any failure happens.
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
static
int CFOriginal(normalized_IS_t *G) {
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

		const int ret = ScaleCF(&A, NULL);
		if ((ret == 1) && (compare_matrices(&A, &M, K) < 0)) {
			touched = 1;
			normalized_copy(&M, &A);
		}
	}

    if (!touched) { return 0; }
	normalized_copy(G, &M);
	return 1;
} /* CFOriginal */

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// This is the second-fastest implementation of the canonical form function.
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
static
int ImprovedCFBase(normalized_IS_t *G) {
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
	    return CFOriginal(G);
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

        if (ScaleCFSubPreprocessBase(&B, z, L)) {
            for (uint32_t row2 = 0; row2 < K; row2++) {
                row_mul3(B.values[row2], G->values[row2], row_inv_data);
            }

		    const int ret = ScaleCF(&B, L);
		    if ((ret == 1) && (compare_matrices(&B, &M, K) < 0)) {
                sort(L, B.values[0], N-K);
		    	touched = 1;
		    	normalized_copy(&M, &B);
		    }
        }
	}

	normalized_copy(G, &M);
    return touched;
} /* ImprovedCFBase */

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// This is the fastest implementation of the canonical form function.
/// This function scans the matrix for the rows which probably lead
/// to the shortest canonical form. If this process fails it falls
/// back to either `CFOriginal` or `ImprovedCF`. The very slow `CFOriginal`
/// is only called if every row in the matrix contains zeros.
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
static
int ImprovedCF(normalized_IS_t *G) {
    /// track the rows with the most zeros.
    uint32_t J[K];
    uint32_t z = 0;
    uint32_t smallest_scaling_row = 0;

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
        return CFOriginal(G);
    }

    static normalized_IS_t B __attribute__((aligned(32))) = {0};
    FQ_ELEM row_inv_data[N_K_pad] = {0};

    /// NOTE: this is already "sorted"
    FQ_ELEM L[Q_pad] __attribute__((aligned(32))) = {0};

    /// Check smallest rows of all matricies to find smallest candidate
    for (uint32_t row = 0; row < K; row++) {
        if (Z[row]) { continue; }

        row_inv2(row_inv_data, G->values[row]);
        for (uint32_t row2 = 0; row2 < z; row2++) {
            row_mul3(B.values[row2], G->values[J[row2]], row_inv_data);
        }

        if (ScaleCFSubPreprocess(&B, z, L)) {
            smallest_scaling_row = row;
        }
    }

    /// Calculate CF for best candidate
    row_inv2(row_inv_data, G->values[smallest_scaling_row]);
    for (uint32_t row2 = 0; row2 < K; row2++) {
        row_mul3(B.values[row2], G->values[row2], row_inv_data);
    }

    const int ret = ScaleCF(&B, L);
    /// If candidate was not valid, fall back to regular approach
    if (ret != 1) {
        // Best scaled row was not valid;
        return ImprovedCFBase(G);
    }

    normalized_copy(G, &B);
    return ret;
} /* ImprovedCF */

/// samples to random monomial matrices (A, B) and comptes A*G*B
/// \param G[in/out]: non IS part of a generator matrix
/// \param prng[in/out]: pointer to an already initialized PRNG.
///     in total two random monomials are sampled from the PRNG.
void blind(normalized_IS_t *G,
           SHAKE_STATE_STRUCT *prng) {
    /// NOTE: the type `monomial` allocates a N elements, which is a little
    /// bit a too much, as we only need K or N-K.
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
    normalized_monomial_right(&B, G, &right);

    // apply the left multiplication
    for (uint32_t i = 0; i < K; i++) {
        const FQ_ELEM a = left.coefficients[i];
        const POSITION_T pos = left.permutation[i];

		row_mul2(G->values[i], B.values[pos], a);
    }
} /* end blind */

/// NOTE: non-constant time implementation
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int CF(normalized_IS_t *G) {
#if defined(LESS_CF_PREPROC_PASS_ENABLE)
    return ImprovedCF(G);
#else
    return ImprovedCFBase(G);
#endif
} /* end CF */
