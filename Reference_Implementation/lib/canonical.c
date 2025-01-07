#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include "sort.h"


/// NOTE: non ct
/// NOTE: computes the result inplace
/// first sort the rows, then the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G) {
    if (row_quick_sort(G, K) == 0) {
	    return 0;
    }
#ifdef LESS_USE_HISTOGRAM
    col_quicksort_transpose(G, K_pad);
#else
    col_lex_quicksort(G, 0, N-K-1);
#endif
    return 1;
}

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4(normalized_IS_t *G) {
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

	return compute_canonical_form_type3(G);
}

/// NOTE: non-constant time
/// \param G[in] sub matrix with only z rows
/// \param z[in] number of rows in G
/// \param min_multiset[in]
/// \return 0: if no multiset was found < `min_multiset`
///         1: if one of the  z rows is < `min_multiset`
int compute_canonical_form_type4_sub(normalized_IS_t *G,
                                     const uint32_t z,
                                     const FQ_ELEM *min_multiset) {
#ifdef LESS_USE_HISTOGRAM
    FQ_ELEM tmp[Q_pad];
#else
    FQ_ELEM tmp[N_K_pad] = {0};
#endif

    for (uint32_t i = 0; i < z; i++) {
		FQ_ELEM s = row_acc(G->values[i]);

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = row_acc_inv(G->values[i]);
			if (s == 0) { return 0; }
		}

		row_mul(G->values[i], s);
        row_sort(tmp, G->values[i], N-K);
        if (compare_rows(tmp, min_multiset) < 0) {
            return 1;
        }
    }

    return 0;
}

/// NOTE: non-constant time
/// implements a total order on matrices
/// we simply compare the columns lexicographically
/// \param V1[in]
/// \param V2[in]
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
	normalized_IS_t Aj = {0}, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, sizeof(normalized_IS_t));

	FQ_ELEM row_inv_data[N_K_pad] = {0};
	for (uint32_t row = 0; row < K; row++) {
        if (row_contains_zero(G->values[row])) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < K; row2++) {
			row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
		}

		const int ret = compute_canonical_form_type4(&Aj);
		if ((ret == 1) && (compare_matrices(&Aj, &smallest, K) < 0)) {
			touched = 1;
			normalized_copy(&smallest, &Aj);
		}
	}

    if (!touched) { return 0; }
	normalized_copy(G, &smallest);
	return 1;
}

/// NOTE: non-constant time
/// NOTE: computes the result inplace
/// \param G[in/out] non IS part of a generator matrix
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_popcnt(normalized_IS_t *G) {
	normalized_IS_t Aj = {0}, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, sizeof(normalized_IS_t));

	uint32_t J[N-K];
    uint32_t z = 0;

    // count zeros in each row
    uint32_t max_zeros = 0;
	uint8_t row_has_zero[K] = {0};
	for (uint32_t row = 0; row < K; row++) {
        const uint32_t num_zeros = row_count_zero(G->values[row]);

        if (num_zeros > 0) {
            row_has_zero[row] = 1;
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

    /// NOTE: is this always correct?
    if (z == (N-K)) {
	    return compute_canonical_form_type5(G);
    }

    normalized_IS_t scaled_sub_G __attribute__((aligned(32))) = {0};
	FQ_ELEM row_inv_data[N_K_pad] = {0};

	/// NOTE: this is already "sorted"
#ifdef LESS_USE_HISTOGRAM
	FQ_ELEM min_multiset[Q_pad];
	memset(min_multiset, 0, Q_pad);
#else
	FQ_ELEM min_multiset[N_K_pad];
	memset(min_multiset, Q-1, N_K_pad);
#endif
	for (uint32_t row = 0; row < K; row++) {
        if (row_has_zero[row]) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < z; row2++) {
			row_mul3(scaled_sub_G.values[row2], G->values[J[row2]], row_inv_data);
		}

        if (compute_canonical_form_type4_sub(&scaled_sub_G, z, min_multiset)) {
            for (uint32_t row2 = 0; row2 < K; row2++) {
                row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
            }

		    const int ret = compute_canonical_form_type4(&Aj);
		    if ((ret == 1) && (compare_matrices(&Aj, &smallest, K) < 0)) {
#ifdef LESS_USE_HISTOGRAM
            row_sort(min_multiset, Aj.values[0], N-K);
#else
            memcpy(min_multiset, Aj.values[0], N_K_pad);
#endif
		    	touched = 1;
		    	normalized_copy(&smallest, &Aj);
		    }
        }
	}

    if (!touched) { return 0; }
	
	normalized_copy(G, &smallest);
	return 1;
}

/// samples to random monomial matrices (A, B) and comptes A*G*B
/// \param G[in/out] non IS part of a generator matrix
/// \param prng[in/out]:
void blind(normalized_IS_t *G,
           SHAKE_STATE_STRUCT *prng) {
    /// note;
    monomial_t left, right;
    normalized_IS_t B;

    // We compute the following matrix multiplication G = left * G * right
    // where `left` and `right` are randomly sampled monomials
    fq_star_rnd_state_elements(prng, left.coefficients, N-K);
    fq_star_rnd_state_elements(prng, right.coefficients, N-K);
    for(uint32_t i = 0; i < N-K; i++) {
        left.permutation[i] = i;
        right.permutation[i] = i;
    }

    yt_shuffle_state_limit(prng, left.permutation, N-K);
    yt_shuffle_state_limit(prng, right.permutation, N-K);

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
int cf5_nonct(normalized_IS_t *G) {
    return compute_canonical_form_type5_popcnt(G);
}
