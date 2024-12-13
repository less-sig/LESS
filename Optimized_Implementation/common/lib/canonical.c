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
    if (row_quick_sort(G) == 0) { return 0; }
    col_quicksort_transpose(G);
    // TODO
    // lex_sort_cols(G);
    return 1;
}

/// NOTE: constant time impl
/// computes the result inplace
/// first sort the rows, then the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3_ct(normalized_IS_t *G) {
    if (row_bitonic_sort(G) == 0) { return 0; }
    col_bitonic_sort_transpose(G);
    // TODO something breaks in release mode?
    // lex_sort_cols(G);
    return 1;
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
			s = fq_inv(s);
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
			if (s == 0) { return 0; }
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
    
    for (uint32_t i = 0; i < z; i++) {
		FQ_ELEM s = row_acc(G->values[i]);

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = row_acc_inv(G->values[i]);
			if (s == 0) { return 0; }
		}

		row_mul(G->values[i], s);
        row_sort(G->values[i], N-K);
        if (compare_rows(G->values[i], min_multiset) < 0) {
            return 1;
        }
    }
    return 0;
}

/// implements a total order on matrices
/// we simply compare the columns lexicographically
/// \return -x if V2 > V1
///			 0 if V2 == V1
///			+x
int compare_matrices(const normalized_IS_t *V1,
                     const normalized_IS_t *V2) {
	for (uint32_t row = 0; row < K; row++) {
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
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5(normalized_IS_t *G) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, K*(N-K));

	FQ_ELEM row_inv_data[N_K_pad];
	for (uint32_t row = 0; row < K; row++) {
        if (row_contains_zero(G->values[row])) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < K; row2++) {
			row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
		}

		const int ret = compute_canonical_form_type4(&Aj);
		if ((ret == 1) && (compare_matrices(&Aj, &smallest) < 0)) {
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
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_popcnt(normalized_IS_t *G) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, K*(N-K));

	FQ_ELEM J[N-K];
    uint32_t z = 0;

    // count zeros in each row
    uint32_t max_zeros = 0;
	FQ_ELEM row_has_zero[K] = {0};
	for (uint32_t row = 0; row < K; row++) {
        uint32_t num_zeros = 0;
	    for (uint32_t col = 0; col < N-K; col++) {
            num_zeros += G->values[row][col] == 0; 
        }

        if (num_zeros > 0) {
            row_has_zero[row] = 1;
        }

        if (num_zeros == max_zeros) {
            J[z++] = row; 
        }

        if (num_zeros > max_zeros) {
            z = 1;
            J[0] = row;
            max_zeros = num_zeros;
        }
    }

    /// NOTE: is this always correct?
    if (z == (N-K)) {
	    return compute_canonical_form_type5(G);
    }

    normalized_IS_t scaled_sub_G;
	FQ_ELEM min_multiset[N_K_pad];
	FQ_ELEM row_inv_data[N_K_pad];

	memset(min_multiset, Q-1, N-K);
	for (uint32_t row = 0; row < K; row++) {
        if (row_has_zero[row]) { continue; }

		row_inv2(row_inv_data, G->values[row]);
		for (uint32_t row2 = 0; row2 < z; row2++) {
			row_mul3(scaled_sub_G.values[row2], G->values[J[row2]], row_inv_data);
		}

        if (compute_canonical_form_type4_sub(&scaled_sub_G, z, min_multiset)) {
            // row_inv2(row_inv_data, G->values[row]);
            for (uint32_t row2 = 0; row2 < K; row2++) {
                row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
            }

		    const int ret = compute_canonical_form_type4(&Aj);
		    if ((ret == 1) && (compare_matrices(&Aj, &smallest) < 0)) {
		    	touched = 1;
		    	normalized_copy(&smallest, &Aj);

				memcpy(min_multiset, smallest.values[0], N-K);
				// should be already sorted
                // row_sort(min_multiset, N-K);
		    }
        }
	}

    if (!touched) { return 0; }
	
	normalized_copy(G, &smallest);
	return 1;
}

/// note: constant time
/// note: computes the result inplace
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_ct(normalized_IS_t *G) {
    normalized_IS_t Aj, smallest;
    int touched = 0;

    // init the output matrix to some `invalid` data
    memset(&smallest.values, Q-1, K*(N-K));

    FQ_ELEM row_inv_data[N_K_pad];
    for (uint32_t row = 0; row < K; row++) {
        if (row_contains_zero(G->values[row])) { continue; }

        row_inv2(row_inv_data, G->values[row]);
        for (uint32_t row2 = 0; row2 < K; row2++) {
            row_mul3(Aj.values[row2], G->values[row2], row_inv_data);
        }

        const int ret = compute_canonical_form_type4_ct(&Aj);
        if ((ret == 1) && (compare_matrices(&Aj, &smallest) < 0)) {
            touched = 1;
            normalized_copy(&smallest, &Aj);
        }
    }

    if (!touched) { return 0; }

    normalized_copy(G, &smallest);
    return 1;
}

/// constant time implementation
/// \param G
/// \return
int cf5(normalized_IS_t *G) {
	return compute_canonical_form_type5_ct(G);
}

/// non-constant time implementation
/// \param G
/// \return
int cf5_nonct(normalized_IS_t *G) {
    return compute_canonical_form_type5(G);
    // return compute_canonical_form_type5_popcnt(G);
}
