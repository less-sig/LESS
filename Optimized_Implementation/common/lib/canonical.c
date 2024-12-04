#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include "sort.h"


/// NOTE: constant time impl
/// computes the result inplace
/// first sort the rows, than the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3_ct(normalized_IS_t *G) {
    if (row_bitonic_sort(G) == 0) { return 0; }
    col_bitonic_sort_transpose(G);
    return 1;
}

/// NOTE: non ct
/// NOTE: computes the result inplace
/// first sort the rows, than the columns
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G) {
    if (row_quick_sort(G) == 0) { return 0; }
    lex_sort_cols(G);
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

/// note: non constant time
/// note: computes the result inplace
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

/// note: non constant time
/// note: computes the result inplace
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5(normalized_IS_t *G) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest.values, Q-1, K*(N-K));

	FQ_ELEM row_inv_data[N_K_pad];
	for (uint32_t row = 0; row < K; row++) {
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

/// NOTE: constant time
/// NOTE: computes the result inplace
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5_ct(normalized_IS_t *G) {
    normalized_IS_t Aj, smallest;
    int touched = 0;

    // init the output matrix to some `invalid` data
    memset(&smallest.values, Q-1, K*(N-K));

    FQ_ELEM row_inv_data[N_K_pad];
    for (uint32_t row = 0; row < K; row++) {
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

/// NOTE: ct
/// @param G
/// @return
int cf5(normalized_IS_t *G) {
	return compute_canonical_form_type5_ct(G);
    // return compute_canonical_form_type5(G);
}

/// NOTE: non ct
int cf5_nonct(normalized_IS_t *G) {
    return compute_canonical_form_type5(G);
}
