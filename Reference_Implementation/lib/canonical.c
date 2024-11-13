#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"
#include "sort.h"


/// computes the result inplace
/// NOTE assumes D_c and P_c are identity matrices
/// \return 0 on failure (zero column)
/// 		1 on success
int compute_canonical_form_type2(normalized_IS_t *G,
                                 diagonal_t *D_c,
                                 permutation_t *P_c) {
	// first iterate over all columns: find the first non-zero value and scale
	// it to 0.
	for (uint32_t col = 0; col < N-K; col++) {
		// find the first non-zero entry in the current column
		uint32_t row = 0;
		for (; row < K; row++) {
			if (G->values[row][col] != 0) {
				break;
			}
		}
		
		// we fail if a zero column occur
		if (row >= K) { return 0; }

        // early exit if the value is already 1
        if (G->values[row][col] == 1) { continue; }

		// get the scaling factor
		const FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[row][col]);
        D_c->coefficients[col] = scaling_factor;

		// rescale the whole column
		for (uint32_t i = 0; i < K; i++) {
			G->values[i][col] = fq_red(scaling_factor * 
					(FQ_DOUBLEPREC)G->values[i][col]);
		}
	}

	// next sort the columns
    canonical_col_lex_quicksort(G, 0, N-K-1, P_c);
	return 1;
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G,
								 permutation_t *P_r,
								 permutation_t *P_c) {
    // first sort the rows
    if (row_quick_sort(G, P_r) == 0) { return 0; }
    canonical_col_lex_quicksort(G, 0, N-K-1, P_c);
    return 1;
}

/// computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4(normalized_IS_t *G,
                                 permutation_t *P_r, diagonal_t *D_c,
                                 permutation_t *P_c) {
	for (uint32_t row = 0; row < K; row++) {
        // if we cant find a power
		FQ_ELEM s = 0, sp = 0;
		for (uint32_t col = 0; col < (N-K); col++) {
			s = fq_add(s, G->values[row][col]);
			const FQ_ELEM tmp = fq_inv(G->values[row][col]);
			sp = fq_add(sp, tmp);
		}

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = sp;
			if (s == 0) {
				return -1;
			}
		}

		if (D_c != NULL) {D_c->coefficients[row] = s;}
		for (uint32_t col = 0; col < (N-K); col++) {
			G->values[row][col] = fq_mul(s, G->values[row][col]);
		}

	}

	return compute_canonical_form_type3(G, P_r, P_c);
}

/// NOTE: non constant time
/// computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column fails.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4_non_ct(normalized_IS_t *G,
									    permutation_t *P_r,
									    diagonal_t *D_c,
										permutation_t *P_c) {
	for (uint32_t row = 0; row < K; row++) {
		FQ_ELEM s = 0;
		for (uint32_t col = 0; col < (N-K); col++) {
			s = fq_add(s, G->values[row][col]);
		}

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = 0;
			for (uint32_t col = 0; col < (N-K); col++) {
				s = fq_add(s, fq_inv(G->values[row][col]));
			}
			if (s == 0) { return -1; }
		}

		if (D_c != NULL) {D_c->coefficients[row] = s;}
		for (uint32_t col = 0; col < (N-K); col++) {
			G->values[row][col] = fq_mul(s, G->values[row][col]);
		}

	}

	return compute_canonical_form_type3(G, P_r, P_c);
}

/// implements a total order on matrices
/// we simply compare the columns lexicographically
int compare_matrices(const normalized_IS_t *V1,
                     const normalized_IS_t *V2) {

	for (uint32_t col = 0; col < N-K; col++) {

        uint32_t i=0;
        while((i < K) &&
             ((V1->values[i][col] - V2->values[i][col]) == 0)) {
            i++;
        }

        if (i >= K) { continue; }

        return (int)(V1->values[i][col]) - (int)(V2->values[i][col]);
	}
	
	// if we are here the two matrices are equal
	return 0;
}

/// computes the result inplace
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5(normalized_IS_t *G,
                                 diagonal_t *D_r,
                                 permutation_t *P_r,
                                 diagonal_t *D_c,
                                 permutation_t *P_c) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest, -1, K*(N-K));

	for (uint32_t row = 0; row < K; row++) {
		for (uint32_t col = 0; col < N-K; col++) {
			const FQ_ELEM a_inv = fq_inv(G->values[row][col]);
			for (uint32_t row2 = 0; row2 < N-K; row2++) {
				Aj.values[row2][col] = fq_mul(G->values[row2][col], a_inv);
			}
		}

		compute_canonical_form_type4(&Aj, P_r, D_c, P_c);
		if (compare_matrices(&Aj, &smallest) < 0) {
			touched = 1;
			memcpy((void *)smallest.values, (void*)Aj.values, sizeof(normalized_IS_t));
		}
	}

    if (!touched) { return 0; }
	
	memcpy((void *)G->values, (void*)smallest.values, sizeof(normalized_IS_t));
	return 1;
}

int cf5(normalized_IS_t *G) {
	return compute_canonical_form_type5(G, NULL, NULL, NULL, NULL);
}
