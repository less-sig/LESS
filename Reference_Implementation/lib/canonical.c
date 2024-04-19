#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"

/// computes the result inplace
/// \return 0 on failure (zero column)
/// 		1 on success
int compute_canonical_form_type2(generator_mat_t *G) {
	// first iterate over all colums: find the first non zero value and scale
	// the it to 0.
	for (uint32_t col = 0; col < N; col++) {
		// find the first non zero entry in the current column
		uint32_t row = 0;
		for (; row < K; row++) {
			if (G->values[row][col] != 0) {
				break;
			}
		}
		
		// we fail if a zero column occur
		if (row == K) { return 0; }
		
		// get the scalling factor
		FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[row][col]);

		// rescale the whole column
		for (uint32_t i = 0; i < K; i++) {
			G->values[i][col] = fq_red(scaling_factor * 
					(FQ_DOUBLEPREC)G->values[i][col]);
		}
	}

	// next sort the columns
	// TODO, currently its a hack: the sorting function does not take a 
	// generator matrix is input, but a type which is nearly identical to the
	// defenintion of a paritiy check matrix.
	col_lex_quicksort((normalized_IS_t *)G, 0, N); 
	return 1;
}

int fqcmp(const void *a, const void *b) {
   return (*(FQ_ELEM *)a) - (*(FQ_ELEM *)b);
}

/// helper function for type3 canonical form. 
/// \input: row1
/// \input: row2
/// \return: 0 if multiset(row1) == multiset(row2)
int compare_rows(const generator_mat_t *G,
		const uint32_t row1, const uint32_t row2) {
	FQ_ELEM sort_row_tmp1[N];
	FQ_ELEM sort_row_tmp2[N];

	memcpy(sort_row_tmp1, G->values[row1], N);
	memcpy(sort_row_tmp2, G->values[row2], N);

	qsort(sort_row_tmp1, N, sizeof(FQ_ELEM), fqcmp);
	qsort(sort_row_tmp2, N, sizeof(FQ_ELEM), fqcmp);


	uint32_t i = 0;
	while(sort_row_tmp1[i] == sort_row_tmp2[i]) {
		i += 1;
	}

	// if they are the same, they generate the same multiset
	if (i == N) {
		return 0;
	}

	return sort_row_tmp1[i] - sort_row_tmp2[i];
}


/// simple bubble sort implementation. Yeah Yeah I know, use quick sort or 
/// radix sort. True. But this is just a demo implementation.
/// \input generator matrix
/// \return the sorting algorithm works inplace
/// 		0 on failure
/// 		1 on success
int row_bubble_sort(generator_mat_t *G) {
	uint32_t swapped;
	do {
		swapped = 0;

		for (uint32_t i = 0; i < N - 1; i++) {
			const int tmp = compare_rows(G, i, i+1);

			// if tmp==0, then row i,i+1 create the same multiset.
			if (tmp == 0) {
				return 0;
			}
			
			if (tmp) {
				row_swap(G, i, i+1);
			}
		}

	} while(swapped);
	return 1;
}

/// numbers which have the following propertie:
/// 	- small as possible
/// 	- m_i < m_i+1 
///		- gdc(m_i, q-1) = 1
///		- char(F_q) dont divide m_i
#define D 10
/// TODO: I choose the first D primes. I have no idea if this is good.
const FQ_ELEM m_array[D] = { 3, 5, 7, 11, 17, 23, 29, 31, 37, 41 };

/// computes (sum_i=0,...,=k  v_i^m_1, ..., sum_i=0,...,=k  v_i^m_d
/// \input G:
/// \input column:
/// \return 1 on success, 0 else
int compute_power_column(generator_mat_t *G, const uint32_t column) {
	uint32_t j = 0;
	FQ_ELEM sum;
	for (; j < D; j++) {
		sum = 0;
		// for each row in the column
		for (uint32_t i = 0; i < K; i++) {
			sum = fq_red(sum + fq_pow(G->values[i][column], m_array[j]));
		}

		if (sum != 0) { break; }
	}

	// we found no j s.t. sum i=0,...,K (v_i^m_j) != 0;
	// return dot
	if (j == D) { return 0; }
	
	sum = fq_inv(sum);
	sum = fq_pow(sum, fq_inv(m_array[j]));

	// scale the vector
	for (uint32_t i = 0; i < K; i++) {
		G->values[i][column] = fq_red(G->values[i][column] * sum);
	}

	return 1;
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(generator_mat_t *G) {
	// first sort the rows 
	if (row_bubble_sort(G) == 0) {
		return 0;
	}

	// next sort the columns 
	col_lex_quicksort((normalized_IS_t *)G, 0, N);
	return 1;
}

/// computes the result inplace
/// \return 0 on failure:
/// 			- compute_power_column failes.
/// 			- identical rows, which create the same multiset
/// 		1 on success
int compute_canonical_form_type4(generator_mat_t *G) {
	for (uint32_t col = 0; col < N; col++) {
		if (compute_power_column(G, col) == 0) {
			// in this case:
			return 0;
		}
	}

	return compute_canonical_form_type3(G);
}

/// implements a total order on matrices
/// we simply compare the columns lexicographically
int comptare_matrices(const generator_mat_t *G1, const generator_mat_t *G2) { 
	generator_mat_t tmp1, tmp2;
	memcpy((void *)&tmp1, G1, K*N);
	memcpy((void *)&tmp2, G1, K*N);
	col_lex_quicksort((normalized_IS_t *)&tmp1, 0, N);
	col_lex_quicksort((normalized_IS_t *)&tmp2, 0, N);

	for (uint32_t col = 0; col < 0; col++) {
		int tmp = lex_compare_column(G1, G2, col, col);
		if (tmp) return tmp;
	}
	
	// if we are here the two matrices are equal
	return 0;
}

/// computes the result inplace
/// \return 0 on failure
/// 		1 on success
int compute_canonical_form_type5(generator_mat_t *G) {
	generator_mat_t tmp, smallest;

	// init the output matrix to some `invalid` data
	memset(&smallest, -1, K*N);
	for (uint32_t col = 0; col < N; col++) {
		memcpy((void *)&tmp, G, K*N);

		// first scale all rows
		for (uint32_t row = 0; row < 0; row++) {
			// TODO not really correct
			if (tmp.values[row][col] == 0) { return 0; }

			scale_row(&tmp, row, fq_inv(tmp.values[row][col]));
		}

		compute_canonical_form_type4(&tmp);
		
		if (comptare_matrices(&tmp, &smallest) == -1) {
			memcpy(&smallest, &tmp, N*K);
		}
	}
	
	memcpy(G, &smallest, N*K);

	return 1;
}
