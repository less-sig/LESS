#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "monomial_mat.h"
#include "fq_arith.h"
#include "parameters.h"


/// computes the result inplace
/// \return 0 on failure (zero column)
/// 		1 on success
int compute_canonical_form_type2(const monomial_t *M,
		diagonal_t *Dc, permutation_t *Pc) {
	// first iterate over all colums: find the first non zero value and scale
	// the it to 0.
	for (uint32_t col = 0; col < N; col++) {
		// find the first non zero entry in the current column
		Dc->coefficients[col] = fq_inv(M->coefficients[col]);
	}

	// next sort the columns
	for (uint32_t col = 0; col < N; col++) {
		uint32_t pos = M->permutation[col];
		Pc->permutation[pos] = col;
	}

	// TODO: the sorting above is in decreasing lexicographic order and not increasing
	// TODO: apply the sort to M;
	return 1;
}

int fqcmp(const void *a, const void *b) {
   return (*(FQ_ELEM *)a) - (*(FQ_ELEM *)b);
}

/// helper function for type3 canonical form. 
/// \input: row1
/// \input: row2
/// \return: 0 if multiset(row1) == multiset(row2)
int compare_rows(const monomial_t *M,
		const uint32_t row1, const uint32_t row2) {
}


/// simple bubble sort implementation. Yeah Yeah I know, use quick sort or 
/// radix sort. True. But this is just a demo implementation.
/// \input generator matrix
/// \return the sorting algorithm works inplace
/// 		0 on failure
/// 		1 on success
int row_bubble_sort(monomial_t *G) {
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
int compute_power_column(monomial_t *M, const uint32_t column) {

	return 1;
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(monomial_t *M,
		permutation_t *Pc) {
	// first sort the rows 
	if (row_bubble_sort(M) == 0) {
		return 0;
	}

	// next sort the columns 
	col_lex_quicksort((normalized_IS_t *)G, 0, N);
	return 1;
}

