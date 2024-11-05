#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "monomial_mat.h"
#include "canonical.h"
#include "fq_arith.h"
#include "parameters.h"

#define SWAP(a, b, tmp) tmp = a; a = b; b = tmp;

/// tries to find the information set. If found a permutation which
int compute_information_set_monomial(permutation_t *P_is, monomial_t *G) {
    POSITION_T tmp;
    FQ_ELEM q_tmp;
    for (uint32_t k = 0; (k < N) ; ++k) {
        POSITION_T to = G->permutation[k];
        while ((to < K) && (to != k)) {
            const POSITION_T too = G->permutation[to];

            SWAP(G->permutation[k], G->permutation[to], tmp)
            SWAP(G->coefficients[k], G->coefficients[to], q_tmp)
            SWAP(P_is->permutation[k], P_is->permutation[to], tmp)

            if (k <= too) {
                break;
            }

            to = too;
        }
    }

    return 1;
}


/// computes the result inplace
/// \return 0 on failure (zero column)
/// 		1 on success
int monomial_compute_canonical_form_type2(const monomial_t *M,
		diagonal_t *Dc, permutation_t *Pc) {
	// first iterate over all columns: find the first non-zero value and scale
	// it to 0.
	for (uint32_t col = 0; col < N; col++) {
		// find the first non-zero entry in the current column
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

/// helper function for type3 canonical form. 
/// \input: row1
/// \input: row2
/// \return: 0 if multiset(row1) == multiset(row2)
int permutation_compare_rows(const monomial_t *M,
		const uint32_t row1, const uint32_t row2) {
	uint32_t t1 = -1, t2 = -1;
	for (uint32_t i = 0; i < N; i++) {
		if (M->permutation[i] == row1) { t1 = i; }
		if (M->permutation[i] == row2) { t2 = i; }
	}

	return M->coefficients[t1] - M->coefficients[t2];
}

int permutation_compare_columns(const monomial_t *M,
		const uint32_t row1, const uint32_t row2) {
	if (M->permutation[row1] < M->permutation[row2]) { return -1; }
	if (M->permutation[row1] == M->permutation[row2]) {
		return M->coefficients[row1] - M->coefficients[row2];
	}

	return 1;
}

/// simple bubble sort implementation. Yeah Yeah I know, use quick sort or 
/// radix sort. True. But this is just a demo implementation.
/// \input generator matrix
/// \return the sorting algorithm works inplace
/// 		0 on failure
/// 		1 on success
int permutation_row_bubble_sort(permutation_t *G, permutation_t *Pr) {
	uint32_t swapped;
	FQ_ELEM ftmp;
	POSITION_T ptmp;
	do {
		swapped = 0;

		for (uint32_t i = 0; i < N - 1; i++) {
			const int tmp = permutation_compare_rows(G, i, i+1);

			// if tmp==0, then row i,i+1 create the same multiset.
			if (tmp == 0) {
				return 0;
			}
			
			// if row_i < row_i+1
			if (tmp < 0) {
				// SWAP(G->coefficients[i], G->coefficients[i+1], ftmp);
				SWAP(G->permutation[i], G->permutation[i+1], ptmp);
				SWAP(Pr->permutation[i], Pr->permutation[i+1], ptmp);
                swapped = 1;
			}
		}
	} while(swapped);
	return 1;
}

int column_bubble_sort(monomial_t *G, permutation_t *Pc) {
	uint32_t swapped;
	FQ_ELEM ftmp;
	POSITION_T ptmp;
	do {
		swapped = 0;

		for (uint32_t i = 0; i < N - 1; i++) {
			const int tmp = permutation_compare_columns(G, i, i+1);
			if (tmp < 0) {
				SWAP(G->coefficients[i], G->coefficients[i+1], ftmp);
				SWAP(G->permutation[i], G->permutation[i+1], ptmp);
				SWAP(Pc->permutation[i], Pc->permutation[i+1], ptmp);
                swapped = 1;
			}
		}
	} while(swapped);
	return 1;
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int monomial_compute_canonical_form_type3(permutation_t *M,
		permutation_t *Pr, permutation_t *Pc) {
	// first sort the rows 
	if (row_bubble_sort(M, Pr) == 0) {
		return 0;
	}

	// next sort the columns 
	column_bubble_sort(M, Pc);
	return 1;
}

#define WRITE_BIT(a, pos) (a[pos/64u] |= (1u << (pos %64u)))
#define TYPE3_COMPRESSION_LIMBS ((K + 63u) / 64u)

/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int monomial_compress_type3(uint64_t out[TYPE3_COMPRESSION_LIMBS],
		permutation_t *G) {
	permutation_t Pr, Pc;
	monomial_compute_canonical_form_type3(G, &Pr, &Pc);

	for (uint32_t i = 0; i < TYPE3_COMPRESSION_LIMBS; i++) {
		out[i] = 0;
	}

	for (uint32_t i = 0; i < K; i++) {
		WRITE_BIT(out, G->permutation[i]);
	}
	return 1;
}
