#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"

// TODO:
//      - avx/neon djbsort, but its just a bitonic sorter, should be easy
//      - bitonic sorting network, with generic comparison
//      - benchmark which sorting algorithm is the fastest:
//          quick, bitonic, bubble lol


/// swaps a and b if
void cswap(uintptr_t *a, uintptr_t *b, const uint64_t mask) {
    *a ^= (mask & *b);
    *b ^= (mask & *a);
    *a ^= (mask & *b);
}

/// swaps a and b if f==1, if f==0, nothing will happen
void cswap_bit(uintptr_t *a, uintptr_t *b, const uint64_t f) {
	const uint64_t mask = -f;
	cswap(a, b, mask);
}

// taken from djbsort
#define int8_MINMAX(a,b)\
do {                    \
    int8_t ab = b ^ a;  \
    int8_t c = b - a;   \
    c ^= ab & (c ^ b);  \
    c >>= 7;            \
    c &= ab;            \
    a ^= c;             \
    b ^= c;             \
} while(0)

// taken from djbsort
void int8_sort(FQ_ELEM *x,
               const long long n) {
    long long top, p, q, r, i;
    if (n < 2) {
        return;
    }
    top = 1;
    while (top < n - top) {
        top += top;
    }
    
    for (p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                int8_MINMAX(x[i], x[i + p]);
            }
        }
        
        i = 0;
        for (q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    int8_t a = x[i + p];
                    for (r = q; r > p; r >>= 1) {
                        int8_MINMAX(a, x[i + r]);
                    }
                    
                    x[i + p] = a;
                }
            }
        }
    }
}

/// the same as `Hoare_partition` except that we track the permutation made
/// \param V
/// \param col_l
/// \param col_h
/// \return
int canonical_Hoare_partition(normalized_IS_t *V,
                              const POSITION_T col_l,
                              const POSITION_T col_h,
                              permutation_t *P) {
    FQ_ELEM pivot_col[K] = {0};
    for(uint32_t i = 0; i < K; i++){
        pivot_col[i] = V->values[i][col_l];
    }
    POSITION_T i = col_l-1, j = col_h+1;
    while(1) {
        do {
            i++;
        } while(lex_compare_with_pivot(V,i,pivot_col) == 1);

        do {
            j--;
        } while(lex_compare_with_pivot(V,j,pivot_col) == -1);

        if(i >= j){
            return j;
        }

        column_swap(V,i,j);
        permutation_swap(P, i, j);
    }
}

/// In-place quicksort
/// the same as `col_lex_quicksort` except we are tracking the permutations
/// \param V
/// \param start
/// \param end
void canonical_col_lex_quicksort(normalized_IS_t *V,
                                 const int start,
                                 const int end,
                                 permutation_t *P) {
    if(start < end){
        int p = canonical_Hoare_partition(V, start, end, P);
        col_lex_quicksort(V,start,p);
        col_lex_quicksort(V,p+1,end);
    }
}

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
		FQ_DOUBLEPREC scaling_factor = fq_inv(G->values[row][col]);
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

/// helper function:
/// \param a:
/// \param b:
/// \return a - b:
///         -1: b > a;
///          0: b == a;
///          1: b < a;
int fqcmp(const void *a,
          const void *b) {
   return (*(FQ_ELEM *)a) - (*(FQ_ELEM *)b);
}

/// helper function for type3 canonical form. 
/// \input: row1
/// \input: row2
/// \return: 0 if multiset(row1) == multiset(row2)
///          1 if 
///         -1 if 
int compare_rows(const FQ_ELEM rows[K][N-K],
		         const uint32_t row1, 
                 const uint32_t row2) {
    ASSERT(row1 < K);
    ASSERT(row2 < K);

	uint32_t i = 0;
	while((rows[row1][i] == rows[row2][i]) && (i < (N-K))) {
		i += 1;
	}

	// if they are the same, they generate the same multiset
	if (i >= (N-K)) {
		return 0;
	}

	return (int)rows[row1][i] - (int)rows[row2][i];
}

int compare_rows2(FQ_ELEM **rows,
                 const uint32_t row1,
                 const uint32_t row2) {
    ASSERT(row1 < K);
    ASSERT(row2 < K);

    uint32_t i = 0;
    while((rows[row1][i] == rows[row2][i]) && (i < (N-K))) {
        i += 1;
    }

    // if they are the same, they generate the same multiset
    if (i >= (N-K)) {
        return 0;
    }

    return (int)rows[row1][i] - (int)rows[row2][i];
}

/// simple bubble sort implementation. Yeah Yeah I know, use quick sort or 
/// radix sort. True. But this is just a demo implementation.
/// \input generator matrix
/// \return the sorting algorithm works only inplace for the sorting of the columns
/// 		0 on failure: row_i and row_j generate the same multiset
/// 		1 on success
int row_bubble_sort(normalized_IS_t *G, permutation_t *P_r) {
    // first sort each row into a tmp buffer
    FQ_ELEM tmp[K][N-K];
    for (uint32_t i = 0; i < K; ++i) {
        memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
        qsort(tmp[i], N-K, sizeof(FQ_ELEM), fqcmp);
    }

    normalized_pretty_print_v(tmp);
	uint32_t swapped;

	do {
		swapped = 0;

        // for all rows
		for (uint32_t i = 0; i < K-1; i++) {
			const int cmp = compare_rows(tmp, i, i+1);

			// if cmp==0, then row i,i+1 create the same multiset.
			if (cmp == 0) {
				return 0;
			}
			
			if (cmp > 0) {
				row_swap(G, i, i+1);
                permutation_swap(P_r, i, i+1);

                // TODO speedup: probably just move pointers
                for (uint32_t j = 0; j < N-K; ++j) {
                    FQ_ELEM tmp2 = tmp[i][j];
                    tmp[i][j] = tmp[i+1][j];
                    tmp[i+1][j] = tmp2;
                }
                swapped = 1;
			}
		}
	} while(swapped);


    normalized_pretty_print_v(tmp);
	return 1;
}

//
int row_bitonic_sort(normalized_IS_t *G,
                     permutation_t *P_r) {
    // first sort each row into a tmp buffer
    FQ_ELEM  tmp[K][N-K];
    FQ_ELEM* ptr[K];
    for (uint32_t i = 0; i < K; ++i) {
        memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
        int8_sort(tmp[i], N-K);
        ptr[i] = tmp[i];
    }

    normalized_pretty_print_v(tmp);
    uint64_t top, r, i;
    uint64_t n = K;
    top = 1;

    // TODO precompute
    while (top < n - top) {
        top += top;
    }
    
    for (uint64_t p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast
			    const int32_t cmp1 = compare_rows2(ptr, i, i + p);
                const int32_t cmp2 = compare_rows(tmp, i, i + p);
                const uint32_t cmp = cmp1;
                if (cmp == 0) { return 0; }

                //const uint64_t flag = -((1ull + cmp) >> 1);
                const uint64_t mask = -(1ull - (cmp >> 31));
                cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[i+p]), mask);
            }
        }
        
        i = 0;
        for (uint64_t q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    uintptr_t *a = (uintptr_t *)(&ptr[i + p]);
                    for (r = q; r > p; r >>= 1) {
			            const uint32_t cmp = compare_rows2(ptr, i, i + p);
                        if (cmp == 0) { return 0; }

                        const uint64_t flag = -(1ull - (cmp >> 31));
                        cswap(a, (uintptr_t *)(&ptr[i+r]), flag);
                    }
                    
                    ptr[i + p] = (FQ_ELEM *)a;
                }
            }
        }
    }

    normalized_pretty_print_v(ptr);
    return 1;
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G, permutation_t *P_r, permutation_t *P_c) {
    // first sort the rows
    if (row_bubble_sort(G, P_r) == 0) {
        return 0;
    }

    // next sort the columns
    canonical_col_lex_quicksort(G, 0, N-K-1, P_c);
    return 1;
}

/// numbers which have the following propertie:
/// 	- small as possible
/// 	- m_i < m_i+1 
///		- gdc(m_i, q-1) = 1
///		- char(F_q) dont divide m_i
/// NOTE: I choose the first D primes. I have no idea if this is good.
#define D 5
const FQ_ELEM m_array[D] = { 3, 5, 7, 11, 17 };//, 23, 29, 31, 37, 41 };

/// computes (sum_i=0,...,=k  v_i^m_1, ..., sum_i=0,...,=k  v_i^m_d
/// \input G:
/// \input column:
/// \return 1 on success, 0 else
int compute_power_column(normalized_IS_t *G, const uint32_t column,
                         diagonal_t *D_c) {
	uint32_t j = 0;
	FQ_ELEM sum;
	for (; j < D; j++) {
		sum = 0;
		// for each row in the column
		for (uint32_t i = 0; i < K; i++) {
			sum = fq_add(sum, fq_pow(G->values[i][column], m_array[j]));
		}

		if (sum != 0) { break; }
	}

	// we found no j s.t. sum i=0,...,K (v_i^m_j) != 0;
	// return dot
	if (j >= D) { return 0; }
	
	sum = fq_inv(sum);
	sum = fq_pow(sum, fq_inv(m_array[j]));

    D_c->coefficients[column] = sum;

	// scale the vector
	for (uint32_t i = 0; i < K; i++) {
		G->values[i][column] = fq_red(G->values[i][column] * sum);
	}

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
	for (uint32_t col = 0; col < N-K; col++) {
        // if we cant find a power
		if (compute_power_column(G, col, D_c) == 0) {
			return 0;
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
             ((V1->values[i][col] - V2->values[i][col]) == 0)){
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
                                 diagonal_t *D_r, permutation_t *P_r,
                                 diagonal_t *D_c, permutation_t *P_c) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest, -1, K*(N-K));


	for (uint32_t col = 0; col < N-K; col++) {
        if (normalized_is_zero_in_column(G, col)) { continue; }
		memcpy((void *)&Aj, G, K*(N-K));
        touched = 1;

		// first scale all rows
		for (uint32_t row = 0; row < K; row++) {
            FQ_ELEM tmp = fq_inv(Aj.values[row][col]);
			normalized_mat_scale_row(&Aj, row, tmp);
            D_r->coefficients[row] = fq_mul(D_r->coefficients[row], tmp);
		}

        // TODO, need the smallest
		compute_canonical_form_type4(&Aj, P_r, D_c, P_c);
		
		if (compare_matrices(&Aj, &smallest) < 0) {
			memcpy(&smallest, &Aj, K*(N-K));
		}
	}

    if (!touched) { return 0; }
	
	memcpy(G, &smallest, K*(N-K));
	return 1;
}
