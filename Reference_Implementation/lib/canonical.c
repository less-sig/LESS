#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "codes.h"
#include "fq_arith.h"
#include "parameters.h"

// TODO:
//      - avx/neon djbsort, but its just a bitonic sorter, should be easy
//      - bitonic sorting network, with generic comparison
//      - benchmark which sorting algorithm is the fastest:
//          quick, bitonic, bubble lol




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

void counting_sort_u8(FQ_ELEM *arr,
					  const size_t size) {
	uint32_t cnt[128] __attribute__((aligned(64))) = { 0 };
	size_t i;

	for (i = 0 ; i < size ; ++i) {
		cnt[arr[i]]++;
	}

	i = 0;
	for (size_t a = 0 ; a < 128; ++a) {
		while (cnt[a]--) {
			arr[i++] = a;
		}
	}
}

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

                // NOTE speedup: probably just move pointers
                for (uint32_t j = 0; j < N-K; ++j) {
                    FQ_ELEM tmp2 = tmp[i][j];
                    tmp[i][j] = tmp[i+1][j];
                    tmp[i+1][j] = tmp2;
                }
                swapped = 1;
			}
		}
	} while(swapped);
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
                // NOTE: here is a sign cast, this is needed, so the
            	// sign extension needed for the mask is a unsigned one.
			    const int32_t cmp1 = compare_rows2(ptr, i, i + p);
                const uint32_t cmp = cmp1;
                if (cmp == 0) { return 0; }

                const uintptr_t mask = -(1ull - (cmp >> 31));
                cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[i+p]), mask);
            	row_cswap(G, i, i+p, mask);
            	permutation_cswap(P_r, i, i+p, mask);
            }
        }
        
        i = 0;
        for (uint64_t q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1) {
			            const uint32_t cmp = compare_rows2(ptr, i+p, i + r);
                        if (cmp == 0) { return 0; }

                        const uintptr_t mask = -(1ull - (cmp >> 31));
                        cswap((uintptr_t *)(&ptr[i+p]), (uintptr_t *)(&ptr[i+r]), mask);
                    	permutation_cswap(P_r, i+p, i+r, mask);
            			row_cswap(G, i+p, i+r, mask);
                    }
                }
            }
        }
    }

    return 1;
}


/// lexicographic comparison of a row with the pivot
/// returns    1 if the pivot is greater,
/// 	      -1 if it is smaller,
/// 		   0 if it matches
int row_lex_compare_with_pivot(FQ_ELEM V[K][N-K],
                           const POSITION_T row_idx,
                           FQ_ELEM pivot[N-K]){
   uint32_t i=0;
   while((i<(N-K)) && (V[row_idx][i]-pivot[i] == 0)){
       i++;
   }
   if (i==(N-K)) {
	   return 0;
   }

   if ((int)V[row_idx][i]-(int)pivot[i] > 0){
      return -1;
   }

   return 1;
}


void print_kek(FQ_ELEM V[K][N-K]) {
	for (uint32_t i = 0; i < K; i++) {
		for (uint32_t j = 0; j < (N-K); j++) {
			printf("%03u ", V[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

int row_hoare_partition(FQ_ELEM V[K][N-K],
						normalized_IS_t *G,
						const POSITION_T row_l,
						const POSITION_T row_h,
						permutation_t *P_r) {
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = V[row_l][i];
    }

	// print_kek(V);

    POSITION_T i = row_l-1, j = row_h+1;
	int ret;
    while(1){
        do {
            i++;
        	ret = row_lex_compare_with_pivot(V, i, pivot_row);
        } while(ret == 1);

        do {
            j--;
        	ret = row_lex_compare_with_pivot(V, j, pivot_row);
        } while(ret == -1);

    	// if (ret == 0) {
    	// 	return 0;
    	// }

        if(i >= j){
            return j;
        }

        // row_swap(V, i, j);
    	row_swap(G, i, j);
    	SWAP(P_r->permutation[i], P_r->permutation[j]);
		for(uint32_t k = 0; k < N-K; k++){
			SWAP(V[i][k], V[j][k]);
		}
    }
}

int row_lex_quicksort(FQ_ELEM V[K][N-K],
					  normalized_IS_t *G,
                      const uint32_t start,
                      const uint32_t end,
                      permutation_t *P_r) {
    if(start < end){
        const int p = row_hoare_partition(V, G, start, end, P_r);
        row_lex_quicksort(V, G, start, p, P_r);
        row_lex_quicksort(V, G, p+1, end, P_r);
    }

	return 1;
}

// non inplace
int row_quick_sort(normalized_IS_t *G,
				   permutation_t *P_r) {
	// first sort each row into a tmp buffer
	FQ_ELEM  tmp[K][N-K];
	for (uint32_t i = 0; i < K; ++i) {
		memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
		// int8_sort(tmp[i], N-K);
		counting_sort_u8(tmp[i], N-K);
	}

	return row_lex_quicksort(tmp, G, 0, N-K-1, P_r);
}


// NOTE: unstable sort, as we sort on 0
/// sort is applied inplace
/// @param G
/// @param P_c
/// @return
void col_bitonic_sort(normalized_IS_t *G,
                      permutation_t *P_c) {
    uint64_t top, r, i;
    const uint64_t n = N-K;
    top = 1;

    // TODO precompute
    while (top < n - top) {
        top += top;
    }

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast, this is needed, so the
            	// sign extension needed for the mask is a unsigned one.
			    const uint32_t cmp = lex_compare_col(G, i, i + p);
                const uintptr_t mask = -(1ull - (cmp >> 31));
            	column_cswap(G, i, i+p, mask);
            	permutation_cswap(P_c, i, i+p, mask);
            }
        }

        i = 0;
        for (uint64_t q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1) {
						const uint32_t cmp = lex_compare_col(G, i+p, i+r);
            			const uintptr_t mask = -(1ull - (cmp >> 31));
                    	permutation_cswap(P_c, i+p, i+r, mask);
            			column_cswap(G, i+p, i+r, mask);
                    }
                }
            }
        }
    }
}

/// computes the result inplace
/// first sort the rows, than the colums
/// \return 0 on failure (identical rows, which create the same multiset)
/// 		1 on success
int compute_canonical_form_type3(normalized_IS_t *G,
								 permutation_t *P_r,
								 permutation_t *P_c) {
    // first sort the rows
    //if (row_bitonic_sort(G, P_r) == 0) {
    if (row_quick_sort(G, P_r) == 0) {
        return 0;
    }

    // next sort the columns
    canonical_col_lex_quicksort(G, 0, N-K-1, P_c);
	//col_bitonic_sort(G, P_c);
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
		FQ_ELEM s = 0, sp = 0, tmp, q2=Q-2;
		for (uint32_t col = 0; col < (N-K); col++) {
			s = fq_add(s, G->values[row][col]);
			tmp = fq_inv(G->values[row][col]); //fq_pow(G->values[row][col], q2);
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
            // D_r->coefficients[row] = fq_mul(D_r->coefficients[row], tmp);
		}

		compute_canonical_form_type4(&Aj, P_r, D_c, P_c);
		
		if (compare_matrices(&Aj, &smallest) < 0) {
			memcpy(&smallest, &Aj, K*(N-K));
		}
	}

    if (!touched) { return 0; }
	
	memcpy(G, &smallest, K*(N-K));
	return 1;
}
