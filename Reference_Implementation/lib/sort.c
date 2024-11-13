#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "utils.h"
#include "fq_arith.h"
#include "codes.h"

/// taken from djbsort
/// \param a first input
/// \param b second input
/// \returns nothing but: a = min(a, b),
///						  b = max(a, b)
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

/// NOTE: specialised counting sort for Fq. Thus,
/// this implementation assumes that every input element
/// is reduced mod 127
/// \param arr input array
/// \param size length
void counting_sort_u8(FQ_ELEM *arr,
                      const size_t size) {
	/// NOTE: the type `uint32_t` is not completly arbitrary choose.
	/// Floyd did a quick benchmark between `uint16_t`, `uint32_t`, `uint64_t`
	/// and `uint32_t` seemed to be the fastest. But thats only true
	/// on a Ryzen7600X. On your machine thats maybe different.
	/// NOTE: `uint8_t` is not possible as there could be 256 times
	/// the same field element.
	uint32_t cnt[128] __attribute__((aligned(128))) = { 0 };
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

/// NOTE: taken from djbsort
/// \param x input array
/// \param n length
void bitonic_sort_i8(FQ_ELEM *x,
                     const long long n) {
    long long p, q, r, i;
    if (n < 2) {
        return;
    }
    const long long top = 1ul << (32 - __builtin_clz(K/2));
    // while (top < n - top) { top += top; }

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

/// helper function for the libc `qsort`
/// \return a - b:
///         -1: b > a;
///          0: b == a;
///          1: b < a;
int fqcmp(const void *a,
          const void *b) {
   return (*(FQ_ELEM *)a) - (*(FQ_ELEM *)b);
}

/// NOTE: helper function for type3 canonical form.
/// \input: rows: K x (N-K) matrix
/// \input: row1: first row to compare
/// \input: row2: secnd row to compare
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
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

/// NOTE: helper function for type3 canonical form.
/// NOTE: specially made for the bitonic sort, which
///		operates on pointers.
/// \input: rows: K x (N-K) matrix
/// \input: row1: first row to compare
/// \input: row2: secnd row to compare
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
int compare_rows_bitonic_sort(FQ_ELEM **rows,
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
/// \input normalised non IS part of a generator matrix
/// \return the sorting algorithm works only inplace for the sorting of the columns
/// 		0 on failure: row_i and row_j generate the same multiset
/// 		1 on success
int row_bubble_sort(normalized_IS_t *G,
					permutation_t *P_r) {
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
                	SWAP(tmp[i][j], tmp[i+1][j]);
                }
                swapped = 1;
			}
		}
	} while(swapped);
	return 1;
}

/// simple bubble sort implementation. Yeah Yeah I know, use quick sort or
/// radix sort. True. But this is just a demo implementation.
/// \input normalised non IS part of a generator matrix
/// \return the sorting algorithm works only inplace for the sorting of the columns
/// 		0 on failure: row_i and row_j generate the same multiset
/// 		1 on success
int row_bitonic_sort(normalized_IS_t *G,
                     permutation_t *P_r) {
    // first sort each row into a tmp buffer
    FQ_ELEM  tmp[K][N-K];
    FQ_ELEM* ptr[K];
    for (uint32_t i = 0; i < K; ++i) {
        memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
        bitonic_sort_i8(tmp[i], N-K);
        ptr[i] = tmp[i];
    }

    uint64_t top, r, i;
    uint64_t n = K;

    top = 1ul << (32 - __builtin_clz(K/2));
    // while (top < n - top) {top += top; }

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast, this is needed, so the
            	// sign extension needed for the mask is a unsigned one.
			    const int32_t cmp1 = compare_rows_bitonic_sort(ptr, i, i + p);
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
			            const uint32_t cmp = compare_rows_bitonic_sort(ptr, i+p, i + r);
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

/// the same as `Hoare_partition` except that we track the permutation made
/// \param V input matrix
/// \param col_l lower bound of the hoare partition
/// \param col_h upper bound of the hoare partition
/// \return new lower bound
int canonical_col_Hoare_partition(normalized_IS_t *V,
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
        if (P != NULL) { permutation_swap(P, i, j); }
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
        int p = canonical_col_Hoare_partition(V, start, end, P);
        canonical_col_lex_quicksort(V, start, p, P);
        canonical_col_lex_quicksort(V, p+1, end, P);
    }
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

///
/// \param V
/// \param G
/// \param row_l
/// \param row_h
/// \param P_r
/// \return
int row_hoare_partition(FQ_ELEM V[K][N-K],
                        normalized_IS_t *G,
                        const POSITION_T row_l,
                        const POSITION_T row_h,
                        permutation_t *P_r) {
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = V[row_l][i];
    }

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

    	// if (ret == 0) { return -1; }
        if(i >= j){ return j; }

    	row_swap(G, i, j);
    	if (P_r != NULL) { SWAP(P_r->permutation[i], P_r->permutation[j]); }
		for(uint32_t k = 0; k < N-K; k++){
			SWAP(V[i][k], V[j][k]);
		}
    }
}

/// \param V copy of G with sorted rows
/// \param G non IS generator matrix
/// \param start inclusive
/// \param end inclusive
/// \param P_r row permutation applied
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_lex_quicksort(FQ_ELEM V[K][N-K],
                      normalized_IS_t *G,
                      const uint32_t start,
                      const uint32_t end,
                      permutation_t *P_r) {
    if(start < end){
        const int p = row_hoare_partition(V, G, start, end, P_r);
    	if (p == -1) { return 0; }
        row_lex_quicksort(V, G, start, p, P_r);
        row_lex_quicksort(V, G, p+1, end, P_r);
    }

	return 1;
}

///  inplace in g
/// \param G: generator matrix to sort
/// \param P_r: row permutatution applied
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort(normalized_IS_t *G,
                   permutation_t *P_r) {
	// first sort each row into a tmp buffer
	FQ_ELEM  tmp[K][N-K];
	for (uint32_t i = 0; i < K; ++i) {
		memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
		counting_sort_u8(tmp[i], N-K);
	}

	return row_lex_quicksort(tmp, G, 0, N-K-1, P_r);
}


// NOTE: unstable sort, as we swap on equal elements
/// sort is applied inplace in G
/// \param G: generator matrix to sort
/// \param P_c: column permutatution applied
void col_bitonic_sort(normalized_IS_t *G,
                      permutation_t *P_c) {
    uint64_t r, i;
    const uint64_t n = N-K;
	const uint64_t top = 1ul << (32 - __builtin_clz(K/2));
    // while (top < n - top) { top += top; }

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
