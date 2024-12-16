#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "utils.h"
#include "fq_arith.h"
#include "codes.h"
#include "transpose.h"

/// taken from djbsort
/// \param a[in/out] first input
/// \param b[in/out] second input
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

/// NOTE: taken from djbsort
/// \param x[in/out] input array
/// \param n[in] length
void bitonic_sort_i8(FQ_ELEM *x,
                     const long long n) {
    long long p, q, r, i;
    if (n < 2) {
        return;
    }
    const long long top = 1ul << (32 - __builtin_clz(K/2));

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

/// NOTE: specialised counting sort for Fq. Thus,
/// this implementation assumes that every input element
/// is reduced mod 127
/// \param arr[in/out] input array
/// \param size[in] length of the input array
void counting_sort_u8(FQ_ELEM *arr,
                      const size_t size) {
	/// NOTE: the type `uint32_t` is not completly arbitrary choose.
	/// Floyd did a quick benchmark between `uint16_t`, `uint32_t`, `uint64_t`
	/// and `uint32_t` seemed to be the fastest. But thats only true
	/// on a Ryzen7600X. On your machine thats maybe different.
	/// NOTE: `uint8_t` is not possible as there could be 256 times
	/// the same field element. Unlikely but possible.
	uint32_t cnt[128] __attribute__((aligned(512))) = { 0 };
	size_t i;

    /// compute the histogram
	for (i = 0 ; i < size ; ++i) {
		cnt[arr[i]]++;
	}

    /// compute the prefixsum
	i = 0;
	for (size_t a = 0 ; a < Q; ++a) {
		while (cnt[a]--) {
			arr[i++] = a;
		}
	}
}

/// \input row1[in]:
/// \input row2[in]:
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
int compare_rows(const FQ_ELEM *row1,
                 const FQ_ELEM *row2) {
    uint32_t i = 0; 
    while((i < (N-K)) && (row1[i] == row2[i])) {
        i += 1;
    }

    // if they are the same, they generate the same multiset
    if (i >= (N-K)) {
        return 0;
    }

    return (int)row1[i] - (int)row2[i];
}


/// NOTE: helper function for type3 canonical form.
/// NOTE: specially made for the bitonic sort, which
///		operates on pointers.
/// \input: rows[in/out]: K x (N-K) matrix
/// \input: row1[in]: first row to compare
/// \input: row2[in]: secnd row to compare
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
int compare_rows_bitonic_sort(FQ_ELEM **rows,
							  const uint32_t row1,
							  const uint32_t row2) {
    ASSERT(row1 < K);
    ASSERT(row2 < K);

    uint32_t i = 0;
    while((i < (N-K)) && (rows[row1][i] == rows[row2][i])) {
        i += 1;
    }

    // if they are the same, they generate the same multiset
    if (i >= (N-K)) {
        return 0;
    }

    return (int)rows[row1][i] - (int)rows[row2][i];
}

/// lexicographic comparison between a row with the pivot row
/// \input: ptr[in/out]: K x (N-K) matrix
/// \input: row_idx[in]: position of the row to compare in `ptr`
/// \input: pivot[in]: pivot row
/// \returns   1 if the pivot is greater,
/// 	      -1 if it is smaller,
/// 		   0 if it matches
int row_quick_sort_internal_compare_with_pivot(uint8_t *ptr[K],
                                               const POSITION_T row_idx,
                                               const uint8_t pivot[K]){
    uint32_t i=0;
    while((i<(N-K)) && (ptr[row_idx][i]-pivot[i] == 0)){
        i++;
    }
    if (i==(N-K)) {
        return 0;
    }

    if ((int)ptr[row_idx][i]-(int)pivot[i] > 0){
        return -1;
    }

    return 1;
}

/// TODO doc
/// \param ptr
/// \param P permutation: to keep track of the sorting
/// \param row_l
/// \param row_h
/// \return
int row_quick_sort_internal_hoare_partition(FQ_ELEM* ptr[K],
                                            uint32_t P[K],
                                            const POSITION_T row_l,
                                            const POSITION_T row_h) {
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = ptr[row_l][i];
    }

    POSITION_T i = row_l-1, j = row_h+1;
	int ret;
    while(1){
        do {
            i++;
        	ret = row_quick_sort_internal_compare_with_pivot(ptr, i, pivot_row);
        } while(ret == 1);

        do {
            j--;
        	ret = row_quick_sort_internal_compare_with_pivot(ptr, j, pivot_row);
        } while(ret == -1);

    	// if (ret == 0) { return -1; }
        if(i >= j){ return j; }

        SWAP(P[i], P[j]);
        cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[j]), -1ull);
    }
}

/// \param ptr[in/out]:
/// \param P[in/out]: a permutation to keep track of the sorting
/// \param start[in]: inclusive
/// \param end[in]: inclusive
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_quick_sort_internal(FQ_ELEM* ptr[K],
                            uint32_t P[K],
                            const uint32_t start,
                            const uint32_t end) {
    if(start < end){
        const int p = row_quick_sort_internal_hoare_partition(ptr, P, start, end);
    	if (p == -1) { return 0; }
        row_quick_sort_internal(ptr, P, start, p);
        row_quick_sort_internal(ptr, P, p + 1, end);
    }

	return 1;
}

/// \param ptr[in/out]: pointer to the row to sort
/// \param len[in]: length of the row
void row_sort(uint8_t *ptr, 
              const uint32_t len) {
    counting_sort_u8(ptr, len);
}

/// NOTE: only operates on ptrs
/// NOTE: not constant time
/// \param G[in/out]: generator matrix to sort
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort(normalized_IS_t *G) {
	// first sort each row into a tmp buffer
	FQ_ELEM  tmp[K][N-K];
    FQ_ELEM* ptr[K];
    uint32_t P[K];
	for (uint32_t i = 0; i < K; ++i) {
		memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
        // TODO maybe just histogram?
        row_sort(tmp[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
	}

	const int ret = row_quick_sort_internal(ptr, P, 0, N - K - 1);
    if (ret == 0) { return 0; }

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(G, t, ind);
    }

    return 1;
}


/// NOTE: non constant time
/// Sorts the columns of the input matrix, via first transposing
/// the matrix, subsequent sorting rows, and finally transposing
/// it back.
/// \param V[in/out]: non IS-part of a generator matrix
void col_quicksort_transpose(normalized_IS_t *V) {
    normalized_IS_t VT;
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, K);

    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = VT.values[i];
        P[i] = i;
    }

    row_quick_sort_internal(ptr, P, 0, N - K - 1);

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(&VT, t, ind);
    }

    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, K);
}

/// NOTE: internal function, do not call it directly.
/// NOTE: constant time version
/// NOTE: sort pointers, and applies the final resulting permutation
///     afterward to the generator matrix.
/// Sorts the columns of the input matrix, via first transposing
/// the matrix, subsequent sorting rows, and finally transposing
/// it back.
/// \input G[in/out]: normalised non IS part of a generator matrix
int col_bitonic_sort_transposed(normalized_IS_t *G) {
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = G->values[i];
    	P[i] = i;
    }

    uint64_t n = K;

    const uint64_t top = 1ul << (32 - __builtin_clz(K/2));

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (uint64_t i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast, this is needed, so the
            	// sign extension needed for the mask is a unsigned one.
			    const uint32_t cmp = compare_rows_bitonic_sort(ptr, i, i + p);

                const uintptr_t mask = -(1ull - (cmp >> 31));
                cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[i+p]), mask);
				MASKED_SWAP(P[i], P[i+p], mask);
            }
        }

        for (uint64_t q = top; q > p; q >>= 1) {
            for (uint64_t i = 0; i < n - q; ++i) {
                if (!(i & p)) {
                    for (uint64_t r = q; r > p; r >>= 1) {
			            const uint32_t cmp = compare_rows_bitonic_sort(ptr, i+p, i + r);
                        const uintptr_t mask = -(1ull - (cmp >> 31));
                        cswap((uintptr_t *)(&ptr[i+p]), (uintptr_t *)(&ptr[i+r]), mask);
						MASKED_SWAP(P[i+p], P[i+r], mask);
                    }
                }
            }
        }
    }

	// apply the permutation
	for (uint32_t t = 0; t < K; t++) {
		uint32_t ind = P[t];
		while(ind<t) { ind = P[ind]; }

		row_swap(G, t, ind);
	}

    return 1;
}

/// NOTE: constant time version
/// NOTE: sort pointers, and applies the final resulting permutation
///     afterward to the generator matrix.
/// Sorts the columns of the input matrix, via first transposing
/// the matrix, subsequent sorting rows, and finally transposing
/// it back.
/// \input G[in/out]: normalised non IS part of a generator matrix
void col_bitonic_sort_transpose(normalized_IS_t *V) {
    normalized_IS_t VT;
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, K);
	col_bitonic_sort_transposed(&VT);
    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, K);
}

/// NOTE: operates on pointers
/// \input G[in/out]: normalised non IS part of a generator matrix
/// \return the sorting algorithm works only inplace for the sorting of the columns
/// 		0 on failure: row_i and row_j generate the same multiset
/// 		1 on success
int row_bitonic_sort(normalized_IS_t *G) {
    // first sort each row into a tmp buffer
    FQ_ELEM  tmp[K][N-K];
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
        row_sort(tmp[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
    }

    uint64_t r, i;
    uint64_t n = K;

    const uint64_t top = 1ul << (32 - __builtin_clz(K/2));

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast, this is needed, so the
                // sign extension needed for the mask is an unsigned one.
                const int32_t cmp1 = compare_rows_bitonic_sort(ptr, i, i + p);
                const uint32_t cmp = cmp1;
                if (cmp == 0) { return 0; }

                const uintptr_t mask = -(1ull - (cmp >> 31));
                cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[i+p]), mask);
                MASKED_SWAP(P[i], P[i+p], mask);
            }
        }

        for (uint64_t q = top; q > p; q >>= 1) {
            for (i = 0; i < n - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1) {
                        const uint32_t cmp = compare_rows_bitonic_sort(ptr, i+p, i + r);
                        if (cmp == 0) { return 0; }

                        const uintptr_t mask = -(1ull - (cmp >> 31));
                        cswap((uintptr_t *)(&ptr[i+p]), (uintptr_t *)(&ptr[i+r]), mask);
                        MASKED_SWAP(P[i+p], P[i+r], mask);
                    }
                }
            }
        }
    }

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(G, t, ind);
    }

    return 1;
}
