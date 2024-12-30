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
                    FQ_ELEM a = x[i + p];
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
	uint32_t cnt[128] __attribute__((aligned(64))) = { 0 };
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

/// NOTE: only needed for `compute_canonical_form_type4_sub`
/// \input row1[in]:
/// \input row2[in]:
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if row1 > row2
///         -x if row1 < row2
int compare_rows(const FQ_ELEM *row1,
                 const FQ_ELEM *row2) {
#ifdef LESS_USE_HISTOGRAM
    uint32_t i=0;
    while((i < (N-K-1)) && (row1[i] == row2[i])) {
        i += 1;
    }
    return -(((int)(row1[i]))-((int)(row2[i])));
#else
    uint32_t i = 0;
    while((i < (N-K)) && (row1[i] == row2[i])) {
        i += 1;
    }

    // if they are the same, they generate the same multiset
    if (i >= (N-K)) {
        return 0;
    }

    return (int)row1[i] - (int)row2[i];
#endif
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

#ifdef LESS_USE_HISTOGRAM
    uint32_t i = 0;
    while((i < (N-K-1)) && (rows[row1][i] == rows[row2][i])) {
        i += 1;
    }
    return -(((int)(rows[row1][i])) - ((int)(rows[row2][i])));

#else
    uint32_t i = 0;
    while((i < (N-K)) && (rows[row1][i] == rows[row2][i])) {
        i += 1;
    }

    // if they are the same, they generate the same multiset
    if (i >= (N-K)) {
        return 0;
    }

    return (int)rows[row1][i] - (int)rows[row2][i];
#endif
}

/// lexicographic comparison between a row with the pivot row
/// \input: ptr[in/out]: K x (N-K) matrix
/// \input: row_idx[in]: position of the row to compare in `ptr`
/// \input: pivot[in]: pivot row
/// \returns   1 if the pivot is greater,
/// 	      -1 if it is smaller,
/// 		   0 if it matches
int row_quick_sort_internal_compare_with_pivot(uint8_t *ptr[K],
                                               const uint32_t row_idx,
                                               const uint8_t pivot[K]){
#ifdef LESS_USE_HISTOGRAM
    uint32_t i=0;
    while((i<(Q-1)) && (ptr[row_idx][i]-pivot[i] == 0)){
        i++;
    }
    return ((int)ptr[row_idx][i]-(int)pivot[i]);
#else
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
#endif
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
#ifdef LESS_USE_HISTOGRAM
    FQ_ELEM pivot_row[Q];
    for(uint32_t i = 0; i < Q; i++){
       pivot_row[i] = ptr[row_l][i];
    }
#else
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = ptr[row_l][i];
    }
#endif

    POSITION_T i = row_l-1, j = row_h+1;
	int ret;
    while(1){
        do {
            i++;
        	ret = row_quick_sort_internal_compare_with_pivot(ptr, i, pivot_row);
        } while(ret > 0);

        do {
            j--;
        	ret = row_quick_sort_internal_compare_with_pivot(ptr, j, pivot_row);
        } while(ret < 0);

    	// if (ret == 0) { return -1; }
        if(i >= j){ return j; }

        SWAP(P[i], P[j]);
        cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[j]), -1ull);
    }
}

int row_quick_sort_internal_hoare_partition2(FQ_ELEM* ptr[K],
                                            uint32_t P[K],
                                            const uint32_t start,
                                            const uint32_t stop) {
    int32_t up = start, down = stop-1;
    if (stop <= start) { return start; }


#ifdef LESS_USE_HISTOGRAM
    FQ_ELEM pivot_row[Q];
    for(uint32_t i = 0; i < Q; i++){
        pivot_row[i] = ptr[down][i];
    }
#else
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = ptr[row_l][i];
    }
#endif

    while (1u) {
        while (row_quick_sort_internal_compare_with_pivot(ptr, up, pivot_row) > 0) {
            up++;
        }
        while ((row_quick_sort_internal_compare_with_pivot(ptr, down, pivot_row) > 0) && (up < down)) {
            down--;
        }
        if (up >= down) { break; }

        SWAP(P[up], P[down]);
        cswap((uintptr_t *)(&ptr[up]), (uintptr_t *)(&ptr[down]), -1ull);
        up++; down--;
    }

    SWAP(P[up], P[stop]);
    cswap((uintptr_t *)(&ptr[up]), (uintptr_t *)(&ptr[stop]), -1ull);
    return up;
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

/// \param out[out]: pointer to the row to sort
/// \param in[in]: pointer to the row to sort
/// \param len[in]: length of the row
void row_sort(uint8_t *out,
              const uint8_t *in,
              const uint32_t len) {
#ifdef LESS_USE_HISTOGRAM
    memset(out, 0, Q);
	for (uint32_t i = 0 ; i < len ; ++i) {
	    const uint8_t t = in[i];
	    ASSERT(t < Q);
		out[t]++;
	}
#else
    memcpy(out, in, sizeof(FQ_ELEM) * N-K);
    counting_sort_u8(out, len);
#endif
}

/// NOTE: only operates on ptrs
/// NOTE: not constant time
/// \param G[in/out]: generator matrix to sort
/// \param n[in] number of elements to sort
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort_recursive(normalized_IS_t *G,
                             const uint32_t n) {
    // first sort each row into a tmp buffer
#ifdef LESS_USE_HISTOGRAM
    FQ_ELEM tmp[K][Q];
#else
    FQ_ELEM tmp[K][N-K];
#endif
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < n; ++i) {
        row_sort(tmp[i], G->values[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
    }

    const int ret = row_quick_sort_internal(ptr, P, 0,  n - 1u);
    if (ret == 0) { return 0; }

    // apply the permutation
    for (uint32_t t = 0; t < n; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(G, t, ind);
    }

    return 1;
}

/// NOTE: only operates on ptrs
/// NOTE: not constant time
/// \param G[in/out]: generator matrix to sort
/// \param n[in] number of elements to sort
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort(normalized_IS_t *G,
                   const uint32_t n) {
	// first sort each row into a tmp buffer
#ifdef LESS_USE_HISTOGRAM
	FQ_ELEM tmp[K][Q];
#else
	FQ_ELEM tmp[K][N-K];
#endif
    FQ_ELEM* ptr[K];
    uint32_t P[K];
	for (uint32_t i = 0; i < n; ++i) {
        row_sort(tmp[i], G->values[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
	}

    uint32_t start = 0, stop = n-1;
    uint32_t s = 0, stack[64];
    stack[s++] = start;
    stack[s++] = stop;
    while (s > 0) {
        stop = stack[--s];
        start = stack[--s];
        if (start >= stop) {continue; }
        // i = partition(a, start, stop);
        const uint32_t i = row_quick_sort_internal_hoare_partition2(ptr, P, start, stop);
        // if (i == ) { return 0; }
        if ((i - start) > (stop - i)) {
            stack[s++] = start; stack[s++] = i - 1;
            stack[s++] = i + 1; stack[s++] = stop;
        } else {
            stack[s++] = i + 1; stack[s++] = stop;
            stack[s++] = start; stack[s++] = i - 1;
        }
    }

    // apply the permutation
    for (uint32_t t = 0; t < n; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(G, t, ind);
    }

    return 1;
}


/// NOTE: operates on pointers
/// \input G[in/out]: normalised non IS part of a generator matrix
/// \return the sorting algorithm works only inplace for the sorting of the columns
/// 		0 on failure: row_i and row_j generate the same multiset
/// 		1 on success
int row_bitonic_sort(normalized_IS_t *G) {
    // first sort each row into a tmp buffer
#ifdef LESS_USE_HISTOGRAM
	FQ_ELEM  tmp[K][Q];
#else
	FQ_ELEM  tmp[K][N-K];
#endif
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        row_sort(tmp[i], G->values[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
    }

    const uint64_t n = K;

    const uint64_t top = 1ul << (32 - __builtin_clz(K/2));

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (uint64_t i = 0; i < n - p; ++i) {
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
            for (uint64_t i = 0; i < n - q; ++i) {
                if (!(i & p)) {
                    for (uint64_t r = q; r > p; r >>= 1) {
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




/// lexicographic comparison
/// \return G1[col1] <=> G2[col2]:
///         -1: G1[col1] > G2[col2]
///          0: G1[col1] == G2[col2]
///          1: G1[col1] < G2[col2]
int lex_compare_column(const generator_mat_t *G1,
					   const generator_mat_t *G2,
                       const POSITION_T col1,
                       const POSITION_T col2) {
   uint32_t i=0;
   while((i < K) &&
         (G1->values[i][col1]-G2->values[i][col2] == 0)) {
       i++;
   }

   if (i >= K) return 0;

   if (G1->values[i][col1]-G2->values[i][col2] > 0){
      return -1;
   }

   return 1;
}

/// lexicographic comparison
/// \return G1[col1] <=> G1[col2]:
///           1: G1[col1] >  G1[col2]
///           0: G1[col1] == G1[col2]
///          -1: G1[col1] <  G1[col2]
int lex_compare_col(const normalized_IS_t *G1,
                    const POSITION_T col1,
                    const POSITION_T col2) {
   uint32_t i=0;
   while((i < (K-1)) &&
         (G1->values[i][col1]-G1->values[i][col2] == 0)) {
       i++;
   }
   return G1->values[i][col1]-G1->values[i][col2];
}

/* lexicographic comparison of a column with the pivot
 * returns 1 if the pivot is greater, -1 if it is smaller,
 * 0 if it matches */
int lex_compare_with_pivot(normalized_IS_t *V,
                           const POSITION_T col_idx,
                           FQ_ELEM pivot[K]){
   uint32_t i=0;
   while(i<K && V->values[i][col_idx]-pivot[i] == 0){
       i++;
   }
   if (i==K) return 0;
   if (V->values[i][col_idx]-pivot[i] > 0){
      return -1;
   }
   return 1;
}

///
int Hoare_partition(normalized_IS_t *V,
                    const POSITION_T col_l,
                    const POSITION_T col_h){
    FQ_ELEM pivot_col[K];
    for(uint32_t i = 0; i < K; i++){
       pivot_col[i] = V->values[i][col_l];
    }
    // TODO double comparison

    POSITION_T i = col_l, j = col_h+1;
    do {
        j--;
    } while(lex_compare_with_pivot(V,j,pivot_col) == -1);
    if(i >= j){
        return j;
    }

    column_swap(V,i,j);

    while(1){
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
    }
}

/* In-place quicksort */
void col_lex_quicksort(normalized_IS_t *V,
                       int start,
                       int end){
    if(start < end){
        int p = Hoare_partition(V,start,end);
        col_lex_quicksort(V,start,p);
        col_lex_quicksort(V,p+1,end);
    }
}

/* Sorts the columns of V in lexicographic order */
void lex_sort_cols(normalized_IS_t *V){
   col_lex_quicksort(V,0,(N-K)-1);
}






/// lexicographic comparison between a row with the pivot row
/// \input: ptr[in/out]: K x (N-K) matrix
/// \input: row_idx[in]: position of the row to compare in `ptr`
/// \input: pivot[in]: pivot row
/// \returns   1 if the pivot is greater,
/// 	      -1 if it is smaller,
/// 		   0 if it matches
int row_quick_sort_internal_compare_with_pivot_without_histogram(uint8_t *ptr[K],
                                                            const POSITION_T row_idx,
                                               const uint8_t pivot[K]){
    uint32_t i=0;
    while((i < (N-K-1)) && (ptr[row_idx][i]-pivot[i] == 0)){
        i++;
    }
    return -(((int)(ptr[row_idx][i]))-((int)(pivot[i])));
}

/// TODO doc
/// \param ptr
/// \param P permutation: to keep track of the sorting
/// \param row_l
/// \param row_h
/// \return
int row_quick_sort_internal_hoare_partition_without_histogram(FQ_ELEM* ptr[K],
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
        	ret = row_quick_sort_internal_compare_with_pivot_without_histogram(ptr, i, pivot_row);
        } while(ret > 0);

        do {
            j--;
        	ret = row_quick_sort_internal_compare_with_pivot_without_histogram(ptr, j, pivot_row);
        } while(ret < 0);

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
int row_quick_sort_internal_without_histogram(FQ_ELEM* ptr[K],
                            uint32_t P[K],
                            const uint32_t start,
                            const uint32_t end) {
    if(start < end){
        const int p = row_quick_sort_internal_hoare_partition_without_histogram(ptr, P, start, end);
    	if (p == -1) { return 0; }
        row_quick_sort_internal_without_histogram(ptr, P, start, p);
        row_quick_sort_internal_without_histogram(ptr, P, p + 1, end);
    }

	return 1;
}


/// NOTE: non-constant time
/// Sorts the columns of the input matrix, via first transposing
/// the matrix, subsequent sorting rows, and finally transposing
/// it back.
/// \param V[in/out]: non IS-part of a generator matrix
/// \param z[in]: number of rows within each col to sort
void col_quicksort_transpose(normalized_IS_t *V,
                             const uint32_t z) {
    normalized_IS_t VT;
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, K, z);

    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = VT.values[i];
        P[i] = i;
    }

    row_quick_sort_internal_without_histogram(ptr, P, 0, N - K - 1);

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        row_swap(&VT, t, ind);
    }

    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, z, K);
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
            	// sign extension needed for the mask is an unsigned one.
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
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, K, K);
	col_bitonic_sort_transposed(&VT);
    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, K, K);
}







/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ah still TODO
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/// lexicographic comparison of a column with the pivot
/// returns 1 if the pivot is greater, -1 if it is smaller,
/// 0 if it matches */
int col_quicksort_internal_compare_with_pivot(FQ_ELEM *V[K],
                                              const POSITION_T col_idx,
                                              const FQ_ELEM pivot[K],
                                              const uint32_t z) {
    uint32_t i=0;
    printf("%d\n", col_idx);
    while(i<z && V[col_idx][i*K]-pivot[i] == 0){
        i++;
    }
    if (i==z) return 0;

     /// TODO: optimize the
    if (V[col_idx][i*K]-pivot[i] > 0){
       return -1;
    }

    return 1;
}

void col_quicksort_internal_column_swap(normalized_IS_t *V,
                                        const POSITION_T col1,
                                        const POSITION_T col2,
                                        const uint32_t z){
   for(uint32_t i = 0; i<z;i++ ){
      const POSITION_T tmp = V->values[i][col2];
      V->values[i][col2] = V->values[i][col1];
      V->values[i][col1] = tmp;
   }
}
///
uint32_t col_quicksort_internal_hoare_partition(FQ_ELEM *V[K],
                                                uint32_t P[K],
                                                const POSITION_T col_l,
                                                const POSITION_T col_h,
                                                const uint32_t z){
    FQ_ELEM pivot_col[K];
    for(uint32_t i = 0; i < z; i++){
       pivot_col[i] = V[col_l][i*K];
    }

    POSITION_T i = col_l, j = col_h+1;
    do {
        j--;
    } while(col_quicksort_internal_compare_with_pivot(V, j, pivot_col, z) == -1);

    if(i >= j){
        return j;
    }

    // column_swap(V,i,j);
    SWAP(P[i], P[j]);
    cswap((uintptr_t *)V[i], (uintptr_t *)V[j], -1ull);

    while(1){
        do {
            i++;
        } while(col_quicksort_internal_compare_with_pivot(V, i, pivot_col, z) == 1);
        do {
            j--;
        } while(col_quicksort_internal_compare_with_pivot(V, j, pivot_col, z) == -1);

        if(i >= j){
            return j;
        }

        // column_swap(V,i,j);
        SWAP(P[i], P[j]);
        cswap((uintptr_t *)V[i], (uintptr_t *)V[j], -1ull);
    }
}

/// TODO doc and unfinished
/// \param V
/// \param P
/// \param start
/// \param end
/// \param z
void col_quicksort_internal(FQ_ELEM *V[K],
                            uint32_t P[K],
                            const uint32_t start,
                            const uint32_t end,
                            const uint32_t z) {
    if(start < end){
        const uint32_t p = col_quicksort_internal_hoare_partition(V, P, start, end, z);
        col_quicksort_internal(V, P, start, p, z);
        col_quicksort_internal(V, P, p+1, end, z);
    }
}

/// \param V[in/out]:
/// \param z[in]: number of rows in each column
void col_quicksort(normalized_IS_t *V,
                   const uint32_t z) {
    uint32_t P[K];
    FQ_ELEM* ptr[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = &(V->values[0][i]);
        P[i] = i;
    }

    col_quicksort_internal(ptr, P, 0, K-1, z);

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        column_swap(V, t, ind);
    }
}
