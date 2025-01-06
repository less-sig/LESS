#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "utils.h"
#include "fq_arith.h"
#include "codes.h"
#include "transpose.h"

/// NOTE: specialised counting sort for Fq. Thus,
/// this implementation assumes that every input element
/// is reduced mod 127
/// \param arr[in/out] input array
/// \param size[in] length of the input array
void counting_sort_u8(FQ_ELEM *arr,
                      const uint32_t size) {
	/// NOTE: the type `uint32_t` is not completly arbitrary choose.
	/// Floyd did a quick benchmark between `uint16_t`, `uint32_t`, `uint64_t`
	/// and `uint32_t` seemed to be the fastest. But thats only true
	/// on a Ryzen7600X. On your machine thats maybe different.
	/// NOTE: `uint8_t` is not possible as there could be 256 times
	/// the same field element. Unlikely but possible.
	uint32_t cnt[128] __attribute__((aligned(32))) = { 0 };

    /// compute the histogram
	for (uint32_t i = 0 ; i < size ; ++i) {
		cnt[arr[i]]++;
	}

    /// compute the prefixsum
	uint32_t i = 0;
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
    while((i < (Q-1)) && (row1[i] == row2[i])) {
        i += 1;
    }
    return (((int)(row2[i]))-((int)(row1[i])));
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

/// lexicographic comparison between a row with the pivot row
/// \input: ptr[in/out]: K x (N-K) matrix
/// \input: row_idx[in]: position of the row to compare in `ptr`
/// \input: pivot[in]: pivot row
/// \returns   1 if the pivot is greater,
/// 	      -1 if it is smaller,
/// 		   0 if it matches
///
#ifdef LESS_USE_HISTOGRAM
int row_quick_sort_internal_compare_with_pivot(uint8_t *ptr[Q],
                                               const uint32_t row_idx,
                                               const uint8_t pivot[Q]){
#else 
int row_quick_sort_internal_compare_with_pivot(uint8_t *ptr[K],
                                               const uint32_t row_idx,
                                               const uint8_t pivot[K]){
#endif
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

/// NOTE: only used
int row_quick_sort_internal_hoare_partition2(FQ_ELEM* ptr[K],
                                            uint32_t P[K],
                                            const int32_t l,
                                            const int32_t h) {
    int32_t i = l - 1;
    for (int32_t j = l; j <= h - 1; j++) {
        if (compare_rows(ptr[j], ptr[h]) < 0) {
            i++;
            if (i == j) { continue; }
            SWAP(P[i], P[j]);
            cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[j]), -1ull);
        }
    }
    if (i+1 != h) {
        SWAP(P[i+1], P[h]);
        cswap((uintptr_t *)(&ptr[i+1]), (uintptr_t *)(&ptr[h]), -1ull);
    }
    return i+1;
}

void HISTEND4(uint8_t *cnt,
              uint8_t c[4][128]) {
#if 1
     for(uint32_t i = 0; i < Q; i++) {
         const uint8_t t = c[0][i]+c[1][i]+c[2][i]+c[3][i];
         ASSERT(t < Q);
         cnt[i] = t;
     }
#else
    for(uint32_t i = 0; i < Q_pad; i+=32) {
        __m256i sv =                  _mm256_load_si256((const __m256i *)&c[0][i]);
                sv = _mm256_add_epi8(_mm256_load_si256((const __m256i *)&c[1][i]), sv);
                sv = _mm256_add_epi8(_mm256_load_si256((const __m256i *)&c[2][i]), sv);
                sv = _mm256_add_epi8(_mm256_load_si256((const __m256i *)&c[3][i]), sv);
        _mm256_storeu_si256((__m256i *)&cnt[i], sv);
    }
#endif
}

void HISTEND8(uint8_t *cnt,
              uint8_t c[8][128]) {
#if 1
    for(uint32_t i = 0; i < Q; i++) {
        const uint8_t t = c[0][i]+c[1][i]+c[2][i]+c[3][i]+c[4][i]+c[5][i]+c[6][i]+c[7][i];
        ASSERT(t < Q);
        cnt[i] = t;
    }
#else
    for(uint32_t i = 0; i < Q_pad; i+=32) {
        __m256i v0 = _mm256_load_si256((const __m256i *)&c[0][i]);
        __m256i v1 = _mm256_load_si256((const __m256i *)&c[1][i]);
	    __m256i s0 = _mm256_add_epi8(v0, v1);
                v0 = _mm256_load_si256((const __m256i *)&c[2][i]);
                v1 = _mm256_load_si256((const __m256i *)&c[3][i]);
	    __m256i s1 = _mm256_add_epi8(v0, v1);
                s0 = _mm256_add_epi8(s0, s1);

                v0 = _mm256_load_si256((const __m256i *)&c[4][i]);
                v1 = _mm256_load_si256((const __m256i *)&c[5][i]);
	    		s1 = _mm256_add_epi8(v0, v1);
                v0 = _mm256_load_si256((const __m256i *)&c[6][i]);
                v1 = _mm256_load_si256((const __m256i *)&c[7][i]);
	            s0 = _mm256_add_epi8(s0, v0);
	    		s1 = _mm256_add_epi8(s1, v1);

        _mm256_storeu_si256((__m256i *)&cnt[i], _mm256_add_epi8(s0, s1));
    }
#endif
}

/// \param out[out]: pointer to the row to sort
/// \param in[in]: pointer to the row to sort
/// \param len[in]: length of the row
void row_sort(uint8_t *out,
              const uint8_t *in,
              const uint32_t len) {
#ifdef LESS_USE_HISTOGRAM
#if 1
    memset(out, 0, Q_pad);
	for (uint32_t i = 0 ; i < len ; ++i) {
	    const uint8_t t = in[i];
		out[t]++;
	}
#else
    uint8_t c[4][Q_pad] __attribute__((aligned(32))) = {0};
    const uint8_t *ip = in;
    while(ip != in+(len&~(4-1))) c[0][*ip++]++, c[1][*ip++]++, c[2][*ip++]++, c[3][*ip++]++;
    while(ip != in+ len        ) c[0][*ip++]++;
    HISTEND4(out, c);

    // uint8_t c[8][Q_pad] __attribute__((aligned(32))) = {0};
    // const uint8_t *ip = in;
    // while(ip != in+(len&~(8-1))) c[0][*ip++]++, c[1][*ip++]++, c[2][*ip++]++, c[3][*ip++]++, c[4][*ip++]++, c[5][*ip++]++, c[6][*ip++]++, c[7][*ip++]++;
    // while(ip != in+ len        ) c[0][*ip++]++;
    // HISTEND8(out, c);
#endif
#else
    memcpy(out, in, sizeof(FQ_ELEM) * N-K);
    counting_sort_u8(out, len);
#endif
}

/// internal sorting function for `row_quick_sort`
int row_quick_sort_internal(FQ_ELEM *ptr[K],
                            uint32_t P[K],
                            const uint32_t n) {
    int32_t l = 0, h = (int32_t)n-1;
    int32_t s = -1;

    // NOTE: worst case is 128
    int32_t stack[128] __attribute__((aligned(32)));
    stack[++s] = l;
    stack[++s] = h;
    while (s >= 0) {
        h = stack[s--];
        l = stack[s--];

        const int32_t p = row_quick_sort_internal_hoare_partition2(ptr, P, l, h);
        if (p - 1 > l) {
            stack[++s] = l;
            stack[++s] = p - 1;
        }
        if (p + 1 < h) {
            stack[++s] = p + 1;
            stack[++s] = h;
        }
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
	FQ_ELEM tmp[K][Q_pad] __attribute__((aligned(32)));
#else
	FQ_ELEM tmp[K][N-K];
#endif
    FQ_ELEM* ptr[K] __attribute__((aligned(32)));
    uint32_t P[K];
	for (uint32_t i = 0; i < n; ++i) {
        row_sort(tmp[i], G->values[i], N-K);

        ptr[i] = tmp[i];
        P[i] = i;
	}

   row_quick_sort_internal(ptr, P, n);

    // apply the permutation
    for (uint32_t t = 0; t < n; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        normalized_row_swap(G, t, ind);
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
                           const FQ_ELEM pivot[K]){
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

/// NOTE: only used in `col_lex_quicksort`
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

/// NOTE: only used
int row_quick_sort_internal_hoare_partition_without_histogram2(FQ_ELEM* ptr[K],
                                            uint32_t P[K],
                                            const int32_t l,
                                            const int32_t h) {
    FQ_ELEM pivot_row[N_K_pad];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = ptr[h][i];
    }

    int32_t i = l - 1;
    for (int32_t j = l; j <= h - 1; j++) {
        if (row_quick_sort_internal_compare_with_pivot_without_histogram(ptr, j, pivot_row) > 0) {
            i++;
            if (i == j) { continue; }
            SWAP(P[i], P[j]);
            cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[j]), -1ull);
        }
    }
    if (i+1 != h) {
        SWAP(P[i+1], P[h]);
        cswap((uintptr_t *)(&ptr[i+1]), (uintptr_t *)(&ptr[h]), -1ull);
    }
    return i+1;
}

/// \param ptr[in/out]:
/// \param P[in/out]: a permutation to keep track of the sorting
/// \param l[in]: inclusive
/// \param h[in]: inclusive
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_quick_sort_internal_without_histogram(FQ_ELEM* ptr[K],
                            uint32_t P[K],
                            int32_t l,
                            int32_t h) {
    int32_t s = -1;

    // NOTE: worst case is 128
    int32_t stack[64] __attribute__((aligned(32)));
    stack[++s] = l;
    stack[++s] = h;
    while (s >= 0) {
        h = stack[s--];
        l = stack[s--];

        const int32_t p = row_quick_sort_internal_hoare_partition_without_histogram2(ptr, P, l, h);
        if (p - 1 > l) {
            stack[++s] = l;
            stack[++s] = p - 1;
        }
        if (p + 1 < h) {
            stack[++s] = p + 1;
            stack[++s] = h;
        }
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
    normalized_IS_t VT = {0};
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, z, K_pad);
    for (uint32_t i = 0; i < z; i++) {
        for (uint32_t j = 0; j < K_pad; j++) {
            if (VT.values[j][i] != V->values[i][j]) {
                return;
            }
        }
    }



    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = VT.values[i];
        P[i] = i;
    }

    // row_quick_sort_recursive_internal_without_histogram(ptr, P, 0, K - 1);
    row_quick_sort_internal_without_histogram(ptr, P, 0, K - 1);

    // apply the permutation
    for (uint32_t t = 0; t < K; t++) {
        uint32_t ind = P[t];
        while(ind<t) { ind = P[ind]; }

        // TODO this swapping can be inmproved if z < K
        normalized_row_swap(&VT, t, ind);
    }

    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, K_pad, z);

    for (uint32_t i = 0; i < z; i++) {
        for (uint32_t j = 0; j < K_pad; j++) {
            if (VT.values[j][i] != V->values[i][j]) {
                return;
            }
        }
    }
}
