#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "utils.h"
#include "fq_arith.h"
#include "codes.h"
#include "transpose.h"

// TODO move to opt impl
#ifdef USE_AVX2
#include <immintrin.h>


#define COEX_u8x16(a, b, tmp)             \
	{                                     \
		tmp = a;                 		  \
		a = _mm_min_epu8(a, b);        	  \
		b = _mm_max_epu8(tmp, b);  	      \
	}
#define COEX_u8x32(a, b, tmp)             \
	{                                     \
		tmp = a;                 		  \
		a = _mm256_min_epu8(a, b);        \
		b = _mm256_max_epu8(tmp, b);      \
	}

// TODO remove
static int8_t blend[4][16] = {
        {0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1},
        {0, 0,-1,-1,0,0,-1,-1,0,0,-1,-1,0,0,-1,-1},
        {0,0,0,0,-1,-1,-1,-1,0,0,0,0,-1,-1,-1,-1},
        {0,0,0,0,0,0,0,0, -1,-1,-1,-1,-1,-1,-1,-1},
};
static int8_t sortingnetwork_u8x32_shuffle_masks[6][32] __attribute((aligned(64))) = {
        {1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14},
        {3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12,3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12},
        {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8},
        {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
        {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0},
        {4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11},
};

__m256i sortingnetwork_sort_u8x32_(__m256i v) {
    __m256i t = v, tmp;
    t = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]));
    COEX_u8x32(t, v, tmp);
    t = _mm256_blendv_epi8(t, v, _mm256_set1_epi16(0xFF));

    // Step 2
    v = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[1]));
	COEX_u8x32(v, t, tmp);
    v = _mm256_blendv_epi8(v, t, _mm256_set1_epi32(0xFFFF));

    // Step 3
    t = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]));
	COEX_u8x32(t, v, tmp);
    t = _mm256_blendv_epi8(t, v, _mm256_set1_epi16(0xFF));

    // Step 4
    v = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[2]));
	COEX_u8x32(v, t, tmp);
    v = _mm256_blendv_epi8(v, t, _mm256_set1_epi64x(0xFFFFFFFF));

    // Step 5
    t = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[3]));
	COEX_u8x32(t, v, tmp);
    t = _mm256_blendv_epi8(t, v, _mm256_set1_epi32(0xFFFF));

    // Step 6
    v = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]));
	COEX_u8x32(v, t, tmp);
    v = _mm256_blendv_epi8(v, t, _mm256_set1_epi16(0xFF));

    // Step 7
    t = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[4]));
	COEX_u8x32(t, v, tmp);
    t = _mm256_blendv_epi8(t, v, _mm256_setr_epi64x(0, -1ull, 0, -1ull));

    // Step 8
    v = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[5]));
	COEX_u8x32(v, t, tmp);
    v = _mm256_blendv_epi8(v, t, _mm256_set1_epi64x(0xFFFFFFFF));

    // Step 9
    t = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[3]));
	COEX_u8x32(t, v, tmp);
    t = _mm256_blendv_epi8(t, v, _mm256_set1_epi32(0xFFFF));

    // Step 10
    v = _mm256_shuffle_epi8(t, _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]));
	COEX_u8x32(v, t, tmp);
    v = _mm256_blendv_epi8(v, t, _mm256_set1_epi16(0xFF));

    __m128 tmp_;
    __m128i L1 = _mm256_extractf128_si256(v, 0);
    __m128i H1 = _mm256_extractf128_si256(v, 1);
    H1 = _mm_shuffle_epi8(H1, _mm_load_si128((__m128i *)sortingnetwork_u8x32_shuffle_masks[4]));

	COEX_u8x16(L1, H1, tmp_);
    __m128i L1p = _mm_blendv_epi8(L1, _mm_bslli_si128(H1, 8), _mm_load_si128((__m128i *)blend[3]));
    __m128i H1p = _mm_blendv_epi8(_mm_bsrli_si128(L1, 8), H1, _mm_load_si128((__m128i *)blend[3]));

	COEX_u8x16(L1p, H1p, tmp_);
    __m128i L2p = _mm_blendv_epi8(L1p, _mm_bslli_si128(H1p, 4), _mm_load_si128((__m128i *)blend[2]));
    __m128i H2p = _mm_blendv_epi8(_mm_bsrli_si128(L1p, 4), H1p, _mm_load_si128((__m128i *)blend[2]));

	COEX_u8x16(L2p, H2p, tmp_);
    __m128i L3p = _mm_blendv_epi8(L2p, _mm_bslli_si128(H2p, 2), _mm_load_si128((__m128i *)blend[1]));
    __m128i H3p = _mm_blendv_epi8(_mm_bsrli_si128(L2p, 2), H2p, _mm_load_si128((__m128i *)blend[1]));

	COEX_u8x16(L3p, H3p, tmp_);
    __m128i L4p = _mm_blendv_epi8(L3p, _mm_bslli_si128(H3p, 1), _mm_load_si128((__m128i *)blend[0]));
    __m128i H4p = _mm_blendv_epi8(_mm_bsrli_si128(L3p, 1), H3p, _mm_load_si128((__m128i *)blend[0]));

	COEX_u8x16(L4p, H4p, tmp_);
    const __m128i kl = _mm_unpacklo_epi8(L4p, H4p);
    const __m128i kh = _mm_unpackhi_epi8(L4p, H4p);
    return _mm256_set_m128i(kh, kl);
}


/// implementation of 8 parrallel `simd_aftermerge_1V`
static inline void sortingnetwork_aftermergesort_u8x64(__m256i *a,
													   __m256i *b) {
	__m256i L0 = *a;
	__m256i H0 = *b;
	__m256i mask, tmp;

	// 3
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[2]);
	__m256i L3p = _mm256_shuffle_epi8(L0, mask);
	__m256i H3p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L3p, tmp);
    COEX_u8x32(H0, H3p, tmp);
	mask = _mm256_set1_epi64x(0xFFFFFFFF);
	L0 = _mm256_blendv_epi8(L3p, L0, mask);
	H0 = _mm256_blendv_epi8(H3p, H0, mask);

	// 4
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[1]);
	__m256i L4p = _mm256_shuffle_epi8(L0, mask);
	__m256i H4p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L4p, tmp);
    COEX_u8x32(H0, H4p, tmp);
	mask = _mm256_set1_epi32(0xFFFF);
	L0 = _mm256_blendv_epi8(L4p, L0, mask);
	H0 = _mm256_blendv_epi8(H4p, H0, mask);

	// 5
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]);
	__m256i L5p = _mm256_shuffle_epi8(L0, mask);
	__m256i H5p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L5p, tmp);
    COEX_u8x32(H0, H5p, tmp);
	mask = _mm256_set1_epi16(0xFF);
	*a = _mm256_blendv_epi8(L5p, L0, mask);
	*b = _mm256_blendv_epi8(H5p, H0, mask);
}
static inline void sortingnetwork_mergesort_u8x64(__m256i *a,
                                                  __m256i *b) {
    __m256i L0 = *a, tmp, mask;
    __m256i H0 = *b;

	// reverse H0
    __m256i H0p = _mm256_permute2x128_si256(H0, H0, 0b00000001);
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[4]);
	H0 =_mm256_shuffle_epi8(H0p, mask);

    // 0
    COEX_u8x32(L0, H0, tmp);

	//1
    __m256i L1p = _mm256_permute2x128_si256(L0, L0, 0b00000001);
    __m256i H1p = _mm256_permute2x128_si256(H0, H0, 0b00000001);
    COEX_u8x32(L0, L1p, tmp);
    COEX_u8x32(H0, H1p, tmp);
	L0 = _mm256_blend_epi32(L0, L1p, 0b11110000);
	H0 = _mm256_blend_epi32(H0, H1p, 0b11110000);

	// 2
	__m256i L2p = _mm256_permute4x64_epi64(L0, 0b10110001);
	__m256i H2p = _mm256_permute4x64_epi64(H0, 0b10110001);
    COEX_u8x32(L0, L2p, tmp);
    COEX_u8x32(H0, H2p, tmp);
	*a = _mm256_blend_epi32(L0, L2p, 0b11001100);
	*b = _mm256_blend_epi32(H0, H2p, 0b11001100);

	sortingnetwork_aftermergesort_u8x64(a, b);
}

static inline void sortingnetwork_sort_u8x64(__m256i *a,
						                     __m256i *b) {
    *a = sortingnetwork_sort_u8x32_(*a);
    *b = sortingnetwork_sort_u8x32_(*b);
    sortingnetwork_mergesort_u8x64(a, b);
}

/// implementation of `simd_aftermerge_8V`
static inline void sortingnetwork_aftermerge_u8x64(__m256i *a,
												   __m256i *b) {
	__m256i tmp;
    COEX_u8x32(*a, *b, tmp);

	__m256i ap = _mm256_permute4x64_epi64(*a, 0b01001110);
	__m256i bp = _mm256_permute4x64_epi64(*b, 0b01001110);
    COEX_u8x32(*a, ap, tmp);
    COEX_u8x32(*b, bp, tmp);
	*a = _mm256_blend_epi32(*a, ap, 0b11110000);
	*b = _mm256_blend_epi32(*b, bp, 0b11110000);

	ap = _mm256_permute4x64_epi64(*a, 0b10110001);
	bp = _mm256_permute4x64_epi64(*b, 0b10110001);
    COEX_u8x32(*a, ap, tmp);
    COEX_u8x32(*b, bp, tmp);
	*a = _mm256_blend_epi32(*a, ap, 0b11001100);
	*b = _mm256_blend_epi32(*b, bp, 0b11001100);

	sortingnetwork_aftermergesort_u8x64(a, b);
}
static inline void sortingnetwork_sort_u8x128(__m256i *a,
                                              __m256i *b,
                                              __m256i *c,
                                              __m256i *d) {
	__m256i tmp, mask;
	sortingnetwork_sort_u8x64(a, b);
	sortingnetwork_sort_u8x64(c, d);

	// reverse c and d
    __m256i cp = _mm256_permute2x128_si256(*c, *c, 0b00000001);
    __m256i dp = _mm256_permute2x128_si256(*d, *d, 0b00000001);
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[4]);
	cp =_mm256_shuffle_epi8(cp, mask);
    COEX_u8x32(*b, cp, tmp);
	dp =_mm256_shuffle_epi8(dp, mask);
    COEX_u8x32(*a, dp, tmp);

	*c = cp;
	*d = dp;
	sortingnetwork_aftermerge_u8x64(a, b);
	sortingnetwork_aftermerge_u8x64(c, d);
}

static inline void sortingnetwork(FQ_ELEM *arr,
                                  const size_t s) {

    __m256i i1 = _mm256_loadu_si256((const __m256i *)(arr +  0));
    __m256i i2 = _mm256_loadu_si256((const __m256i *)(arr + 32));
    __m256i i3 = _mm256_loadu_si256((const __m256i *)(arr + 64));
    __m256i i4 = _mm256_loadu_si256((const __m256i *)(arr + 96));
    sortingnetwork_sort_u8x128(&i1, &i2, &i3, &i4);
    _mm256_storeu_si256((__m256i *)(arr +  0), i1);
    _mm256_storeu_si256((__m256i *)(arr + 32), i2);
    _mm256_storeu_si256((__m256i *)(arr + 64), i3);
    _mm256_storeu_si256((__m256i *)(arr + 96), i4);
}
#endif

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
	/// the same field element. Unlikely but possible.
	uint32_t cnt[128] __attribute__((aligned(512))) = { 0 };
	size_t i;

	for (i = 0 ; i < size ; ++i) {
		cnt[arr[i]]++;
	}

	i = 0;
	for (size_t a = 0 ; a < Q; ++a) {
		while (cnt[a]--) {
			arr[i++] = a;
		}
	}
}

/// @param arr
/// @param in
/// @param size
void counting_sort_u8_v2(FQ_ELEM *arr,
                         const FQ_ELEM *in,
                         const size_t size) {
	/// NOTE: the type `uint32_t` is not completly arbitrary choose.
	/// Floyd did a quick benchmark between `uint16_t`, `uint32_t`, `uint64_t`
	/// and `uint32_t` seemed to be the fastest. But thats only true
	/// on a Ryzen7600X. On your machine thats maybe different.
	/// NOTE: `uint8_t` is not possible as there could be 256 times
	/// the same field element. Unlikely but possible.
	uint32_t cnt[128] __attribute__((aligned(512))) = { 0 };
	size_t i;

	for (i = 0 ; i < size ; ++i) {
		cnt[in[i]]++;
	}

	i = 0;
	for (size_t a = 0 ; a < Q; ++a) {
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
///         -x: b > a;
///          0: b == a;
///          x: b < a;
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
int row_bubble_sort(normalized_IS_t *G) {
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

/// TODO merge with  `void canonical_col_lex_quicksort_transpose`
/// \input normalised non IS part of a generator matrix
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
        counting_sort_u8(tmp[i], N-K);
    	//sortingnetwork(tmp[i], N-K);

        // counting_sort_u8_v2(tmp[i], G->values, N-K);
        ptr[i] = tmp[i];
    	P[i] = i;
    }

    uint64_t top, r, i;
    uint64_t n = K;

    top = 1ul << (32 - __builtin_clz(K/2));

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
				MASKED_SWAP(P[i], P[i+p], mask);
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

/// the same as `Hoare_partition` except that we track the permutation made
/// \param V input matrix
/// \param col_l lower bound of the hoare partition
/// \param col_h upper bound of the hoare partition
/// \return new lower bound
int canonical_col_Hoare_partition(normalized_IS_t *V,
                                  const POSITION_T col_l,
								  const POSITION_T col_h) {
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
    }
}

/// TODO remove
/// In-place quicksort
/// the same as `col_lex_quicksort` except we are tracking the permutations
/// \param V
/// \param start
/// \param end
void canonical_col_lex_quicksort(normalized_IS_t *V,
                                 const int start,
                                 const int end) {
    if(start < end){
        int p = canonical_col_Hoare_partition(V, start, end);
        canonical_col_lex_quicksort(V, start, p);
        canonical_col_lex_quicksort(V, p+1, end);
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
/// \return
int row_hoare_partition(FQ_ELEM V[K][N-K],
                        normalized_IS_t *G,
                        const POSITION_T row_l,
                        const POSITION_T row_h) {
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
		for(uint32_t k = 0; k < N-K; k++){
			SWAP(V[i][k], V[j][k]);
		}
    }
}

/// \param V copy of G with sorted rows
/// \param G non IS generator matrix
/// \param start inclusive
/// \param end inclusive
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_lex_quicksort(FQ_ELEM V[K][N-K],
                      normalized_IS_t *G,
                      const uint32_t start,
                      const uint32_t end) {
    if(start < end){
        const int p = row_hoare_partition(V, G, start, end);
    	if (p == -1) { return 0; }
        row_lex_quicksort(V, G, start, p);
        row_lex_quicksort(V, G, p+1, end);
    }

	return 1;
}

///  inplace in g
/// \param G: generator matrix to sort
/// \return 1 on success
///			0 if two rows generate the same multiset
int row_quick_sort(normalized_IS_t *G) {
	// first sort each row into a tmp buffer
	FQ_ELEM  tmp[K][N-K];
	for (uint32_t i = 0; i < K; ++i) {
		counting_sort_u8_v2(tmp[i], G->values[i], N-K);

		// memcpy(tmp[i], G->values[i], sizeof(FQ_ELEM) * N-K);
		// counting_sort_u8(tmp[i], N-K);
		// qsort(tmp[i], N-K, 1, fqcmp);
	}

	return row_lex_quicksort(tmp, G, 0, N-K-1);
}


// NOTE: unstable sort, as we swap on equal elements
/// sort is applied inplace in G
/// \param G: generator matrix to sort
void col_bitonic_sort(normalized_IS_t *G) {
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
            }
        }

        i = 0;
        for (uint64_t q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1) {
						const uint32_t cmp = lex_compare_col(G, i+p, i+r);
            			const uintptr_t mask = -(1ull - (cmp >> 31));
            			column_cswap(G, i+p, i+r, mask);
                    }
                }
            }
        }
    }
}

/// \param G
/// \param row_l
/// \param row_h
/// \return
int row_hoare_partition_transposed(normalized_IS_t *G,
                                   const POSITION_T row_l,
                                   const POSITION_T row_h) {
    FQ_ELEM pivot_row[N-K];
    for(uint32_t i = 0; i < N-K; i++){
       pivot_row[i] = G->values[row_l][i];
    }

    POSITION_T i = row_l-1, j = row_h+1;
	int ret;
    while(1){
        do {
            i++;
        	ret = row_lex_compare_with_pivot(G->values, i, pivot_row);
        } while(ret == 1);

        do {
            j--;
        	ret = row_lex_compare_with_pivot(G->values, j, pivot_row);
        } while(ret == -1);

    	// if (ret == 0) { return -1; }
        if(i >= j){ return j; }

    	row_swap(G, i, j);
    }
}

/// \param G non IS generator matrix
/// \param start inclusive
/// \param end inclusive
/// \return 1 on success
///			0 if two rows generate the same multi set
int row_lex_quicksort_transposed(normalized_IS_t *G,
                                 const uint32_t start,
                                 const uint32_t end) {
    if(start < end){
        const int p = row_hoare_partition_transposed(G, start, end);
    	if (p == -1) { return 0; }
        row_lex_quicksort_transposed(G, start, p);
        row_lex_quicksort_transposed(G, p+1, end);
    }

	return 1;
}


/// \input normalised non IS part of a generator matrix
int col_bitonic_sort_transposed(normalized_IS_t *G) {
    FQ_ELEM* ptr[K];
    uint32_t P[K];
    for (uint32_t i = 0; i < K; ++i) {
        ptr[i] = G->values[i];
    	P[i] = i;
    }

    uint64_t top, r, i;
    uint64_t n = K;

    top = 1ul << (32 - __builtin_clz(K/2));

    for (uint64_t p = top; p > 0; p >>= 1) {
        for (i = 0; i < n - p; ++i) {
            if (!(i & p)) {
                // NOTE: here is a sign cast, this is needed, so the
            	// sign extension needed for the mask is a unsigned one.
			    const uint32_t cmp = compare_rows_bitonic_sort(ptr, i, i + p);

                const uintptr_t mask = -(1ull - (cmp >> 31));
                cswap((uintptr_t *)(&ptr[i]), (uintptr_t *)(&ptr[i+p]), mask);
				MASKED_SWAP(P[i], P[i+p], mask);
            }
        }

        i = 0;
        for (uint64_t q = top; q > p; q >>= 1) {
            for (; i < n - q; ++i) {
                if (!(i & p)) {
                    for (r = q; r > p; r >>= 1) {
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

/// In-place quicksort
/// the same as `col_lex_quicksort` except we are tracking the permutations
void canonical_col_lex_quicksort_transpose(normalized_IS_t *V) {
    normalized_IS_t VT;
    matrix_transpose_opt((uint8_t *)VT.values, (uint8_t *)V->values, K);
	col_bitonic_sort_transposed(&VT);
    matrix_transpose_opt((uint8_t *)V->values, (uint8_t *)VT.values, K);
}
