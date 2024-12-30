
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

const static int8_t sortingnetwork_u8x32_shuffle_masks[8][32] __attribute((aligned(64))) = {
        {1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14,1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14},
        {3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12,3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12},
        {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8},
        {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
        {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0},
        {4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11},
	// NOTE: this is a special mask, currently only needed in `sortingnetwork_aftermergesort_u8x64`
        {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
	// TODO: this seems to be a bug in the implementation. Currently only needed for `sortingnetwork_mergesort_u8x64`
        {0,2,1,3,4,6,5,7, 8,10,9,11, 12,14,13,15, 0,2,1,3,4,6,5,7, 8,10,9,11, 12,14,13,15},
};


/// \param  a = [a0, ..., a31], u8 elements
/// \return a = [a31, ..., a0]
static inline void sortingnetwork_reverse_u8x32(__m256i *a) {
    const __m256i p = _mm256_permute2x128_si256(*a, *a, 0b00000001);
	const __m256i mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[4]);
	*a =_mm256_shuffle_epi8(p, mask);
}

/// TODO doc
__m256i sortingnetwork_sort_u8x32(__m256i v) {
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

    __m128i tmp_;
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

///
static inline void sortingnetwork_aftermergesort_u8x32(__m256i *a)  {
	__m256i mask, tmp, L0;
    L0 = (__m256i)_mm256_permute_ps((__m256)*a, 0b10110001);

	// 3
    COEX_u8x32(*a, L0, tmp);
	mask = _mm256_set1_epi64x(0xFFFFFFFF);
	L0 = _mm256_blendv_epi8(L0, *a, mask);

	// 4
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]);
	__m256i L4p = _mm256_shuffle_epi8(L0, mask);
    COEX_u8x32(L0, L4p, tmp);
	mask = _mm256_set1_epi16(0x00FF);
	L0 = _mm256_blendv_epi8(L4p, L0, mask);

	// 5
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[6]);
	__m256i L5p = _mm256_shuffle_epi8(L0, mask);
    COEX_u8x32(L0, L5p, tmp);
	mask = _mm256_set1_epi32(0x0000FFFF);
	L0 = _mm256_blendv_epi8(L5p, L0, mask);

	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[7]);
	__m256i L6p = _mm256_shuffle_epi8(L0, mask);
    COEX_u8x32(L0, L6p, tmp);
	mask = _mm256_set1_epi32(0x0000FF00);
	*a = _mm256_blendv_epi8(L6p, L0, mask);
}

// implementation of `simd_aftermerge_4V`
static inline void sortingnetwork_aftermerge_u8x32(__m256i *a) {
	__m256i tmp, L0 = *a;

    __m256i L1p = _mm256_permute2x128_si256(L0, L0, 0b00000001);
    COEX_u8x32(L0, L1p, tmp);
	L0 = _mm256_blend_epi32(L0, L1p, 0b11110000);

	// 2 (a b, a c)
	__m256i L2p = _mm256_permute4x64_epi64(L0, 0b10110001);
    COEX_u8x32(L0, L2p, tmp);
	*a = _mm256_blend_epi32(L0, L2p, 0b11001100);

    sortingnetwork_aftermergesort_u8x32(a);
}

/// implementation of 8 parallel "simd_aftermerge_1V",  NOT the implementation of 2 parallel `simd_aftermerge_2V`
static inline void sortingnetwork_aftermergesort_u8x64(__m256i *a,
													   __m256i *b) {
	__m256i mask, tmp;
    __m256i L0 = (__m256i)_mm256_permute_ps((__m256)*a, 0b10110001);
    __m256i H0 = (__m256i)_mm256_permute_ps((__m256)*b, 0b10110001);

	// 3
    COEX_u8x32(*a, L0, tmp);
    COEX_u8x32(*b, H0, tmp);
	mask = _mm256_set1_epi64x(0xFFFFFFFF);
	L0 = _mm256_blendv_epi8(L0, *a, mask);
	H0 = _mm256_blendv_epi8(H0, *b, mask);

	// 4
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[0]);
	__m256i L4p = _mm256_shuffle_epi8(L0, mask);
	__m256i H4p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L4p, tmp);
    COEX_u8x32(H0, H4p, tmp);
	mask = _mm256_set1_epi16(0x00FF);
	L0 = _mm256_blendv_epi8(L4p, L0, mask);
	H0 = _mm256_blendv_epi8(H4p, H0, mask);

	// 5
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[6]);
	__m256i L5p = _mm256_shuffle_epi8(L0, mask);
	__m256i H5p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L5p, tmp);
    COEX_u8x32(H0, H5p, tmp);
	mask = _mm256_set1_epi32(0x0000FFFF);
	L0 = _mm256_blendv_epi8(L5p, L0, mask);
	H0 = _mm256_blendv_epi8(H5p, H0, mask);

	// 6, TODO somewhere is a bug, for reasons I dont understand this shuffle is missing
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[7]);
	__m256i L6p = _mm256_shuffle_epi8(L0, mask);
	__m256i H6p = _mm256_shuffle_epi8(H0, mask);
    COEX_u8x32(L0, L6p, tmp);
    COEX_u8x32(H0, H6p, tmp);
	mask = _mm256_set1_epi32(0x0000FF00);
	*a = _mm256_blendv_epi8(L6p, L0, mask);
	*b = _mm256_blendv_epi8(H6p, H0, mask);
}

// implementation of `simd_sort_4V_sorted`
static inline void sortingnetwork_mergesort_u8x64(__m256i *a,
                                                  __m256i *b) {

    __m256i L0 = *a, tmp, mask;
    __m256i H0 = *b;

	// reverse H0
    __m256i H0p = _mm256_permute2x128_si256(H0, H0, 0b00000001);
	mask = _mm256_load_si256((__m256i *)sortingnetwork_u8x32_shuffle_masks[4]);
	H0p =_mm256_shuffle_epi8(H0p, mask);
	H0 = H0p;

    // 0
    COEX_u8x32(L0, H0, tmp);

	// 1 (a c, b d)
    __m256i L1p = _mm256_permute2x128_si256(L0, L0, 0b00000001);
    __m256i H1p = _mm256_permute2x128_si256(H0, H0, 0b00000001);
    COEX_u8x32(L0, L1p, tmp);
    COEX_u8x32(H0, H1p, tmp);
    // __m256i L1p2= _mm256_permute2x128_si256(L1p, L1p, 0b00000000); // TODO optimize
    // __m256i H1p2= _mm256_permute2x128_si256(H1p, H1p, 0b00000000);
	L0 = _mm256_blend_epi32(L0, L1p, 0b11110000);
	H0 = _mm256_blend_epi32(H0, H1p, 0b11110000);

	// 2 (a b, a c)
	__m256i L2p = _mm256_permute4x64_epi64(L0, 0b10110001);
	__m256i H2p = _mm256_permute4x64_epi64(H0, 0b10110001);
    COEX_u8x32(L0, L2p, tmp);
    COEX_u8x32(H0, H2p, tmp);
	// __m256i L2p2= _mm256_permute4x64_epi64(L2p, 0b10110001);
	// __m256i H2p2= _mm256_permute4x64_epi64(H2p, 0b10110001);
	*a = _mm256_blend_epi32(L0, L2p, 0b11001100);
	*b = _mm256_blend_epi32(H0, H2p, 0b11001100);

	sortingnetwork_aftermergesort_u8x64(a, b);
}

// implementation of "simd_sort_8V"
static inline void sortingnetwork_sort_u8x64(__m256i *a,
						                     __m256i *b) {
    *a = sortingnetwork_sort_u8x32(*a);
    *b = sortingnetwork_sort_u8x32(*b);
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

/// implementation of `simd_sort_16V`
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


/// implementation of `simd_aftermerge_16V`
static inline void sortingnetwork_aftermerge_u8x128(__m256i *a,
												    __m256i *b,
                                                    __m256i *c,
                                                    __m256i *d) {
	__m256i tmp;
    COEX_u8x32(*a, *c, tmp);
    COEX_u8x32(*b, *d, tmp);
    sortingnetwork_aftermerge_u8x64(a, b);
    sortingnetwork_aftermerge_u8x64(c, d);
}

///
static inline void sortingnetwork_sort_u8x256(__m256i *a,
											  __m256i *b,
                                              __m256i *c,
                                              __m256i *d,
                                              __m256i *e,
                                              __m256i *f,
                                              __m256i *g,
                                              __m256i *h) {
    __m256i tmp;
	sortingnetwork_sort_u8x128(a, b, c, d);
	sortingnetwork_sort_u8x128(e, f, g, h);

    sortingnetwork_reverse_u8x32(e);
    sortingnetwork_reverse_u8x32(f);
    sortingnetwork_reverse_u8x32(g);
    sortingnetwork_reverse_u8x32(h);

    COEX_u8x32(*d, *e, tmp);
    COEX_u8x32(*c, *f, tmp);
    COEX_u8x32(*b, *g, tmp);
    COEX_u8x32(*a, *h, tmp);

    sortingnetwork_aftermerge_u8x128(a, b, c, d);
    sortingnetwork_aftermerge_u8x128(e, f, g, h);
}

///
static inline void sortingnetwork_aftermerge_u8x256(__m256i *a,
												    __m256i *b,
                                                    __m256i *c,
                                                    __m256i *d,
                                                    __m256i *e,
                                                    __m256i *f,
                                                    __m256i *g,
                                                    __m256i *h) {
	__m256i tmp;
    COEX_u8x32(*a, *e, tmp);
    COEX_u8x32(*b, *f, tmp);
    COEX_u8x32(*c, *g, tmp);
    COEX_u8x32(*d, *h, tmp);
    sortingnetwork_aftermerge_u8x128(a, b, c, d);
    sortingnetwork_aftermerge_u8x128(e, f, g, h);
}

/// implementation of `simd_sort_33V`
static inline void sortingnetwork_sort_u8x288(__m256i *a,
											  __m256i *b,
											  __m256i *c,
											  __m256i *d,
											  __m256i *e,
											  __m256i *f,
											  __m256i *g,
											  __m256i *h,
											  __m256i *i) {
	__m256i tmp;
	sortingnetwork_sort_u8x256(a, b, c, d, e, f, g, h);
	*i = sortingnetwork_sort_u8x32(*i);
	sortingnetwork_reverse_u8x32(i);
	COEX_u8x32(*h, *i, tmp);
	sortingnetwork_aftermerge_u8x256(a, b, c, d, e, f, g, h);
	sortingnetwork_aftermerge_u8x32(i);
}


///
static inline void sortingnetwork_sort_u8x512(__m256i *a,
                                              __m256i *b,
                                              __m256i *c,
                                              __m256i *d,
                                              __m256i *e,
                                              __m256i *f,
                                              __m256i *g,
                                              __m256i *h,
                                              __m256i *i,
											  __m256i *j,
                                              __m256i *k,
                                              __m256i *l,
                                              __m256i *m,
                                              __m256i *n,
                                              __m256i *o,
                                              __m256i *p) {
    __m256i tmp;
	sortingnetwork_sort_u8x256(a, b, c, d, e, f, g, h);
	sortingnetwork_sort_u8x256(i, j, k, l, m, n, o, p);

    sortingnetwork_reverse_u8x32(i);
    sortingnetwork_reverse_u8x32(j);
    sortingnetwork_reverse_u8x32(k);
    sortingnetwork_reverse_u8x32(l);
    sortingnetwork_reverse_u8x32(m);
    sortingnetwork_reverse_u8x32(n);
    sortingnetwork_reverse_u8x32(o);
    sortingnetwork_reverse_u8x32(p);

    COEX_u8x32(*a, *p, tmp);
    COEX_u8x32(*b, *o, tmp);
    COEX_u8x32(*c, *n, tmp);
    COEX_u8x32(*d, *m, tmp);
    COEX_u8x32(*e, *l, tmp);
    COEX_u8x32(*f, *k, tmp);
    COEX_u8x32(*g, *j, tmp);
    COEX_u8x32(*h, *i, tmp);

    sortingnetwork_aftermerge_u8x256(a, b, c, d, e, f, g, h);
    sortingnetwork_aftermerge_u8x256(i, j, k, l, m, n, o, p);
}


void sortingnetwork(FQ_ELEM *arr,
                    const size_t s) {
	uint8_t tmp[32] __attribute__((aligned(32)));
#ifdef CATEGORY_1
    __m256i i1 = _mm256_loadu_si256((const __m256i *)(arr +  0));
    __m256i i2 = _mm256_loadu_si256((const __m256i *)(arr + 32));
    __m256i i3 = _mm256_loadu_si256((const __m256i *)(arr + 64));

	uint32_t i = 0;
	for (; i < (s%32); i++) { tmp[i] = arr[96 + i];}
	for (; i < 32; i++) { tmp[i] = -1;}
    __m256i i4 = _mm256_loadu_si256((const __m256i *)tmp);
    sortingnetwork_sort_u8x128(&i1, &i2, &i3, &i4);
    _mm256_storeu_si256((__m256i *)(arr +  0), i1);
    _mm256_storeu_si256((__m256i *)(arr + 32), i2);
    _mm256_storeu_si256((__m256i *)(arr + 64), i3);
    _mm256_storeu_si256((__m256i *)tmp, i4);
	for (uint32_t i = 0; i < (s%32); i++) { arr[96 + i] =  tmp[i];}
#elif defined(CATEGORY_3)
    __m256i i1 = _mm256_loadu_si256((const __m256i *)(arr +  0));
    __m256i i2 = _mm256_loadu_si256((const __m256i *)(arr + 32));
    __m256i i3 = _mm256_loadu_si256((const __m256i *)(arr + 64));
    __m256i i4 = _mm256_loadu_si256((const __m256i *)(arr + 96));
    __m256i i5 = _mm256_loadu_si256((const __m256i *)(arr +128));
    __m256i i6 = _mm256_loadu_si256((const __m256i *)(arr +160));
	__m256i i8 = _mm256_set1_epi8(-1);

	uint32_t i = 0;
	for (; i < (s%32); i++) { tmp[i] = arr[192 + i];}
	for (; i < 32; i++) { tmp[i] = -1;}
    __m256i i7 = _mm256_loadu_si256((const __m256i *)tmp);
    sortingnetwork_sort_u8x256(&i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8);
    _mm256_storeu_si256((__m256i *)(arr +  0), i1);
    _mm256_storeu_si256((__m256i *)(arr + 32), i2);
    _mm256_storeu_si256((__m256i *)(arr + 64), i3);
    _mm256_storeu_si256((__m256i *)(arr + 96), i4);
    _mm256_storeu_si256((__m256i *)(arr +128), i5);
    _mm256_storeu_si256((__m256i *)(arr +160), i6);
    _mm256_storeu_si256((__m256i *)tmp, i7);
	for (uint32_t i = 0; i < (s%32); i++) { arr[192 + i] =  tmp[i];}
#elif defined(CATEGORY_5)
    __m256i  i1 = _mm256_loadu_si256((const __m256i *)(arr +  0));
    __m256i  i2 = _mm256_loadu_si256((const __m256i *)(arr + 32));
    __m256i  i3 = _mm256_loadu_si256((const __m256i *)(arr + 64));
    __m256i  i4 = _mm256_loadu_si256((const __m256i *)(arr + 96));
    __m256i  i5 = _mm256_loadu_si256((const __m256i *)(arr +128));
    __m256i  i6 = _mm256_loadu_si256((const __m256i *)(arr +160));
    __m256i  i7 = _mm256_loadu_si256((const __m256i *)(arr +192));
    __m256i  i8 = _mm256_loadu_si256((const __m256i *)(arr +224));

	uint32_t i = 0;
	for (; i < (s%32); i++) { tmp[i] = arr[256 + i];}
	for (; i < 32; i++) { tmp[i] = -1;}
    __m256i i9 = _mm256_loadu_si256((const __m256i *)tmp);

    sortingnetwork_sort_u8x288(&i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9);

    _mm256_storeu_si256((__m256i *)(arr +  0), i1);
    _mm256_storeu_si256((__m256i *)(arr + 32), i2);
    _mm256_storeu_si256((__m256i *)(arr + 64), i3);
    _mm256_storeu_si256((__m256i *)(arr + 96), i4);
    _mm256_storeu_si256((__m256i *)(arr +128), i5);
    _mm256_storeu_si256((__m256i *)(arr +160), i6);
    _mm256_storeu_si256((__m256i *)(arr +192), i7);
    _mm256_storeu_si256((__m256i *)(arr +224), i8);
    _mm256_storeu_si256((__m256i *)tmp, i9);
	for (uint32_t k = 0; k < (s%32); k++) { arr[256 + k] =  tmp[k];}
#endif
}
