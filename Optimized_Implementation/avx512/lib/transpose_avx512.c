#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

static const uint32_t matrix_transpose_table[] __attribute__((aligned(32))) = {
    0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15
};

/// \param dst_origin[out]: output matrix
/// \param src_origin[in]: input matrix
/// \param src_stride[in]: number of bytes between two rows
/// \param dst_stride[in]: number of bytes between two cols
void matrix_transpose_64x64(uint8_t* dst_origin,
                            const uint8_t* src_origin,
                            const size_t src_stride,
                            const size_t dst_stride) {
    const __m512i m1 = _mm512_setr_epi64(0b0000, 0b0001, 0b1000, 0b1001, 0b0100, 0b0101, 0b1100, 0b1101);
    const __m512i m2 = _mm512_setr_epi64(0b0010, 0b0011, 0b1010, 0b1011, 0b0110, 0b0111, 0b1110, 0b1111);

    const __m512i m3 = _mm512_setr_epi64(0b0000, 0b0001, 0b0010, 0b0011, 0b1000, 0b1001, 0b1010, 0b1011);
    const __m512i m4 = _mm512_setr_epi64(0b0100, 0b0101, 0b0110, 0b0111, 0b1100, 0b1101, 0b1110, 0b1111);

    __m512i t[64];
    for (uint32_t i = 0; i < 64; i++) {
        t[i] = _mm512_loadu_si512((const __m512i *)(src_origin + i*src_stride));
    }

    #pragma unroll
    for (uint32_t i = 0; i < 64; i+=2) {
        const __m512i t0 = _mm512_unpacklo_epi8(t[i+0], t[i+1]);
        const __m512i t1 = _mm512_unpackhi_epi8(t[i+0], t[i+1]);
        t[i+0] = t0;
        t[i+1] = t1;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 64; i+=4) {
        const __m512i t0 = _mm512_unpacklo_epi16(t[i+0], t[i+2]);
        const __m512i t1 = _mm512_unpacklo_epi16(t[i+1], t[i+3]);
        const __m512i t2 = _mm512_unpackhi_epi16(t[i+0], t[i+2]);
        const __m512i t3 = _mm512_unpackhi_epi16(t[i+1], t[i+3]);
        t[i+0] = t0;
        t[i+1] = t1;
        t[i+2] = t2;
        t[i+3] = t3;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 64; i+=8) {
        const __m512i t0 = _mm512_unpacklo_epi32(t[i+0], t[i+4]);
        const __m512i t1 = _mm512_unpacklo_epi32(t[i+1], t[i+5]);
        const __m512i t2 = _mm512_unpacklo_epi32(t[i+2], t[i+6]);
        const __m512i t3 = _mm512_unpacklo_epi32(t[i+3], t[i+7]);
        const __m512i t4 = _mm512_unpackhi_epi32(t[i+0], t[i+4]);
        const __m512i t5 = _mm512_unpackhi_epi32(t[i+1], t[i+5]);
        const __m512i t6 = _mm512_unpackhi_epi32(t[i+2], t[i+6]);
        const __m512i t7 = _mm512_unpackhi_epi32(t[i+3], t[i+7]);
        t[i+0] = t0;
        t[i+1] = t1;
        t[i+2] = t2;
        t[i+3] = t3;
        t[i+4] = t4;
        t[i+5] = t5;
        t[i+6] = t6;
        t[i+7] = t7;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 8; i++) {
        const __m512i t0 = _mm512_unpacklo_epi64(t[i+ 0], t[i+ 8]);
        const __m512i t1 = _mm512_unpackhi_epi64(t[i+ 0], t[i+ 8]);
        const __m512i t2 = _mm512_unpacklo_epi64(t[i+16], t[i+24]);
        const __m512i t3 = _mm512_unpackhi_epi64(t[i+16], t[i+24]);
        const __m512i t4 = _mm512_unpacklo_epi64(t[i+32], t[i+40]);
        const __m512i t5 = _mm512_unpackhi_epi64(t[i+32], t[i+40]);
        const __m512i t6 = _mm512_unpacklo_epi64(t[i+48], t[i+56]);
        const __m512i t7 = _mm512_unpackhi_epi64(t[i+48], t[i+56]);
        t[i+ 0] = t0;
        t[i+ 8] = t1;
        t[i+16] = t2;
        t[i+24] = t3;
        t[i+32] = t4;
        t[i+40] = t5;
        t[i+48] = t6;
        t[i+56] = t7;
    }

    // swap 128 bit limbs
    #pragma unroll
    for (uint32_t i = 0; i < 16; i++) {
        const __m512i t0 = _mm512_permutex2var_epi64(t[i+ 0], m1, t[i+16]);
        const __m512i t1 = _mm512_permutex2var_epi64(t[i+ 0], m2, t[i+16]);
        const __m512i t2 = _mm512_permutex2var_epi64(t[i+32], m1, t[i+48]);
        const __m512i t3 = _mm512_permutex2var_epi64(t[i+32], m2, t[i+48]);
        t[i+ 0] = t0;
        t[i+16] = t1;
        t[i+32] = t2;
        t[i+48] = t3;
    }

    // swap 256 bit limbs
    #pragma unroll
    for (uint32_t i = 0; i < 32; i++) {
        const __m512i t0 = _mm512_permutex2var_epi64(t[i+0], m3, t[i+32]);
        const __m512i t1 = _mm512_permutex2var_epi64(t[i+0], m4, t[i+32]);
        t[i+ 0] = t0;
        t[i+32] = t1;
    }


    #pragma unroll
    for (uint32_t j = 0; j < 4; j++) {
        const uint32_t off = j*16;

        #pragma unroll
        for (uint32_t i = 0; i < 16; i++) {
            const uint32_t pos = matrix_transpose_table[i];
            _mm512_storeu_si512((__m512i *)(dst_origin + (off+i)*dst_stride), t[off+pos]);
        }
    }
}
