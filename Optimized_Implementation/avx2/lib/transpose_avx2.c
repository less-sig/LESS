/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (May 2025)
 *
 * @author Floyd Zweydinger <zweydfg8+github@rub.de>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

static const uint32_t matrix_transpose_table[] __attribute__((aligned(32))) = {
    0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15
};

/// \param dst_origin[out]: output matrix
/// \param src_origin[in]: input matrix
/// \param prf_origin[in]: lookahead pointer to prefetch it
/// \param src_stride[in]:
/// \param dst_stride[in]:
void matrix_transpose_32x32(uint8_t* dst_origin,
                            const uint8_t* src_origin,
                            const uint8_t* prf_origin,
                            const size_t src_stride,
                            const size_t dst_stride) {
    (void)prf_origin;
    __m256i t[32];
    for (uint32_t i = 0; i < 32; i++) {
        t[i] = _mm256_loadu_si256((const __m256i *)(src_origin + i*src_stride));
    }

    #pragma unroll
    for (uint32_t i = 0; i < 32; i+=2) {
        const __m256i t0 = _mm256_unpacklo_epi8(t[i+0], t[i+1]);
        const __m256i t1 = _mm256_unpackhi_epi8(t[i+0], t[i+1]);
        t[i+0] = t0;
        t[i+1] = t1;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 32; i+=4) {
        const __m256i t0 = _mm256_unpacklo_epi16(t[i+0], t[i+2]);
        const __m256i t1 = _mm256_unpacklo_epi16(t[i+1], t[i+3]);
        const __m256i t2 = _mm256_unpackhi_epi16(t[i+0], t[i+2]);
        const __m256i t3 = _mm256_unpackhi_epi16(t[i+1], t[i+3]);
        t[i+0] = t0;
        t[i+1] = t1;
        t[i+2] = t2;
        t[i+3] = t3;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 32; i+=8) {
        const __m256i t0 = _mm256_unpacklo_epi32(t[i+0], t[i+4]);
        const __m256i t1 = _mm256_unpacklo_epi32(t[i+1], t[i+5]);
        const __m256i t2 = _mm256_unpacklo_epi32(t[i+2], t[i+6]);
        const __m256i t3 = _mm256_unpacklo_epi32(t[i+3], t[i+7]);
        const __m256i t4 = _mm256_unpackhi_epi32(t[i+0], t[i+4]);
        const __m256i t5 = _mm256_unpackhi_epi32(t[i+1], t[i+5]);
        const __m256i t6 = _mm256_unpackhi_epi32(t[i+2], t[i+6]);
        const __m256i t7 = _mm256_unpackhi_epi32(t[i+3], t[i+7]);
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
        const __m256i t0 = _mm256_unpacklo_epi64(t[i+ 0], t[i+ 8]);
        const __m256i t1 = _mm256_unpackhi_epi64(t[i+ 0], t[i+ 8]);
        const __m256i t2 = _mm256_unpacklo_epi64(t[i+16], t[i+24]);
        const __m256i t3 = _mm256_unpackhi_epi64(t[i+16], t[i+24]);
        t[i+ 0] = t0;
        t[i+ 8] = t1;
        t[i+16] = t2;
        t[i+24] = t3;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 16; i++) {
        const __m256i t0 = _mm256_permute2x128_si256(t[i+0], t[i+16], 0b100000);
        const __m256i t1 = _mm256_permute2x128_si256(t[i+0], t[i+16], 0b110001);
        t[i+ 0] = t0;
        t[i+16] = t1;
    }

    #pragma unroll
    for (uint32_t i = 0; i < 16; i++) {
        const uint32_t pos = matrix_transpose_table[i];
        _mm256_storeu_si256((__m256i *)(dst_origin + i*dst_stride), t[pos]);
    }

    for (uint32_t i = 0; i < 16; i++) {
        const uint32_t pos = matrix_transpose_table[i];
        _mm256_storeu_si256((__m256i *)(dst_origin + (i+16)*dst_stride), t[pos+16]);
    }
}


// TODO adapt
static void gf16_transpose_64x64_avx2(uint8_t *B,
                                      const uint8_t *const A,
                                      const uint32_t src_stride,
                                      const uint32_t dst_stride) {
    __m256i M[64];
    const __m256i mask1 = _mm256_set1_epi8  (0x0F),
                  mask2 = _mm256_set1_epi16 (0x00FF),
                  mask3 = _mm256_set1_epi32 (0x0000FFFF),
                  mask4 = _mm256_set1_epi64x(0x00000000FFFFFFFF),
                  mask5 = _mm256_setr_epi64x(0xFFFFFFFFFFFFFFFF, 0, 0xFFFFFFFFFFFFFFFF, 0),
                  mask6 = _mm256_setr_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0, 0);
    for (uint32_t i = 0; i < 64; i++) {
        M[i] = _mm256_loadu_si256((__m256i *)(A+i*src_stride));
    }

    for (uint32_t i = 0; i < 64; i+=2) {
        const __m256i t = (_mm256_srli_epi64(M[i], 4) ^ M[i+1]) & mask1;
        M[i+0] ^= _mm256_slli_epi64(t, 4);
        M[i+1] ^= t;
    }

    for (uint32_t i = 0; i < 64; i+=4) {
        const __m256i t0 = (_mm256_srli_epi64(M[i+0], 8) ^ M[i+2]) & mask2;
        const __m256i t1 = (_mm256_srli_epi64(M[i+1], 8) ^ M[i+3]) & mask2;
        M[i+0] ^= _mm256_slli_epi64(t0, 8);
        M[i+1] ^= _mm256_slli_epi64(t1, 8);
        M[i+2] ^= t0;
        M[i+3] ^= t1;
    }

    for (uint32_t i = 0; i < 64; i+=8) {
        const __m256i t0 = (_mm256_srli_epi64(M[i+0], 16) ^ M[i+4]) & mask3;
        const __m256i t1 = (_mm256_srli_epi64(M[i+1], 16) ^ M[i+5]) & mask3;
        const __m256i t2 = (_mm256_srli_epi64(M[i+2], 16) ^ M[i+6]) & mask3;
        const __m256i t3 = (_mm256_srli_epi64(M[i+3], 16) ^ M[i+7]) & mask3;
        M[i+0] ^= _mm256_slli_epi64(t0, 16);
        M[i+1] ^= _mm256_slli_epi64(t1, 16);
        M[i+2] ^= _mm256_slli_epi64(t2, 16);
        M[i+3] ^= _mm256_slli_epi64(t3, 16);
        M[i+4] ^= t0;
        M[i+5] ^= t1;
        M[i+6] ^= t2;
        M[i+7] ^= t3;
    }

    for (uint32_t i = 0; i < 64; i+=16) {
        const __m256i t0 = (_mm256_srli_epi64(M[i+0], 32) ^ M[i+ 8]) & mask4;
        const __m256i t1 = (_mm256_srli_epi64(M[i+1], 32) ^ M[i+ 9]) & mask4;
        const __m256i t2 = (_mm256_srli_epi64(M[i+2], 32) ^ M[i+10]) & mask4;
        const __m256i t3 = (_mm256_srli_epi64(M[i+3], 32) ^ M[i+11]) & mask4;
        const __m256i t4 = (_mm256_srli_epi64(M[i+4], 32) ^ M[i+12]) & mask4;
        const __m256i t5 = (_mm256_srli_epi64(M[i+5], 32) ^ M[i+13]) & mask4;
        const __m256i t6 = (_mm256_srli_epi64(M[i+6], 32) ^ M[i+14]) & mask4;
        const __m256i t7 = (_mm256_srli_epi64(M[i+7], 32) ^ M[i+15]) & mask4;
        M[i+ 0] ^= _mm256_slli_epi64(t0, 32);
        M[i+ 1] ^= _mm256_slli_epi64(t1, 32);
        M[i+ 2] ^= _mm256_slli_epi64(t2, 32);
        M[i+ 3] ^= _mm256_slli_epi64(t3, 32);
        M[i+ 4] ^= _mm256_slli_epi64(t4, 32);
        M[i+ 5] ^= _mm256_slli_epi64(t5, 32);
        M[i+ 6] ^= _mm256_slli_epi64(t6, 32);
        M[i+ 7] ^= _mm256_slli_epi64(t7, 32);
        M[i+ 8] ^= t0;
        M[i+ 9] ^= t1;
        M[i+10] ^= t2;
        M[i+11] ^= t3;
        M[i+12] ^= t4;
        M[i+13] ^= t5;
        M[i+14] ^= t6;
        M[i+15] ^= t7;
    }

    for (uint32_t i = 0; i < 16; i++) {
        const __m256i t0 = (_mm256_srli_si256(M[i+ 0], 8) ^ M[i+16]) & mask5;
        const __m256i t1 = (_mm256_srli_si256(M[i+32], 8) ^ M[i+48]) & mask5;
        M[i   ] ^= _mm256_slli_si256(t0, 8);
        M[i+32] ^= _mm256_slli_si256(t1, 8);
        M[i+16] ^= t0;
        M[i+48] ^= t1;
    }

    for (uint32_t i = 0; i < 32; i++) {
        const __m256i t = (_mm256_permute2x128_si256(M[i+0], M[i+0], 0b10000001) ^ M[i+32]) & mask6;
        M[i   ] ^= _mm256_permute2x128_si256(t, t, 0b01000); //
        M[i+32] ^= t;
    }

    // write out
    for (uint32_t i = 0; i < 64; i++) {
        _mm256_storeu_si256((__m256i *)(B + i*dst_stride), M[i]);
    }
}

