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
/// \param src_stride[in]: number of bytes between two input rows
/// \param dst_stride[in]: number of bytes between two input cols
void matrix_transpose_32x32(uint8_t* dst_origin,
                            const uint8_t* src_origin,
                            const size_t src_stride,
                            const size_t dst_stride) {
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

    #pragma unroll
    for (uint32_t i = 0; i < 16; i++) {
        const uint32_t pos = matrix_transpose_table[i];
        _mm256_storeu_si256((__m256i *)(dst_origin + (i+16)*dst_stride), t[pos+16]);
    }
}
