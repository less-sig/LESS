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

/// needed for matrix_transpose_32x32_avx2
static const uint8_t BLENDV_MASK[5][32] __attribute__((aligned(32)))= {
    { 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff },
    { 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0xff, 0xff },
    { 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff },
    { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff },
    { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff },
};

// needed for matrix_transpose_32x32_avx2
static const uint8_t SHUFFLE_MASK[4][32] __attribute__((aligned(32))) = {
    { 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30 },
    { 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13, 18, 19, 16, 17, 22, 23, 20, 21, 26, 27, 24, 25, 30, 31, 28, 29 },
    { 4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 20, 21, 22, 23, 16, 17, 18, 19, 28, 29, 30, 31, 24, 25, 26, 27 },
    { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22, 23 },
};

typedef __m256i_u LOAD_TYPE;
typedef __m256i_u STORE_TYPE;

/// \param dst_origin[out]: output matrix
/// \param src_origin[in]: input matrix
/// \param prf_origin[in]: lookahead pointer to prefetch it
/// \param src_stride[in]:
/// \param dst_stride[in]:
void matrix_transpose_32x32_(uint8_t* dst_origin,
                                 const uint8_t* src_origin,
                                 const uint8_t* prf_origin,
                                 const size_t src_stride,
                                 const size_t dst_stride) {
    const __m256i shm_1 = _mm256_load_si256((const __m256i *)SHUFFLE_MASK[0]);
    const __m256i blm_1 = _mm256_load_si256((const __m256i *)BLENDV_MASK[0]);
    __m256i rnd_0_0 = *(const LOAD_TYPE*)(src_origin + 0*src_stride);
    __m256i rnd_0_1 = *(const LOAD_TYPE*)(src_origin + 1*src_stride);
    __m256i shf_1_0 = _mm256_shuffle_epi8(rnd_0_0, shm_1);
    __m256i shf_1_1 = _mm256_shuffle_epi8(rnd_0_1, shm_1);
    __m256i rnd_1_0 = _mm256_blendv_epi8(rnd_0_0, shf_1_1, blm_1);
    __m256i rnd_1_1 = _mm256_blendv_epi8(shf_1_0, rnd_0_1, blm_1);
    __m256i rnd_1_0_= _mm256_unpacklo_epi8(rnd_0_0, rnd_0_1);
    __m256i rnd_0_2 = *(const LOAD_TYPE*)(src_origin + 2*src_stride);
    __m256i rnd_0_3 = *(const LOAD_TYPE*)(src_origin + 3*src_stride);
    __m256i shf_1_2 = _mm256_shuffle_epi8(rnd_0_2, shm_1);
    __m256i shf_1_3 = _mm256_shuffle_epi8(rnd_0_3, shm_1);
    __m256i rnd_1_2 = _mm256_blendv_epi8(rnd_0_2, shf_1_3, blm_1);
    __m256i rnd_1_3 = _mm256_blendv_epi8(shf_1_2, rnd_0_3, blm_1);
    __m256i rnd_0_4 = *(const LOAD_TYPE*)(src_origin + 4*src_stride);
    __m256i rnd_0_5 = *(const LOAD_TYPE*)(src_origin + 5*src_stride);
    __m256i shf_1_4 = _mm256_shuffle_epi8(rnd_0_4, shm_1);
    __m256i shf_1_5 = _mm256_shuffle_epi8(rnd_0_5, shm_1);
    __m256i rnd_1_4 = _mm256_blendv_epi8(rnd_0_4, shf_1_5, blm_1);
    __m256i rnd_1_5 = _mm256_blendv_epi8(shf_1_4, rnd_0_5, blm_1);
    __m256i rnd_0_6 = *(const LOAD_TYPE*)(src_origin + 6*src_stride);
    __m256i rnd_0_7 = *(const LOAD_TYPE*)(src_origin + 7*src_stride);
    __m256i shf_1_6 = _mm256_shuffle_epi8(rnd_0_6, shm_1);
    __m256i shf_1_7 = _mm256_shuffle_epi8(rnd_0_7, shm_1);
    _mm_prefetch(prf_origin+0*src_stride, _MM_HINT_NTA);
    __m256i rnd_1_6 = _mm256_blendv_epi8(rnd_0_6, shf_1_7, blm_1);
    __m256i rnd_1_7 = _mm256_blendv_epi8(shf_1_6, rnd_0_7, blm_1);
    __m256i rnd_0_8 = *(const LOAD_TYPE*)(src_origin + 8*src_stride);
    __m256i rnd_0_9 = *(const LOAD_TYPE*)(src_origin + 9*src_stride);
    __m256i shf_1_8 = _mm256_shuffle_epi8(rnd_0_8, shm_1);
    __m256i shf_1_9 = _mm256_shuffle_epi8(rnd_0_9, shm_1);
    __m256i rnd_1_8 = _mm256_blendv_epi8(rnd_0_8, shf_1_9, blm_1);
    __m256i rnd_1_9 = _mm256_blendv_epi8(shf_1_8, rnd_0_9, blm_1);
    __m256i rnd_0_10 = *(const LOAD_TYPE*)(src_origin + 10*src_stride);
    __m256i rnd_0_11 = *(const LOAD_TYPE*)(src_origin + 11*src_stride);
    __m256i shf_1_10 = _mm256_shuffle_epi8(rnd_0_10, shm_1);
    __m256i shf_1_11 = _mm256_shuffle_epi8(rnd_0_11, shm_1);
    __m256i rnd_1_10 = _mm256_blendv_epi8(rnd_0_10, shf_1_11, blm_1);
    __m256i rnd_1_11 = _mm256_blendv_epi8(shf_1_10, rnd_0_11, blm_1);
    __m256i rnd_0_12 = *(const LOAD_TYPE*)(src_origin + 12*src_stride);
    __m256i rnd_0_13 = *(const LOAD_TYPE*)(src_origin + 13*src_stride);
    __m256i shf_1_12 = _mm256_shuffle_epi8(rnd_0_12, shm_1);
    __m256i shf_1_13 = _mm256_shuffle_epi8(rnd_0_13, shm_1);
    __m256i rnd_1_12 = _mm256_blendv_epi8(rnd_0_12, shf_1_13, blm_1);
    __m256i rnd_1_13 = _mm256_blendv_epi8(shf_1_12, rnd_0_13, blm_1);
    __m256i rnd_0_14 = *(const LOAD_TYPE*)(src_origin + 14*src_stride);
    __m256i rnd_0_15 = *(const LOAD_TYPE*)(src_origin + 15*src_stride);
    __m256i shf_1_14 = _mm256_shuffle_epi8(rnd_0_14, shm_1);
    __m256i shf_1_15 = _mm256_shuffle_epi8(rnd_0_15, shm_1);
    _mm_prefetch(prf_origin+1*src_stride, _MM_HINT_NTA);
    __m256i rnd_1_14 = _mm256_blendv_epi8(rnd_0_14, shf_1_15, blm_1);
    __m256i rnd_1_15 = _mm256_blendv_epi8(shf_1_14, rnd_0_15, blm_1);
    __m256i rnd_0_16 = *(const LOAD_TYPE*)(src_origin + 16*src_stride);
    __m256i rnd_0_17 = *(const LOAD_TYPE*)(src_origin + 17*src_stride);
    __m256i shf_1_16 = _mm256_shuffle_epi8(rnd_0_16, shm_1);
    __m256i shf_1_17 = _mm256_shuffle_epi8(rnd_0_17, shm_1);
    __m256i rnd_1_16 = _mm256_blendv_epi8(rnd_0_16, shf_1_17, blm_1);
    __m256i rnd_1_17 = _mm256_blendv_epi8(shf_1_16, rnd_0_17, blm_1);
    __m256i rnd_0_18 = *(const LOAD_TYPE*)(src_origin + 18*src_stride);
    __m256i rnd_0_19 = *(const LOAD_TYPE*)(src_origin + 19*src_stride);
    __m256i shf_1_18 = _mm256_shuffle_epi8(rnd_0_18, shm_1);
    __m256i shf_1_19 = _mm256_shuffle_epi8(rnd_0_19, shm_1);
    __m256i rnd_1_18 = _mm256_blendv_epi8(rnd_0_18, shf_1_19, blm_1);
    __m256i rnd_1_19 = _mm256_blendv_epi8(shf_1_18, rnd_0_19, blm_1);
    __m256i rnd_0_20 = *(const LOAD_TYPE*)(src_origin + 20*src_stride);
    __m256i rnd_0_21 = *(const LOAD_TYPE*)(src_origin + 21*src_stride);
    __m256i shf_1_20 = _mm256_shuffle_epi8(rnd_0_20, shm_1);
    __m256i shf_1_21 = _mm256_shuffle_epi8(rnd_0_21, shm_1);
    __m256i rnd_1_20 = _mm256_blendv_epi8(rnd_0_20, shf_1_21, blm_1);
    __m256i rnd_1_21 = _mm256_blendv_epi8(shf_1_20, rnd_0_21, blm_1);
    __m256i rnd_0_22 = *(const LOAD_TYPE*)(src_origin + 22*src_stride);
    __m256i rnd_0_23 = *(const LOAD_TYPE*)(src_origin + 23*src_stride);
    __m256i shf_1_22 = _mm256_shuffle_epi8(rnd_0_22, shm_1);
    __m256i shf_1_23 = _mm256_shuffle_epi8(rnd_0_23, shm_1);
    _mm_prefetch(prf_origin+2*src_stride, _MM_HINT_NTA);
    __m256i rnd_1_22 = _mm256_blendv_epi8(rnd_0_22, shf_1_23, blm_1);
    __m256i rnd_1_23 = _mm256_blendv_epi8(shf_1_22, rnd_0_23, blm_1);
    __m256i rnd_0_24 = *(const LOAD_TYPE*)(src_origin + 24*src_stride);
    __m256i rnd_0_25 = *(const LOAD_TYPE*)(src_origin + 25*src_stride);
    __m256i shf_1_24 = _mm256_shuffle_epi8(rnd_0_24, shm_1);
    __m256i shf_1_25 = _mm256_shuffle_epi8(rnd_0_25, shm_1);
    __m256i rnd_1_24 = _mm256_blendv_epi8(rnd_0_24, shf_1_25, blm_1);
    __m256i rnd_1_25 = _mm256_blendv_epi8(shf_1_24, rnd_0_25, blm_1);
    __m256i rnd_0_26 = *(const LOAD_TYPE*)(src_origin + 26*src_stride);
    __m256i rnd_0_27 = *(const LOAD_TYPE*)(src_origin + 27*src_stride);
    __m256i shf_1_26 = _mm256_shuffle_epi8(rnd_0_26, shm_1);
    __m256i shf_1_27 = _mm256_shuffle_epi8(rnd_0_27, shm_1);
    __m256i rnd_1_26 = _mm256_blendv_epi8(rnd_0_26, shf_1_27, blm_1);
    __m256i rnd_1_27 = _mm256_blendv_epi8(shf_1_26, rnd_0_27, blm_1);
    __m256i rnd_0_28 = *(const LOAD_TYPE*)(src_origin + 28*src_stride);
    __m256i rnd_0_29 = *(const LOAD_TYPE*)(src_origin + 29*src_stride);
    __m256i shf_1_28 = _mm256_shuffle_epi8(rnd_0_28, shm_1);
    __m256i shf_1_29 = _mm256_shuffle_epi8(rnd_0_29, shm_1);
    __m256i rnd_1_28 = _mm256_blendv_epi8(rnd_0_28, shf_1_29, blm_1);
    __m256i rnd_1_29 = _mm256_blendv_epi8(shf_1_28, rnd_0_29, blm_1);
    __m256i rnd_0_30 = *(const LOAD_TYPE*)(src_origin + 30*src_stride);
    __m256i rnd_0_31 = *(const LOAD_TYPE*)(src_origin + 31*src_stride);
    __m256i shf_1_30 = _mm256_shuffle_epi8(rnd_0_30, shm_1);
    __m256i shf_1_31 = _mm256_shuffle_epi8(rnd_0_31, shm_1);
    _mm_prefetch(prf_origin+3*src_stride, _MM_HINT_NTA);
    __m256i rnd_1_30 = _mm256_blendv_epi8(rnd_0_30, shf_1_31, blm_1);
    __m256i rnd_1_31 = _mm256_blendv_epi8(shf_1_30, rnd_0_31, blm_1);
    const __m256i shm_2 = _mm256_load_si256((const __m256i *)SHUFFLE_MASK[1]);
    const __m256i blm_2 = _mm256_load_si256((const __m256i *)BLENDV_MASK[1]);
    __m256i shf_2_0 = _mm256_shuffle_epi8(rnd_1_0, shm_2);
    __m256i shf_2_2 = _mm256_shuffle_epi8(rnd_1_2, shm_2);
    __m256i rnd_2_0 = _mm256_blendv_epi8(rnd_1_0, shf_2_2, blm_2);
    __m256i rnd_2_2 = _mm256_blendv_epi8(shf_2_0, rnd_1_2, blm_2);
    __m256i shf_2_1 = _mm256_shuffle_epi8(rnd_1_1, shm_2);
    __m256i shf_2_3 = _mm256_shuffle_epi8(rnd_1_3, shm_2);
    __m256i rnd_2_1 = _mm256_blendv_epi8(rnd_1_1, shf_2_3, blm_2);
    __m256i rnd_2_3 = _mm256_blendv_epi8(shf_2_1, rnd_1_3, blm_2);
    __m256i shf_2_4 = _mm256_shuffle_epi8(rnd_1_4, shm_2);
    __m256i shf_2_6 = _mm256_shuffle_epi8(rnd_1_6, shm_2);
    __m256i rnd_2_4 = _mm256_blendv_epi8(rnd_1_4, shf_2_6, blm_2);
    __m256i rnd_2_6 = _mm256_blendv_epi8(shf_2_4, rnd_1_6, blm_2);
    __m256i shf_2_5 = _mm256_shuffle_epi8(rnd_1_5, shm_2);
    __m256i shf_2_7 = _mm256_shuffle_epi8(rnd_1_7, shm_2);
    __m256i rnd_2_5 = _mm256_blendv_epi8(rnd_1_5, shf_2_7, blm_2);
    __m256i rnd_2_7 = _mm256_blendv_epi8(shf_2_5, rnd_1_7, blm_2);
    __m256i shf_2_8 = _mm256_shuffle_epi8(rnd_1_8, shm_2);
    __m256i shf_2_10 = _mm256_shuffle_epi8(rnd_1_10, shm_2);
    __m256i rnd_2_8 = _mm256_blendv_epi8(rnd_1_8, shf_2_10, blm_2);
    __m256i rnd_2_10 = _mm256_blendv_epi8(shf_2_8, rnd_1_10, blm_2);
    _mm_prefetch(prf_origin+4*src_stride, _MM_HINT_NTA);
    __m256i shf_2_9 = _mm256_shuffle_epi8(rnd_1_9, shm_2);
    __m256i shf_2_11 = _mm256_shuffle_epi8(rnd_1_11, shm_2);
    __m256i rnd_2_9 = _mm256_blendv_epi8(rnd_1_9, shf_2_11, blm_2);
    __m256i rnd_2_11 = _mm256_blendv_epi8(shf_2_9, rnd_1_11, blm_2);
    __m256i shf_2_12 = _mm256_shuffle_epi8(rnd_1_12, shm_2);
    __m256i shf_2_14 = _mm256_shuffle_epi8(rnd_1_14, shm_2);
    __m256i rnd_2_12 = _mm256_blendv_epi8(rnd_1_12, shf_2_14, blm_2);
    __m256i rnd_2_14 = _mm256_blendv_epi8(shf_2_12, rnd_1_14, blm_2);
    __m256i shf_2_13 = _mm256_shuffle_epi8(rnd_1_13, shm_2);
    __m256i shf_2_15 = _mm256_shuffle_epi8(rnd_1_15, shm_2);
    __m256i rnd_2_13 = _mm256_blendv_epi8(rnd_1_13, shf_2_15, blm_2);
    __m256i rnd_2_15 = _mm256_blendv_epi8(shf_2_13, rnd_1_15, blm_2);
    __m256i shf_2_16 = _mm256_shuffle_epi8(rnd_1_16, shm_2);
    __m256i shf_2_18 = _mm256_shuffle_epi8(rnd_1_18, shm_2);
    __m256i rnd_2_16 = _mm256_blendv_epi8(rnd_1_16, shf_2_18, blm_2);
    __m256i rnd_2_18 = _mm256_blendv_epi8(shf_2_16, rnd_1_18, blm_2);
    __m256i shf_2_17 = _mm256_shuffle_epi8(rnd_1_17, shm_2);
    __m256i shf_2_19 = _mm256_shuffle_epi8(rnd_1_19, shm_2);
    __m256i rnd_2_17 = _mm256_blendv_epi8(rnd_1_17, shf_2_19, blm_2);
    __m256i rnd_2_19 = _mm256_blendv_epi8(shf_2_17, rnd_1_19, blm_2);
    __m256i shf_2_20 = _mm256_shuffle_epi8(rnd_1_20, shm_2);
    __m256i shf_2_22 = _mm256_shuffle_epi8(rnd_1_22, shm_2);
    __m256i rnd_2_20 = _mm256_blendv_epi8(rnd_1_20, shf_2_22, blm_2);
    __m256i rnd_2_22 = _mm256_blendv_epi8(shf_2_20, rnd_1_22, blm_2);
    _mm_prefetch(prf_origin+5*src_stride, _MM_HINT_NTA);
    __m256i shf_2_21 = _mm256_shuffle_epi8(rnd_1_21, shm_2);
    __m256i shf_2_23 = _mm256_shuffle_epi8(rnd_1_23, shm_2);
    __m256i rnd_2_21 = _mm256_blendv_epi8(rnd_1_21, shf_2_23, blm_2);
    __m256i rnd_2_23 = _mm256_blendv_epi8(shf_2_21, rnd_1_23, blm_2);
    __m256i shf_2_24 = _mm256_shuffle_epi8(rnd_1_24, shm_2);
    __m256i shf_2_26 = _mm256_shuffle_epi8(rnd_1_26, shm_2);
    __m256i rnd_2_24 = _mm256_blendv_epi8(rnd_1_24, shf_2_26, blm_2);
    __m256i rnd_2_26 = _mm256_blendv_epi8(shf_2_24, rnd_1_26, blm_2);
    __m256i shf_2_25 = _mm256_shuffle_epi8(rnd_1_25, shm_2);
    __m256i shf_2_27 = _mm256_shuffle_epi8(rnd_1_27, shm_2);
    __m256i rnd_2_25 = _mm256_blendv_epi8(rnd_1_25, shf_2_27, blm_2);
    __m256i rnd_2_27 = _mm256_blendv_epi8(shf_2_25, rnd_1_27, blm_2);
    __m256i shf_2_28 = _mm256_shuffle_epi8(rnd_1_28, shm_2);
    __m256i shf_2_30 = _mm256_shuffle_epi8(rnd_1_30, shm_2);
    __m256i rnd_2_28 = _mm256_blendv_epi8(rnd_1_28, shf_2_30, blm_2);
    __m256i rnd_2_30 = _mm256_blendv_epi8(shf_2_28, rnd_1_30, blm_2);
    __m256i shf_2_29 = _mm256_shuffle_epi8(rnd_1_29, shm_2);
    __m256i shf_2_31 = _mm256_shuffle_epi8(rnd_1_31, shm_2);
    __m256i rnd_2_29 = _mm256_blendv_epi8(rnd_1_29, shf_2_31, blm_2);
    __m256i rnd_2_31 = _mm256_blendv_epi8(shf_2_29, rnd_1_31, blm_2);
    const __m256i shm_3 = _mm256_load_si256((const __m256i *)SHUFFLE_MASK[2]);
    const __m256i blm_3 = _mm256_load_si256((const __m256i *)BLENDV_MASK[2]);
    __m256i shf_3_0 = _mm256_shuffle_epi8(rnd_2_0, shm_3);
    __m256i shf_3_4 = _mm256_shuffle_epi8(rnd_2_4, shm_3);
    _mm_prefetch(prf_origin+6*src_stride, _MM_HINT_NTA);
    __m256i rnd_3_0 = _mm256_blendv_epi8(rnd_2_0, shf_3_4, blm_3);
    __m256i rnd_3_4 = _mm256_blendv_epi8(shf_3_0, rnd_2_4, blm_3);
    __m256i shf_3_1 = _mm256_shuffle_epi8(rnd_2_1, shm_3);
    __m256i shf_3_5 = _mm256_shuffle_epi8(rnd_2_5, shm_3);
    __m256i rnd_3_1 = _mm256_blendv_epi8(rnd_2_1, shf_3_5, blm_3);
    __m256i rnd_3_5 = _mm256_blendv_epi8(shf_3_1, rnd_2_5, blm_3);
    __m256i shf_3_2 = _mm256_shuffle_epi8(rnd_2_2, shm_3);
    __m256i shf_3_6 = _mm256_shuffle_epi8(rnd_2_6, shm_3);
    __m256i rnd_3_2 = _mm256_blendv_epi8(rnd_2_2, shf_3_6, blm_3);
    __m256i rnd_3_6 = _mm256_blendv_epi8(shf_3_2, rnd_2_6, blm_3);
    __m256i shf_3_3 = _mm256_shuffle_epi8(rnd_2_3, shm_3);
    __m256i shf_3_7 = _mm256_shuffle_epi8(rnd_2_7, shm_3);
    __m256i rnd_3_3 = _mm256_blendv_epi8(rnd_2_3, shf_3_7, blm_3);
    __m256i rnd_3_7 = _mm256_blendv_epi8(shf_3_3, rnd_2_7, blm_3);
    __m256i shf_3_8 = _mm256_shuffle_epi8(rnd_2_8, shm_3);
    __m256i shf_3_12 = _mm256_shuffle_epi8(rnd_2_12, shm_3);
    __m256i rnd_3_8 = _mm256_blendv_epi8(rnd_2_8, shf_3_12, blm_3);
    __m256i rnd_3_12 = _mm256_blendv_epi8(shf_3_8, rnd_2_12, blm_3);
    __m256i shf_3_9 = _mm256_shuffle_epi8(rnd_2_9, shm_3);
    __m256i shf_3_13 = _mm256_shuffle_epi8(rnd_2_13, shm_3);
    __m256i rnd_3_9 = _mm256_blendv_epi8(rnd_2_9, shf_3_13, blm_3);
    __m256i rnd_3_13 = _mm256_blendv_epi8(shf_3_9, rnd_2_13, blm_3);
    __m256i shf_3_10 = _mm256_shuffle_epi8(rnd_2_10, shm_3);
    __m256i shf_3_14 = _mm256_shuffle_epi8(rnd_2_14, shm_3);
    _mm_prefetch(prf_origin+7*src_stride, _MM_HINT_NTA);
    __m256i rnd_3_10 = _mm256_blendv_epi8(rnd_2_10, shf_3_14, blm_3);
    __m256i rnd_3_14 = _mm256_blendv_epi8(shf_3_10, rnd_2_14, blm_3);
    __m256i shf_3_11 = _mm256_shuffle_epi8(rnd_2_11, shm_3);
    __m256i shf_3_15 = _mm256_shuffle_epi8(rnd_2_15, shm_3);
    __m256i rnd_3_11 = _mm256_blendv_epi8(rnd_2_11, shf_3_15, blm_3);
    __m256i rnd_3_15 = _mm256_blendv_epi8(shf_3_11, rnd_2_15, blm_3);
    __m256i shf_3_16 = _mm256_shuffle_epi8(rnd_2_16, shm_3);
    __m256i shf_3_20 = _mm256_shuffle_epi8(rnd_2_20, shm_3);
    __m256i rnd_3_16 = _mm256_blendv_epi8(rnd_2_16, shf_3_20, blm_3);
    __m256i rnd_3_20 = _mm256_blendv_epi8(shf_3_16, rnd_2_20, blm_3);
    __m256i shf_3_17 = _mm256_shuffle_epi8(rnd_2_17, shm_3);
    __m256i shf_3_21 = _mm256_shuffle_epi8(rnd_2_21, shm_3);
    __m256i rnd_3_17 = _mm256_blendv_epi8(rnd_2_17, shf_3_21, blm_3);
    __m256i rnd_3_21 = _mm256_blendv_epi8(shf_3_17, rnd_2_21, blm_3);
    __m256i shf_3_18 = _mm256_shuffle_epi8(rnd_2_18, shm_3);
    __m256i shf_3_22 = _mm256_shuffle_epi8(rnd_2_22, shm_3);
    __m256i rnd_3_18 = _mm256_blendv_epi8(rnd_2_18, shf_3_22, blm_3);
    __m256i rnd_3_22 = _mm256_blendv_epi8(shf_3_18, rnd_2_22, blm_3);
    __m256i shf_3_19 = _mm256_shuffle_epi8(rnd_2_19, shm_3);
    __m256i shf_3_23 = _mm256_shuffle_epi8(rnd_2_23, shm_3);
    __m256i rnd_3_19 = _mm256_blendv_epi8(rnd_2_19, shf_3_23, blm_3);
    __m256i rnd_3_23 = _mm256_blendv_epi8(shf_3_19, rnd_2_23, blm_3);
    __m256i shf_3_24 = _mm256_shuffle_epi8(rnd_2_24, shm_3);
    __m256i shf_3_28 = _mm256_shuffle_epi8(rnd_2_28, shm_3);
    _mm_prefetch(prf_origin+8*src_stride, _MM_HINT_NTA);
    __m256i rnd_3_24 = _mm256_blendv_epi8(rnd_2_24, shf_3_28, blm_3);
    __m256i rnd_3_28 = _mm256_blendv_epi8(shf_3_24, rnd_2_28, blm_3);
    __m256i shf_3_25 = _mm256_shuffle_epi8(rnd_2_25, shm_3);
    __m256i shf_3_29 = _mm256_shuffle_epi8(rnd_2_29, shm_3);
    __m256i rnd_3_25 = _mm256_blendv_epi8(rnd_2_25, shf_3_29, blm_3);
    __m256i rnd_3_29 = _mm256_blendv_epi8(shf_3_25, rnd_2_29, blm_3);
    __m256i shf_3_26 = _mm256_shuffle_epi8(rnd_2_26, shm_3);
    __m256i shf_3_30 = _mm256_shuffle_epi8(rnd_2_30, shm_3);
    __m256i rnd_3_26 = _mm256_blendv_epi8(rnd_2_26, shf_3_30, blm_3);
    __m256i rnd_3_30 = _mm256_blendv_epi8(shf_3_26, rnd_2_30, blm_3);
    __m256i shf_3_27 = _mm256_shuffle_epi8(rnd_2_27, shm_3);
    __m256i shf_3_31 = _mm256_shuffle_epi8(rnd_2_31, shm_3);
    __m256i rnd_3_27 = _mm256_blendv_epi8(rnd_2_27, shf_3_31, blm_3);
    __m256i rnd_3_31 = _mm256_blendv_epi8(shf_3_27, rnd_2_31, blm_3);
    const __m256i shm_4 = _mm256_load_si256((const __m256i *)SHUFFLE_MASK[3]);
    const __m256i blm_4 = _mm256_load_si256((const __m256i *)BLENDV_MASK[3]);
    __m256i shf_4_0 = _mm256_shuffle_epi8(rnd_3_0, shm_4);
    __m256i shf_4_8 = _mm256_shuffle_epi8(rnd_3_8, shm_4);
    __m256i rnd_4_0 = _mm256_blendv_epi8(rnd_3_0, shf_4_8, blm_4);
    __m256i rnd_4_8 = _mm256_blendv_epi8(shf_4_0, rnd_3_8, blm_4);
    __m256i shf_4_1 = _mm256_shuffle_epi8(rnd_3_1, shm_4);
    __m256i shf_4_9 = _mm256_shuffle_epi8(rnd_3_9, shm_4);
    __m256i rnd_4_1 = _mm256_blendv_epi8(rnd_3_1, shf_4_9, blm_4);
    __m256i rnd_4_9 = _mm256_blendv_epi8(shf_4_1, rnd_3_9, blm_4);
    _mm_prefetch(prf_origin+9*src_stride, _MM_HINT_NTA);
    __m256i shf_4_2 = _mm256_shuffle_epi8(rnd_3_2, shm_4);
    __m256i shf_4_10 = _mm256_shuffle_epi8(rnd_3_10, shm_4);
    __m256i rnd_4_2 = _mm256_blendv_epi8(rnd_3_2, shf_4_10, blm_4);
    __m256i rnd_4_10 = _mm256_blendv_epi8(shf_4_2, rnd_3_10, blm_4);
    __m256i shf_4_3 = _mm256_shuffle_epi8(rnd_3_3, shm_4);
    __m256i shf_4_11 = _mm256_shuffle_epi8(rnd_3_11, shm_4);
    __m256i rnd_4_3 = _mm256_blendv_epi8(rnd_3_3, shf_4_11, blm_4);
    __m256i rnd_4_11 = _mm256_blendv_epi8(shf_4_3, rnd_3_11, blm_4);
    __m256i shf_4_4 = _mm256_shuffle_epi8(rnd_3_4, shm_4);
    __m256i shf_4_12 = _mm256_shuffle_epi8(rnd_3_12, shm_4);
    __m256i rnd_4_4 = _mm256_blendv_epi8(rnd_3_4, shf_4_12, blm_4);
    __m256i rnd_4_12 = _mm256_blendv_epi8(shf_4_4, rnd_3_12, blm_4);
    __m256i shf_4_5 = _mm256_shuffle_epi8(rnd_3_5, shm_4);
    __m256i shf_4_13 = _mm256_shuffle_epi8(rnd_3_13, shm_4);
    __m256i rnd_4_5 = _mm256_blendv_epi8(rnd_3_5, shf_4_13, blm_4);
    __m256i rnd_4_13 = _mm256_blendv_epi8(shf_4_5, rnd_3_13, blm_4);
    __m256i shf_4_6 = _mm256_shuffle_epi8(rnd_3_6, shm_4);
    __m256i shf_4_14 = _mm256_shuffle_epi8(rnd_3_14, shm_4);
    __m256i rnd_4_6 = _mm256_blendv_epi8(rnd_3_6, shf_4_14, blm_4);
    __m256i rnd_4_14 = _mm256_blendv_epi8(shf_4_6, rnd_3_14, blm_4);
    __m256i shf_4_7 = _mm256_shuffle_epi8(rnd_3_7, shm_4);
    __m256i shf_4_15 = _mm256_shuffle_epi8(rnd_3_15, shm_4);
    __m256i rnd_4_7 = _mm256_blendv_epi8(rnd_3_7, shf_4_15, blm_4);
    __m256i rnd_4_15 = _mm256_blendv_epi8(shf_4_7, rnd_3_15, blm_4);
    _mm_prefetch(prf_origin+10*src_stride, _MM_HINT_NTA);
    __m256i shf_4_16 = _mm256_shuffle_epi8(rnd_3_16, shm_4);
    __m256i shf_4_24 = _mm256_shuffle_epi8(rnd_3_24, shm_4);
    __m256i rnd_4_16 = _mm256_blendv_epi8(rnd_3_16, shf_4_24, blm_4);
    __m256i rnd_4_24 = _mm256_blendv_epi8(shf_4_16, rnd_3_24, blm_4);
    __m256i shf_4_17 = _mm256_shuffle_epi8(rnd_3_17, shm_4);
    __m256i shf_4_25 = _mm256_shuffle_epi8(rnd_3_25, shm_4);
    __m256i rnd_4_17 = _mm256_blendv_epi8(rnd_3_17, shf_4_25, blm_4);
    __m256i rnd_4_25 = _mm256_blendv_epi8(shf_4_17, rnd_3_25, blm_4);
    __m256i shf_4_18 = _mm256_shuffle_epi8(rnd_3_18, shm_4);
    __m256i shf_4_26 = _mm256_shuffle_epi8(rnd_3_26, shm_4);
    __m256i rnd_4_18 = _mm256_blendv_epi8(rnd_3_18, shf_4_26, blm_4);
    __m256i rnd_4_26 = _mm256_blendv_epi8(shf_4_18, rnd_3_26, blm_4);
    __m256i shf_4_19 = _mm256_shuffle_epi8(rnd_3_19, shm_4);
    __m256i shf_4_27 = _mm256_shuffle_epi8(rnd_3_27, shm_4);
    __m256i rnd_4_19 = _mm256_blendv_epi8(rnd_3_19, shf_4_27, blm_4);
    __m256i rnd_4_27 = _mm256_blendv_epi8(shf_4_19, rnd_3_27, blm_4);
    __m256i shf_4_20 = _mm256_shuffle_epi8(rnd_3_20, shm_4);
    __m256i shf_4_28 = _mm256_shuffle_epi8(rnd_3_28, shm_4);
    __m256i rnd_4_20 = _mm256_blendv_epi8(rnd_3_20, shf_4_28, blm_4);
    __m256i rnd_4_28 = _mm256_blendv_epi8(shf_4_20, rnd_3_28, blm_4);
    __m256i shf_4_21 = _mm256_shuffle_epi8(rnd_3_21, shm_4);
    __m256i shf_4_29 = _mm256_shuffle_epi8(rnd_3_29, shm_4);
    __m256i rnd_4_21 = _mm256_blendv_epi8(rnd_3_21, shf_4_29, blm_4);
    __m256i rnd_4_29 = _mm256_blendv_epi8(shf_4_21, rnd_3_29, blm_4);
    _mm_prefetch(prf_origin+11*src_stride, _MM_HINT_NTA);
    __m256i shf_4_22 = _mm256_shuffle_epi8(rnd_3_22, shm_4);
    __m256i shf_4_30 = _mm256_shuffle_epi8(rnd_3_30, shm_4);
    __m256i rnd_4_22 = _mm256_blendv_epi8(rnd_3_22, shf_4_30, blm_4);
    __m256i rnd_4_30 = _mm256_blendv_epi8(shf_4_22, rnd_3_30, blm_4);
    __m256i shf_4_23 = _mm256_shuffle_epi8(rnd_3_23, shm_4);
    __m256i shf_4_31 = _mm256_shuffle_epi8(rnd_3_31, shm_4);
    __m256i rnd_4_23 = _mm256_blendv_epi8(rnd_3_23, shf_4_31, blm_4);
    __m256i rnd_4_31 = _mm256_blendv_epi8(shf_4_23, rnd_3_31, blm_4);
    const __m256i blm_5 = _mm256_load_si256((const __m256i *)BLENDV_MASK[4]);
    __m256i shf_5_0 = _mm256_permute2x128_si256(rnd_4_0, rnd_4_0, 0x01);
    __m256i shf_5_16 = _mm256_permute2x128_si256(rnd_4_16, rnd_4_16, 0x01);
    __m256i rnd_5_0 = _mm256_blendv_epi8(rnd_4_0, shf_5_16, blm_5);
    __m256i rnd_5_16 = _mm256_blendv_epi8(shf_5_0, rnd_4_16, blm_5);
    *(STORE_TYPE*)(dst_origin + 0*dst_stride) = rnd_5_0;
    *(STORE_TYPE*)(dst_origin + 16*dst_stride) = rnd_5_16;
    __m256i shf_5_1 = _mm256_permute2x128_si256(rnd_4_1, rnd_4_1, 0x01);
    __m256i shf_5_17 = _mm256_permute2x128_si256(rnd_4_17, rnd_4_17, 0x01);
    __m256i rnd_5_1 = _mm256_blendv_epi8(rnd_4_1, shf_5_17, blm_5);
    __m256i rnd_5_17 = _mm256_blendv_epi8(shf_5_1, rnd_4_17, blm_5);
    *(STORE_TYPE*)(dst_origin + 1*dst_stride) = rnd_5_1;
    *(STORE_TYPE*)(dst_origin + 17*dst_stride) = rnd_5_17;
    __m256i shf_5_2 = _mm256_permute2x128_si256(rnd_4_2, rnd_4_2, 0x01);
    __m256i shf_5_18 = _mm256_permute2x128_si256(rnd_4_18, rnd_4_18, 0x01);
    __m256i rnd_5_2 = _mm256_blendv_epi8(rnd_4_2, shf_5_18, blm_5);
    _mm_prefetch(prf_origin+12*src_stride, _MM_HINT_NTA);
    __m256i rnd_5_18 = _mm256_blendv_epi8(shf_5_2, rnd_4_18, blm_5);
    *(STORE_TYPE*)(dst_origin + 2*dst_stride) = rnd_5_2;
    *(STORE_TYPE*)(dst_origin + 18*dst_stride) = rnd_5_18;
    __m256i shf_5_3 = _mm256_permute2x128_si256(rnd_4_3, rnd_4_3, 0x01);
    __m256i shf_5_19 = _mm256_permute2x128_si256(rnd_4_19, rnd_4_19, 0x01);
    __m256i rnd_5_3 = _mm256_blendv_epi8(rnd_4_3, shf_5_19, blm_5);
    __m256i rnd_5_19 = _mm256_blendv_epi8(shf_5_3, rnd_4_19, blm_5);
    *(STORE_TYPE*)(dst_origin + 3*dst_stride) = rnd_5_3;
    *(STORE_TYPE*)(dst_origin + 19*dst_stride) = rnd_5_19;
    __m256i shf_5_4 = _mm256_permute2x128_si256(rnd_4_4, rnd_4_4, 0x01);
    __m256i shf_5_20 = _mm256_permute2x128_si256(rnd_4_20, rnd_4_20, 0x01);
    __m256i rnd_5_4 = _mm256_blendv_epi8(rnd_4_4, shf_5_20, blm_5);
    __m256i rnd_5_20 = _mm256_blendv_epi8(shf_5_4, rnd_4_20, blm_5);
    *(STORE_TYPE*)(dst_origin + 4*dst_stride) = rnd_5_4;
    *(STORE_TYPE*)(dst_origin + 20*dst_stride) = rnd_5_20;
    __m256i shf_5_5 = _mm256_permute2x128_si256(rnd_4_5, rnd_4_5, 0x01);
    __m256i shf_5_21 = _mm256_permute2x128_si256(rnd_4_21, rnd_4_21, 0x01);
    __m256i rnd_5_5 = _mm256_blendv_epi8(rnd_4_5, shf_5_21, blm_5);
    __m256i rnd_5_21 = _mm256_blendv_epi8(shf_5_5, rnd_4_21, blm_5);
    *(STORE_TYPE*)(dst_origin + 5*dst_stride) = rnd_5_5;
    *(STORE_TYPE*)(dst_origin + 21*dst_stride) = rnd_5_21;
    __m256i shf_5_6 = _mm256_permute2x128_si256(rnd_4_6, rnd_4_6, 0x01);
    __m256i shf_5_22 = _mm256_permute2x128_si256(rnd_4_22, rnd_4_22, 0x01);
    __m256i rnd_5_6 = _mm256_blendv_epi8(rnd_4_6, shf_5_22, blm_5);
    _mm_prefetch(prf_origin+13*src_stride, _MM_HINT_NTA);
    __m256i rnd_5_22 = _mm256_blendv_epi8(shf_5_6, rnd_4_22, blm_5);
    *(STORE_TYPE*)(dst_origin + 6*dst_stride) = rnd_5_6;
    *(STORE_TYPE*)(dst_origin + 22*dst_stride) = rnd_5_22;
    __m256i shf_5_7 = _mm256_permute2x128_si256(rnd_4_7, rnd_4_7, 0x01);
    __m256i shf_5_23 = _mm256_permute2x128_si256(rnd_4_23, rnd_4_23, 0x01);
    __m256i rnd_5_7 = _mm256_blendv_epi8(rnd_4_7, shf_5_23, blm_5);
    __m256i rnd_5_23 = _mm256_blendv_epi8(shf_5_7, rnd_4_23, blm_5);
    *(STORE_TYPE*)(dst_origin + 7*dst_stride) = rnd_5_7;
    *(STORE_TYPE*)(dst_origin + 23*dst_stride) = rnd_5_23;
    __m256i shf_5_8 = _mm256_permute2x128_si256(rnd_4_8, rnd_4_8, 0x01);
    __m256i shf_5_24 = _mm256_permute2x128_si256(rnd_4_24, rnd_4_24, 0x01);
    __m256i rnd_5_8 = _mm256_blendv_epi8(rnd_4_8, shf_5_24, blm_5);
    __m256i rnd_5_24 = _mm256_blendv_epi8(shf_5_8, rnd_4_24, blm_5);
    *(STORE_TYPE*)(dst_origin + 8*dst_stride) = rnd_5_8;
    *(STORE_TYPE*)(dst_origin + 24*dst_stride) = rnd_5_24;
    __m256i shf_5_9 = _mm256_permute2x128_si256(rnd_4_9, rnd_4_9, 0x01);
    __m256i shf_5_25 = _mm256_permute2x128_si256(rnd_4_25, rnd_4_25, 0x01);
    __m256i rnd_5_9 = _mm256_blendv_epi8(rnd_4_9, shf_5_25, blm_5);
    __m256i rnd_5_25 = _mm256_blendv_epi8(shf_5_9, rnd_4_25, blm_5);
    *(STORE_TYPE*)(dst_origin + 9*dst_stride) = rnd_5_9;
    *(STORE_TYPE*)(dst_origin + 25*dst_stride) = rnd_5_25;
    __m256i shf_5_10 = _mm256_permute2x128_si256(rnd_4_10, rnd_4_10, 0x01);
    __m256i shf_5_26 = _mm256_permute2x128_si256(rnd_4_26, rnd_4_26, 0x01);
    __m256i rnd_5_10 = _mm256_blendv_epi8(rnd_4_10, shf_5_26, blm_5);
    _mm_prefetch(prf_origin+14*src_stride, _MM_HINT_NTA);
    __m256i rnd_5_26 = _mm256_blendv_epi8(shf_5_10, rnd_4_26, blm_5);
    *(STORE_TYPE*)(dst_origin + 10*dst_stride) = rnd_5_10;
    *(STORE_TYPE*)(dst_origin + 26*dst_stride) = rnd_5_26;
    __m256i shf_5_11 = _mm256_permute2x128_si256(rnd_4_11, rnd_4_11, 0x01);
    __m256i shf_5_27 = _mm256_permute2x128_si256(rnd_4_27, rnd_4_27, 0x01);
    __m256i rnd_5_11 = _mm256_blendv_epi8(rnd_4_11, shf_5_27, blm_5);
    __m256i rnd_5_27 = _mm256_blendv_epi8(shf_5_11, rnd_4_27, blm_5);
    *(STORE_TYPE*)(dst_origin + 11*dst_stride) = rnd_5_11;
    *(STORE_TYPE*)(dst_origin + 27*dst_stride) = rnd_5_27;
    __m256i shf_5_12 = _mm256_permute2x128_si256(rnd_4_12, rnd_4_12, 0x01);
    __m256i shf_5_28 = _mm256_permute2x128_si256(rnd_4_28, rnd_4_28, 0x01);
    __m256i rnd_5_12 = _mm256_blendv_epi8(rnd_4_12, shf_5_28, blm_5);
    __m256i rnd_5_28 = _mm256_blendv_epi8(shf_5_12, rnd_4_28, blm_5);
    *(STORE_TYPE*)(dst_origin + 12*dst_stride) = rnd_5_12;
    *(STORE_TYPE*)(dst_origin + 28*dst_stride) = rnd_5_28;
    __m256i shf_5_13 = _mm256_permute2x128_si256(rnd_4_13, rnd_4_13, 0x01);
    __m256i shf_5_29 = _mm256_permute2x128_si256(rnd_4_29, rnd_4_29, 0x01);
    __m256i rnd_5_13 = _mm256_blendv_epi8(rnd_4_13, shf_5_29, blm_5);
    __m256i rnd_5_29 = _mm256_blendv_epi8(shf_5_13, rnd_4_29, blm_5);
    *(STORE_TYPE*)(dst_origin + 13*dst_stride) = rnd_5_13;
    *(STORE_TYPE*)(dst_origin + 29*dst_stride) = rnd_5_29;
    __m256i shf_5_14 = _mm256_permute2x128_si256(rnd_4_14, rnd_4_14, 0x01);
    __m256i shf_5_30 = _mm256_permute2x128_si256(rnd_4_30, rnd_4_30, 0x01);
    __m256i rnd_5_14 = _mm256_blendv_epi8(rnd_4_14, shf_5_30, blm_5);
    _mm_prefetch(prf_origin+15*src_stride, _MM_HINT_NTA);
    __m256i rnd_5_30 = _mm256_blendv_epi8(shf_5_14, rnd_4_30, blm_5);
    *(STORE_TYPE*)(dst_origin + 14*dst_stride) = rnd_5_14;
    *(STORE_TYPE*)(dst_origin + 30*dst_stride) = rnd_5_30;
    __m256i shf_5_15 = _mm256_permute2x128_si256(rnd_4_15, rnd_4_15, 0x01);
    __m256i shf_5_31 = _mm256_permute2x128_si256(rnd_4_31, rnd_4_31, 0x01);
    __m256i rnd_5_15 = _mm256_blendv_epi8(rnd_4_15, shf_5_31, blm_5);
    __m256i rnd_5_31 = _mm256_blendv_epi8(shf_5_15, rnd_4_31, blm_5);
    *(STORE_TYPE*)(dst_origin + 15*dst_stride) = rnd_5_15;
    *(STORE_TYPE*)(dst_origin + 31*dst_stride) = rnd_5_31;
}

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
    
    for (uint32_t i = 0; i < 32; i+=2) {
        const __m256i t0 = _mm256_unpacklo_epi8(t[i+0], t[i+1]);
        const __m256i t1 = _mm256_unpackhi_epi8(t[i+0], t[i+1]);
        t[i+0] = t0;
        t[i+1] = t1;
    }


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

    for (uint32_t i = 0; i < 16; i++) {
        const __m256i t0 = _mm256_permute2x128_si256(t[i+0], t[i+16], 0b100000);
        const __m256i t1 = _mm256_permute2x128_si256(t[i+0], t[i+16], 0b110001);
        t[i+ 0] = t0;
        t[i+16] = t1;
    }

    for (uint32_t i = 0; i < 32; i++) {
        _mm256_storeu_si256((__m256i *)(dst_origin + i*dst_stride), t[i]);
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

