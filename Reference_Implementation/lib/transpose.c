/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
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

#include "transpose.h"

/// \param dst[out]: out data
/// \param src[in] input bytes 8x8 matrix
/// \param src_stride[in]: number of bytes between two rows in `src`
/// \param dst_stride[in]: number of bytes between two cols in `dst`
void matrix_transpose_8x8(uint8_t* dst,
                          const uint8_t* src,
                          const size_t src_stride,
                          const size_t dst_stride) {
    // load rows of src matrix
    const uint64_t a0 = *((uint64_t*)(src+0*src_stride));
    const uint64_t a1 = *((uint64_t*)(src+1*src_stride));
    const uint64_t a2 = *((uint64_t*)(src+2*src_stride));
    const uint64_t a3 = *((uint64_t*)(src+3*src_stride));
    const uint64_t a4 = *((uint64_t*)(src+4*src_stride));
    const uint64_t a5 = *((uint64_t*)(src+5*src_stride));
    const uint64_t a6 = *((uint64_t*)(src+6*src_stride));
    const uint64_t a7 = *((uint64_t*)(src+7*src_stride));

    // 2x2 block matrices
    const uint64_t b0 = (a0 & 0x00ff00ff00ff00ffULL) | ((a1 << 8) & 0xff00ff00ff00ff00ULL);
    const uint64_t b1 = (a1 & 0xff00ff00ff00ff00ULL) | ((a0 >> 8) & 0x00ff00ff00ff00ffULL);
    const uint64_t b2 = (a2 & 0x00ff00ff00ff00ffULL) | ((a3 << 8) & 0xff00ff00ff00ff00ULL);
    const uint64_t b3 = (a3 & 0xff00ff00ff00ff00ULL) | ((a2 >> 8) & 0x00ff00ff00ff00ffULL);
    const uint64_t b4 = (a4 & 0x00ff00ff00ff00ffULL) | ((a5 << 8) & 0xff00ff00ff00ff00ULL);
    const uint64_t b5 = (a5 & 0xff00ff00ff00ff00ULL) | ((a4 >> 8) & 0x00ff00ff00ff00ffULL);
    const uint64_t b6 = (a6 & 0x00ff00ff00ff00ffULL) | ((a7 << 8) & 0xff00ff00ff00ff00ULL);
    const uint64_t b7 = (a7 & 0xff00ff00ff00ff00ULL) | ((a6 >> 8) & 0x00ff00ff00ff00ffULL);

    // 4x4 block matrices
    const uint64_t c0 = (b0 & 0x0000ffff0000ffffULL) | ((b2 << 16) & 0xffff0000ffff0000ULL);
    const uint64_t c1 = (b1 & 0x0000ffff0000ffffULL) | ((b3 << 16) & 0xffff0000ffff0000ULL);
    const uint64_t c2 = (b2 & 0xffff0000ffff0000ULL) | ((b0 >> 16) & 0x0000ffff0000ffffULL);
    const uint64_t c3 = (b3 & 0xffff0000ffff0000ULL) | ((b1 >> 16) & 0x0000ffff0000ffffULL);
    const uint64_t c4 = (b4 & 0x0000ffff0000ffffULL) | ((b6 << 16) & 0xffff0000ffff0000ULL);
    const uint64_t c5 = (b5 & 0x0000ffff0000ffffULL) | ((b7 << 16) & 0xffff0000ffff0000ULL);
    const uint64_t c6 = (b6 & 0xffff0000ffff0000ULL) | ((b4 >> 16) & 0x0000ffff0000ffffULL);
    const uint64_t c7 = (b7 & 0xffff0000ffff0000ULL) | ((b5 >> 16) & 0x0000ffff0000ffffULL);

    // 8x8 block matrix
    const uint64_t d0 = (c0 & 0x00000000ffffffffULL) | ((c4 << 32) & 0xffffffff00000000ULL);
    const uint64_t d1 = (c1 & 0x00000000ffffffffULL) | ((c5 << 32) & 0xffffffff00000000ULL);
    const uint64_t d2 = (c2 & 0x00000000ffffffffULL) | ((c6 << 32) & 0xffffffff00000000ULL);
    const uint64_t d3 = (c3 & 0x00000000ffffffffULL) | ((c7 << 32) & 0xffffffff00000000ULL);
    const uint64_t d4 = (c4 & 0xffffffff00000000ULL) | ((c0 >> 32) & 0x00000000ffffffffULL);
    const uint64_t d5 = (c5 & 0xffffffff00000000ULL) | ((c1 >> 32) & 0x00000000ffffffffULL);
    const uint64_t d6 = (c6 & 0xffffffff00000000ULL) | ((c2 >> 32) & 0x00000000ffffffffULL);
    const uint64_t d7 = (c7 & 0xffffffff00000000ULL) | ((c3 >> 32) & 0x00000000ffffffffULL);

    // write to dst matrix
    *(uint64_t*)(dst + 0*dst_stride) = d0;
    *(uint64_t*)(dst + 1*dst_stride) = d1;
    *(uint64_t*)(dst + 2*dst_stride) = d2;
    *(uint64_t*)(dst + 3*dst_stride) = d3;
    *(uint64_t*)(dst + 4*dst_stride) = d4;
    *(uint64_t*)(dst + 5*dst_stride) = d5;
    *(uint64_t*)(dst + 6*dst_stride) = d6;
    *(uint64_t*)(dst + 7*dst_stride) = d7;
}

#if !(defined(USE_AVX2) || defined(USE_NEON) || defined(USE_AVX512))
/// \param dst[out]: out data
/// \param src[in] input bytes 8x8 matrix
/// \param src_stride[in]: number of bytes between two rows in `src`
/// \param dst_stride[in]: number of bytes between two cols in `dst`
void matrix_transpose_32x32(uint8_t* dst,
                            const uint8_t* src,
                            const size_t src_stride,
                            const size_t dst_stride) {
    for (size_t rw = 0; rw < 4; rw++) {
        for (size_t cw = 0; cw < 4; cw++) {
            const uint8_t *srcw_origin = src + (cw*src_stride + rw) * 8;
                  uint8_t *dstw_origin = dst + (rw*dst_stride + cw) * 8;
            matrix_transpose_8x8(dstw_origin, srcw_origin, src_stride, dst_stride);
        }
    }
}
#endif

/// \param dst[out]: output non-IS matrix: K \times N-K
/// \param src[in]: input non-IS matrix: K \times N-K
/// \param r[in]: number of rows in `src`
/// \param c[in]: number of cols in `src`
/// \param src_stride[in]: number of elements between two rows in `src`
/// \param dst_stride[in]: number of elements between two cols in `dst`
void matrix_transpose(uint8_t *dst,
                      const uint8_t *src,
                      const size_t r,
                      const size_t c,
                      const size_t dst_stride,
                      const size_t src_stride) {
#if !(defined(USE_AVX2) || defined(USE_NEON) || defined(USE_AVX512))
    for (uint64_t row = 0; row < r; row++) {
        for (uint64_t col = 0; col < c; col++) {
            dst[col*dst_stride + row] = src[row*src_stride + col];
        }
    }
#else
    const size_t bsize = 32;

    for (uint64_t rb = 0; rb < r / bsize; rb++) {
        for (uint64_t cb = 0; cb < c / bsize; cb++) {
            const uint8_t* src_origin = src + (rb*src_stride+cb)*bsize;
                  uint8_t* dst_origin = dst + (cb*dst_stride+rb)*bsize;
#if defined(USE_AVX512)
            // turns out: on some machines the avx512 is slower than
            // the avx2 implementation. Choose whatever is faster for you.
            // matrix_transpose_64x64(dst_origin, src_origin, n, n);
            matrix_transpose_32x32(dst_origin, src_origin, src_stride, dst_stride);
#else
            matrix_transpose_32x32(dst_origin, src_origin, src_stride, dst_stride);
#endif
        }
    }
#endif
}
