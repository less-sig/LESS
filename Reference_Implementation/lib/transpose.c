/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Floyd Zweydinge <zweydfg8+github@rub.de>
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
#include "parameters.h"

/// \param dst[out]: out  data
/// \param src[in] input bytes 8x8 matrix
/// \param src_stride[in] in bytes
/// \param dst_stride[in] in bytes
void matrix_transpose8x8(uint8_t* dst,
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

/// assumes max 8 rows in the input matrix
/// assumes that the output matrix has n columns
/// \param dst
/// \param src
/// \param n   number of columns
static inline void matrix_transpose8xN(uint8_t *dst, 
                                       const uint8_t *src,
                                       const uint32_t n) {
    const uint32_t bsize = 8;
    uint64_t rb = 0;
    for (; rb < n / bsize; rb++) {
            const uint8_t *srcb_origin = src + ( 0 * n + rb) * bsize;
                  uint8_t *dstb_origin = dst + (rb * n +  0) * bsize;
        matrix_transpose8x8(dstb_origin, srcb_origin, n, n);
    }

    rb *= bsize;
    for (; rb < n; rb++) {
        for (uint32_t j = 0; j < 8; j++) {
            const uint8_t t = src[j*n + rb];
            dst[rb*n + j] = t;
        }
    }
}

/// assumes max 8 cols in the input matrix
/// assumes that the output matrix has n columns
/// \param dst
/// \param src
/// \param n   number of columns
static inline void matrix_transposeNx8(uint8_t *dst,
                                       const uint8_t *src,
                                       const uint32_t n) {
    const uint32_t bsize = 8;
    uint64_t cb = 0;
    for (; cb < n / bsize; cb++) {
            const uint8_t *srcb_origin = src + (cb * n +  0) * bsize;
                  uint8_t *dstb_origin = dst + ( 0 * n + cb) * bsize;
        matrix_transpose8x8(dstb_origin, srcb_origin, n, n);
    }

    cb *= bsize;
    for (; cb < n; cb++) {
        for (uint32_t j = 0; j < 8; j++) {
            const uint8_t t = src[cb*n + j];
            dst[j*n + cb] = t;
        }
    }
}


/// Compute origin of the 64-block next to (rb, cb) in row-major order
/// NOTE: internal function. Do no call directly.
const uint8_t* next_block(const uint8_t *src,
                          uint64_t rb,
                          uint64_t cb,
                          const size_t n) {
    uint64_t cb1 = cb + 1;
    uint64_t rb1 = rb;
    if (cb1 == n/64) {
        rb1 += 1;
        cb1 = 0;
    }

    return src + (rb1*n + cb1) * 64;
}

/// \param dst[out]: output non-IS matrix: K \times N-K
/// \param src[in]: input non-IS matrix: K \times N-K
/// \param r nr cols
/// \param c nr rows
void matrix_transpose_opt(uint8_t *dst,
                          const uint8_t *src,
                          const uint32_t r,
                          const uint32_t c) {
    // small block size
    const size_t small = 8;

#if defined(USE_AVX2) || defined(USE_NEON)
    // big block size
    const size_t bsize = 32;
#else
    const size_t bsize = 64;
#endif

    const size_t src_stride = K_pad;
    const size_t dst_stride = N_K_pad;

    if ((c < bsize) || (r < bsize)) {
        if (c <= small) {
            matrix_transposeNx8(dst, src, r);
            return;
        }

        if (r <= small) {
            matrix_transpose8xN(dst, src, c);
            return;
        }
        for (uint32_t i = 0; i < r; i++) {
            for (uint32_t j = 0; j < c; j++) {
                dst[j*dst_stride + i] = src[i*src_stride + j];
            }
        }

        return ;
    }

    uint64_t rb = 0;
    for (; rb < c / bsize; rb++) {
        for (uint64_t cb = 0; cb < c / bsize; cb++) {
#if defined(USE_AVX2) || defined(USE_NEON)
            const uint8_t* prf_origin = next_block(src, rb, cb, src_stride);
            const uint8_t* src_origin = src + (rb*src_stride+cb)*bsize;
                  uint8_t* dst_origin = dst + (cb*dst_stride+rb)*bsize;
            const uint32_t n = src_stride;

            matrix_transpose_32x32(dst_origin,                  src_origin,                  prf_origin,               n, n);
#else
            const uint8_t *srcb_origin = src + (rb*src_stride + cb) * bsize;
                  uint8_t *dstb_origin = dst + (cb*dst_stride + rb) * bsize;
            for (size_t rw = 0; rw < bsize / 8; rw++) {
                for (size_t cw = 0; cw < bsize / 8; cw++) {
                    const uint8_t *srcw_origin = srcb_origin + (cw*src_stride + rw) * 8;
                          uint8_t *dstw_origin = dstb_origin + (rw*dst_stride + cw) * 8;
                    matrix_transpose8x8(dstw_origin, srcw_origin, src_stride, dst_stride);
                }
            }
#endif
        }
    }

    const uint32_t rem = c % bsize;
    if (rem) {
        rb *= (64 / bsize);

        // solve the last columns
        for (uint32_t i = rb*bsize; i < c; i++) {
            for(uint32_t j = 0; j < c; j++) {
                dst[j*dst_stride + i] = src[i*src_stride + j];
            }
        }
        // solve the last rows
        for (uint32_t i = 0; i < c; i++) {
            for(uint32_t j = rb*bsize; j < c; j++) {
                dst[j*dst_stride + i] = src[i*src_stride + j];
            }
        }
    }
}
