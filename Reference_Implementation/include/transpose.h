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

#pragma once
#include <stdint.h>

#if defined(USE_AVX2) || defined(USE_NEON)
#include <stdlib.h>

/// Only avx2 and neon implement this function
/// transposes a 32x32 matrix
/// \param dst_origin[out]: pointer to the out buffer
/// \param src_origin[in]: pointer to the src buffer
/// \param src_stride[in]: number of bytes between two rows.
/// \param dst_stride[in]: number of bytes between two cols.
void matrix_transpose_32x32(uint8_t* dst_origin,
                            const uint8_t* src_origin,
                            size_t src_stride,
                            size_t dst_stride);
#endif

#if defined(USE_AVX512)
#include <stdlib.h>

/// only avx512 implements this function
/// transposes a 64x64 matrix
/// \param dst_origin[out]: pointer to the out buffer
/// \param src_origin[in]: pointer to the src buffer
/// \param src_stride[in]: number of bytes between two rows.
/// \param dst_stride[in]: number of bytes between two cols.
void matrix_transpose_64x64(uint8_t* dst_origin,
                            const uint8_t* src_origin,
                            size_t src_stride,
                            size_t dst_stride);
#endif

/// transpose `src` int `dst`
/// \param dst[out]: pointer to the output matrix
/// \param src[in]: pointer to the input matrix
/// \param r[in]: number of rows in `src`
/// \param c[in]: number of cols in `src`
void matrix_transpose_opt(uint8_t *dst,
                          const uint8_t *src,
                          uint32_t r,
                          uint32_t c);
