#pragma once

#ifdef USE_AVX2

void matrix_transpose_32x32_avx2(uint8_t* dst_origin,
                                 const uint8_t* src_origin,
                                 const uint8_t* prf_origin,
                                 const size_t src_stride,
                                 const size_t dst_stride);

#endif

void matrix_transpose_opt(uint8_t *dst,
                          const uint8_t *src,
                          const uint32_t r,
                          const uint32_t c);
