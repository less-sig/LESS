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

/// NOTE: needed for `vec256_t` definition.
#include "fq_arith.h"

/// \param dst_origin[out]: output matrix
/// \param src_origin[in]: input matrix
/// \param prf_origin[in]: lookahead pointer to prefetch it (unused in neon impl.)
/// \param src_stride[in]: number of bytes (including alignment) in each row for source matrix
/// \param dst_stride[in]: number of bytes (including alignment) in each row for destination matrix
void matrix_transpose_32x32(uint8_t* dst_origin,
                           const uint8_t* src_origin,
                           const uint8_t* prf_origin,
                           const size_t src_stride,
                           const size_t dst_stride) {
    (void)prf_origin;
    const vec256_t rnd_0_0 = *((vec256_t *)(src_origin + 0*src_stride));
    const vec256_t rnd_0_1 = *((vec256_t *)(src_origin + 1*src_stride));
    vec256_t rnd_1_0; rnd_1_0.v[0] = vtrn1q_u8(rnd_0_0.v[0], rnd_0_1.v[0]); rnd_1_0.v[1] = vtrn1q_u8(rnd_0_0.v[1], rnd_0_1.v[1]);
    vec256_t rnd_1_1; rnd_1_1.v[0] = vtrn2q_u8(rnd_0_0.v[0], rnd_0_1.v[0]); rnd_1_1.v[1] = vtrn2q_u8(rnd_0_0.v[1], rnd_0_1.v[1]);
    const vec256_t rnd_0_2 = *((vec256_t *)(src_origin + 2*src_stride));
    const vec256_t rnd_0_3 = *((vec256_t *)(src_origin + 3*src_stride));
    vec256_t rnd_1_2; rnd_1_2.v[0] = vtrn1q_u8(rnd_0_2.v[0], rnd_0_3.v[0]); rnd_1_2.v[1] = vtrn1q_u8(rnd_0_2.v[1], rnd_0_3.v[1]);
    vec256_t rnd_1_3; rnd_1_3.v[0] = vtrn2q_u8(rnd_0_2.v[0], rnd_0_3.v[0]); rnd_1_3.v[1] = vtrn2q_u8(rnd_0_2.v[1], rnd_0_3.v[1]);
    const vec256_t rnd_0_4 = *((vec256_t *)(src_origin + 4*src_stride));
    const vec256_t rnd_0_5 = *((vec256_t *)(src_origin + 5*src_stride));
    vec256_t rnd_1_4; rnd_1_4.v[0] = vtrn1q_u8(rnd_0_4.v[0], rnd_0_5.v[0]); rnd_1_4.v[1] = vtrn1q_u8(rnd_0_4.v[1], rnd_0_5.v[1]);
    vec256_t rnd_1_5; rnd_1_5.v[0] = vtrn2q_u8(rnd_0_4.v[0], rnd_0_5.v[0]); rnd_1_5.v[1] = vtrn2q_u8(rnd_0_4.v[1], rnd_0_5.v[1]);
    const vec256_t rnd_0_6 = *((vec256_t *)(src_origin + 6*src_stride));
    const vec256_t rnd_0_7 = *((vec256_t *)(src_origin + 7*src_stride));
    vec256_t rnd_1_6; rnd_1_6.v[0] = vtrn1q_u8(rnd_0_6.v[0], rnd_0_7.v[0]); rnd_1_6.v[1] = vtrn1q_u8(rnd_0_6.v[1], rnd_0_7.v[1]);
    vec256_t rnd_1_7; rnd_1_7.v[0] = vtrn2q_u8(rnd_0_6.v[0], rnd_0_7.v[0]); rnd_1_7.v[1] = vtrn2q_u8(rnd_0_6.v[1], rnd_0_7.v[1]);
    const vec256_t rnd_0_8 = *((vec256_t *)(src_origin + 8*src_stride));
    const vec256_t rnd_0_9 = *((vec256_t *)(src_origin + 9*src_stride));
    vec256_t rnd_1_8; rnd_1_8.v[0] = vtrn1q_u8(rnd_0_8.v[0], rnd_0_9.v[0]); rnd_1_8.v[1] = vtrn1q_u8(rnd_0_8.v[1], rnd_0_9.v[1]);
    vec256_t rnd_1_9; rnd_1_9.v[0] = vtrn2q_u8(rnd_0_8.v[0], rnd_0_9.v[0]); rnd_1_9.v[1] = vtrn2q_u8(rnd_0_8.v[1], rnd_0_9.v[1]);
    const vec256_t rnd_0_10 = *((vec256_t *)(src_origin + 10*src_stride));
    const vec256_t rnd_0_11 = *((vec256_t *)(src_origin + 11*src_stride));
    vec256_t rnd_1_10; rnd_1_10.v[0] = vtrn1q_u8(rnd_0_10.v[0], rnd_0_11.v[0]); rnd_1_10.v[1] = vtrn1q_u8(rnd_0_10.v[1], rnd_0_11.v[1]);
    vec256_t rnd_1_11; rnd_1_11.v[0] = vtrn2q_u8(rnd_0_10.v[0], rnd_0_11.v[0]); rnd_1_11.v[1] = vtrn2q_u8(rnd_0_10.v[1], rnd_0_11.v[1]);
    const vec256_t rnd_0_12 = *((vec256_t *)(src_origin + 12*src_stride));
    const vec256_t rnd_0_13 = *((vec256_t *)(src_origin + 13*src_stride));
    vec256_t rnd_1_12; rnd_1_12.v[0] = vtrn1q_u8(rnd_0_12.v[0], rnd_0_13.v[0]); rnd_1_12.v[1] = vtrn1q_u8(rnd_0_12.v[1], rnd_0_13.v[1]);
    vec256_t rnd_1_13; rnd_1_13.v[0] = vtrn2q_u8(rnd_0_12.v[0], rnd_0_13.v[0]); rnd_1_13.v[1] = vtrn2q_u8(rnd_0_12.v[1], rnd_0_13.v[1]);
    const vec256_t rnd_0_14 = *((vec256_t *)(src_origin + 14*src_stride));
    const vec256_t rnd_0_15 = *((vec256_t *)(src_origin + 15*src_stride));
    vec256_t rnd_1_14; rnd_1_14.v[0] = vtrn1q_u8(rnd_0_14.v[0], rnd_0_15.v[0]); rnd_1_14.v[1] = vtrn1q_u8(rnd_0_14.v[1], rnd_0_15.v[1]);
    vec256_t rnd_1_15; rnd_1_15.v[0] = vtrn2q_u8(rnd_0_14.v[0], rnd_0_15.v[0]); rnd_1_15.v[1] = vtrn2q_u8(rnd_0_14.v[1], rnd_0_15.v[1]);
    const vec256_t rnd_0_16 = *((vec256_t *)(src_origin + 16*src_stride));
    const vec256_t rnd_0_17 = *((vec256_t *)(src_origin + 17*src_stride));
    vec256_t rnd_1_16; rnd_1_16.v[0] = vtrn1q_u8(rnd_0_16.v[0], rnd_0_17.v[0]); rnd_1_16.v[1] = vtrn1q_u8(rnd_0_16.v[1], rnd_0_17.v[1]);
    vec256_t rnd_1_17; rnd_1_17.v[0] = vtrn2q_u8(rnd_0_16.v[0], rnd_0_17.v[0]); rnd_1_17.v[1] = vtrn2q_u8(rnd_0_16.v[1], rnd_0_17.v[1]);
    const vec256_t rnd_0_18 = *((vec256_t *)(src_origin + 18*src_stride));
    const vec256_t rnd_0_19 = *((vec256_t *)(src_origin + 19*src_stride));
    vec256_t rnd_1_18; rnd_1_18.v[0] = vtrn1q_u8(rnd_0_18.v[0], rnd_0_19.v[0]); rnd_1_18.v[1] = vtrn1q_u8(rnd_0_18.v[1], rnd_0_19.v[1]);
    vec256_t rnd_1_19; rnd_1_19.v[0] = vtrn2q_u8(rnd_0_18.v[0], rnd_0_19.v[0]); rnd_1_19.v[1] = vtrn2q_u8(rnd_0_18.v[1], rnd_0_19.v[1]);
    const vec256_t rnd_0_20 = *((vec256_t *)(src_origin + 20*src_stride));
    const vec256_t rnd_0_21 = *((vec256_t *)(src_origin + 21*src_stride));
    vec256_t rnd_1_20; rnd_1_20.v[0] = vtrn1q_u8(rnd_0_20.v[0], rnd_0_21.v[0]); rnd_1_20.v[1] = vtrn1q_u8(rnd_0_20.v[1], rnd_0_21.v[1]);
    vec256_t rnd_1_21; rnd_1_21.v[0] = vtrn2q_u8(rnd_0_20.v[0], rnd_0_21.v[0]); rnd_1_21.v[1] = vtrn2q_u8(rnd_0_20.v[1], rnd_0_21.v[1]);
    const vec256_t rnd_0_22 = *((vec256_t *)(src_origin + 22*src_stride));
    const vec256_t rnd_0_23 = *((vec256_t *)(src_origin + 23*src_stride));
    vec256_t rnd_1_22; rnd_1_22.v[0] = vtrn1q_u8(rnd_0_22.v[0], rnd_0_23.v[0]); rnd_1_22.v[1] = vtrn1q_u8(rnd_0_22.v[1], rnd_0_23.v[1]);
    vec256_t rnd_1_23; rnd_1_23.v[0] = vtrn2q_u8(rnd_0_22.v[0], rnd_0_23.v[0]); rnd_1_23.v[1] = vtrn2q_u8(rnd_0_22.v[1], rnd_0_23.v[1]);
    const vec256_t rnd_0_24 = *((vec256_t *)(src_origin + 24*src_stride));
    const vec256_t rnd_0_25 = *((vec256_t *)(src_origin + 25*src_stride));
    vec256_t rnd_1_24; rnd_1_24.v[0] = vtrn1q_u8(rnd_0_24.v[0], rnd_0_25.v[0]); rnd_1_24.v[1] = vtrn1q_u8(rnd_0_24.v[1], rnd_0_25.v[1]);
    vec256_t rnd_1_25; rnd_1_25.v[0] = vtrn2q_u8(rnd_0_24.v[0], rnd_0_25.v[0]); rnd_1_25.v[1] = vtrn2q_u8(rnd_0_24.v[1], rnd_0_25.v[1]);
    const vec256_t rnd_0_26 = *((vec256_t *)(src_origin + 26*src_stride));
    const vec256_t rnd_0_27 = *((vec256_t *)(src_origin + 27*src_stride));
    vec256_t rnd_1_26; rnd_1_26.v[0] = vtrn1q_u8(rnd_0_26.v[0], rnd_0_27.v[0]); rnd_1_26.v[1] = vtrn1q_u8(rnd_0_26.v[1], rnd_0_27.v[1]);
    vec256_t rnd_1_27; rnd_1_27.v[0] = vtrn2q_u8(rnd_0_26.v[0], rnd_0_27.v[0]); rnd_1_27.v[1] = vtrn2q_u8(rnd_0_26.v[1], rnd_0_27.v[1]);
    const vec256_t rnd_0_28 = *((vec256_t *)(src_origin + 28*src_stride));
    const vec256_t rnd_0_29 = *((vec256_t *)(src_origin + 29*src_stride));
    vec256_t rnd_1_28; rnd_1_28.v[0] = vtrn1q_u8(rnd_0_28.v[0], rnd_0_29.v[0]); rnd_1_28.v[1] = vtrn1q_u8(rnd_0_28.v[1], rnd_0_29.v[1]);
    vec256_t rnd_1_29; rnd_1_29.v[0] = vtrn2q_u8(rnd_0_28.v[0], rnd_0_29.v[0]); rnd_1_29.v[1] = vtrn2q_u8(rnd_0_28.v[1], rnd_0_29.v[1]);
    const vec256_t rnd_0_30 = *((vec256_t *)(src_origin + 30*src_stride));
    const vec256_t rnd_0_31 = *((vec256_t *)(src_origin + 31*src_stride));
    vec256_t rnd_1_30; rnd_1_30.v[0] = vtrn1q_u8(rnd_0_30.v[0], rnd_0_31.v[0]); rnd_1_30.v[1] = vtrn1q_u8(rnd_0_30.v[1], rnd_0_31.v[1]);
    vec256_t rnd_1_31; rnd_1_31.v[0] = vtrn2q_u8(rnd_0_30.v[0], rnd_0_31.v[0]); rnd_1_31.v[1] = vtrn2q_u8(rnd_0_30.v[1], rnd_0_31.v[1]);

    vec256_t rnd_2_0; rnd_2_0.v[0] = vtrn1q_u16(rnd_1_0.v[0], rnd_1_2.v[0]); rnd_2_0.v[1] = vtrn1q_u16(rnd_1_0.v[1], rnd_1_2.v[1]);
    vec256_t rnd_2_2; rnd_2_2.v[0] = vtrn2q_u16(rnd_1_0.v[0], rnd_1_2.v[0]); rnd_2_2.v[1] = vtrn2q_u16(rnd_1_0.v[1], rnd_1_2.v[1]);
    vec256_t rnd_2_1; rnd_2_1.v[0] = vtrn1q_u16(rnd_1_1.v[0], rnd_1_3.v[0]); rnd_2_1.v[1] = vtrn1q_u16(rnd_1_1.v[1], rnd_1_3.v[1]);
    vec256_t rnd_2_3; rnd_2_3.v[0] = vtrn2q_u16(rnd_1_1.v[0], rnd_1_3.v[0]); rnd_2_3.v[1] = vtrn2q_u16(rnd_1_1.v[1], rnd_1_3.v[1]);
    vec256_t rnd_2_4; rnd_2_4.v[0] = vtrn1q_u16(rnd_1_4.v[0], rnd_1_6.v[0]); rnd_2_4.v[1] = vtrn1q_u16(rnd_1_4.v[1], rnd_1_6.v[1]);
    vec256_t rnd_2_6; rnd_2_6.v[0] = vtrn2q_u16(rnd_1_4.v[0], rnd_1_6.v[0]); rnd_2_6.v[1] = vtrn2q_u16(rnd_1_4.v[1], rnd_1_6.v[1]);
    vec256_t rnd_2_5; rnd_2_5.v[0] = vtrn1q_u16(rnd_1_5.v[0], rnd_1_7.v[0]); rnd_2_5.v[1] = vtrn1q_u16(rnd_1_5.v[1], rnd_1_7.v[1]);
    vec256_t rnd_2_7; rnd_2_7.v[0] = vtrn2q_u16(rnd_1_5.v[0], rnd_1_7.v[0]); rnd_2_7.v[1] = vtrn2q_u16(rnd_1_5.v[1], rnd_1_7.v[1]);
    vec256_t rnd_2_8 ; rnd_2_8.v[0]  = vtrn1q_u16(rnd_1_8.v[0], rnd_1_10.v[0]); rnd_2_8.v[1]  = vtrn1q_u16(rnd_1_8.v[1], rnd_1_10.v[1]);
    vec256_t rnd_2_10; rnd_2_10.v[0] = vtrn2q_u16(rnd_1_8.v[0], rnd_1_10.v[0]); rnd_2_10.v[1] = vtrn2q_u16(rnd_1_8.v[1], rnd_1_10.v[1]);
    vec256_t rnd_2_9 ; rnd_2_9.v[0]  = vtrn1q_u16(rnd_1_9.v[0], rnd_1_11.v[0]); rnd_2_9.v[1]  = vtrn1q_u16(rnd_1_9.v[1], rnd_1_11.v[1]);
    vec256_t rnd_2_11; rnd_2_11.v[0] = vtrn2q_u16(rnd_1_9.v[0], rnd_1_11.v[0]); rnd_2_11.v[1] = vtrn2q_u16(rnd_1_9.v[1], rnd_1_11.v[1]);
    vec256_t rnd_2_12; rnd_2_12.v[0] = vtrn1q_u16(rnd_1_12.v[0], rnd_1_14.v[0]); rnd_2_12.v[1] = vtrn1q_u16(rnd_1_12.v[1], rnd_1_14.v[1]);
    vec256_t rnd_2_14; rnd_2_14.v[0] = vtrn2q_u16(rnd_1_12.v[0], rnd_1_14.v[0]); rnd_2_14.v[1] = vtrn2q_u16(rnd_1_12.v[1], rnd_1_14.v[1]);
    vec256_t rnd_2_13; rnd_2_13.v[0] = vtrn1q_u16(rnd_1_13.v[0], rnd_1_15.v[0]); rnd_2_13.v[1] = vtrn1q_u16(rnd_1_13.v[1], rnd_1_15.v[1]);
    vec256_t rnd_2_15; rnd_2_15.v[0] = vtrn2q_u16(rnd_1_13.v[0], rnd_1_15.v[0]); rnd_2_15.v[1] = vtrn2q_u16(rnd_1_13.v[1], rnd_1_15.v[1]);
    vec256_t rnd_2_16; rnd_2_16.v[0] = vtrn1q_u16(rnd_1_16.v[0], rnd_1_18.v[0]); rnd_2_16.v[1] = vtrn1q_u16(rnd_1_16.v[1], rnd_1_18.v[1]);
    vec256_t rnd_2_18; rnd_2_18.v[0] = vtrn2q_u16(rnd_1_16.v[0], rnd_1_18.v[0]); rnd_2_18.v[1] = vtrn2q_u16(rnd_1_16.v[1], rnd_1_18.v[1]);
    vec256_t rnd_2_17; rnd_2_17.v[0] = vtrn1q_u16(rnd_1_17.v[0], rnd_1_19.v[0]); rnd_2_17.v[1] = vtrn1q_u16(rnd_1_17.v[1], rnd_1_19.v[1]);
    vec256_t rnd_2_19; rnd_2_19.v[0] = vtrn2q_u16(rnd_1_17.v[0], rnd_1_19.v[0]); rnd_2_19.v[1] = vtrn2q_u16(rnd_1_17.v[1], rnd_1_19.v[1]);
    vec256_t rnd_2_20; rnd_2_20.v[0] = vtrn1q_u16(rnd_1_20.v[0], rnd_1_22.v[0]); rnd_2_20.v[1] = vtrn1q_u16(rnd_1_20.v[1], rnd_1_22.v[1]);
    vec256_t rnd_2_22; rnd_2_22.v[0] = vtrn2q_u16(rnd_1_20.v[0], rnd_1_22.v[0]); rnd_2_22.v[1] = vtrn2q_u16(rnd_1_20.v[1], rnd_1_22.v[1]);
    vec256_t rnd_2_21; rnd_2_21.v[0] = vtrn1q_u16(rnd_1_21.v[0], rnd_1_23.v[0]); rnd_2_21.v[1] = vtrn1q_u16(rnd_1_21.v[1], rnd_1_23.v[1]);
    vec256_t rnd_2_23; rnd_2_23.v[0] = vtrn2q_u16(rnd_1_21.v[0], rnd_1_23.v[0]); rnd_2_23.v[1] = vtrn2q_u16(rnd_1_21.v[1], rnd_1_23.v[1]);
    vec256_t rnd_2_24; rnd_2_24.v[0] = vtrn1q_u16(rnd_1_24.v[0], rnd_1_26.v[0]); rnd_2_24.v[1] = vtrn1q_u16(rnd_1_24.v[1], rnd_1_26.v[1]);
    vec256_t rnd_2_26; rnd_2_26.v[0] = vtrn2q_u16(rnd_1_24.v[0], rnd_1_26.v[0]); rnd_2_26.v[1] = vtrn2q_u16(rnd_1_24.v[1], rnd_1_26.v[1]);
    vec256_t rnd_2_25; rnd_2_25.v[0] = vtrn1q_u16(rnd_1_25.v[0], rnd_1_27.v[0]); rnd_2_25.v[1] = vtrn1q_u16(rnd_1_25.v[1], rnd_1_27.v[1]);
    vec256_t rnd_2_27; rnd_2_27.v[0] = vtrn2q_u16(rnd_1_25.v[0], rnd_1_27.v[0]); rnd_2_27.v[1] = vtrn2q_u16(rnd_1_25.v[1], rnd_1_27.v[1]);
    vec256_t rnd_2_28; rnd_2_28.v[0] = vtrn1q_u16(rnd_1_28.v[0], rnd_1_30.v[0]); rnd_2_28.v[1] = vtrn1q_u16(rnd_1_28.v[1], rnd_1_30.v[1]);
    vec256_t rnd_2_30; rnd_2_30.v[0] = vtrn2q_u16(rnd_1_28.v[0], rnd_1_30.v[0]); rnd_2_30.v[1] = vtrn2q_u16(rnd_1_28.v[1], rnd_1_30.v[1]);
    vec256_t rnd_2_29; rnd_2_29.v[0] = vtrn1q_u16(rnd_1_29.v[0], rnd_1_31.v[0]); rnd_2_29.v[1] = vtrn1q_u16(rnd_1_29.v[1], rnd_1_31.v[1]);
    vec256_t rnd_2_31; rnd_2_31.v[0] = vtrn2q_u16(rnd_1_29.v[0], rnd_1_31.v[0]); rnd_2_31.v[1] = vtrn2q_u16(rnd_1_29.v[1], rnd_1_31.v[1]);

    vec256_t rnd_3_0; rnd_3_0.v[0] = vtrn1q_u32(rnd_2_0.v[0], rnd_2_4.v[0]); rnd_3_0.v[1] = vtrn1q_u32(rnd_2_0.v[1], rnd_2_4.v[1]);
    vec256_t rnd_3_4; rnd_3_4.v[0] = vtrn2q_u32(rnd_2_0.v[0], rnd_2_4.v[0]); rnd_3_4.v[1] = vtrn2q_u32(rnd_2_0.v[1], rnd_2_4.v[1]);
    vec256_t rnd_3_1; rnd_3_1.v[0] = vtrn1q_u32(rnd_2_1.v[0], rnd_2_5.v[0]); rnd_3_1.v[1] = vtrn1q_u32(rnd_2_1.v[1], rnd_2_5.v[1]);
    vec256_t rnd_3_5; rnd_3_5.v[0] = vtrn2q_u32(rnd_2_1.v[0], rnd_2_5.v[0]); rnd_3_5.v[1] = vtrn2q_u32(rnd_2_1.v[1], rnd_2_5.v[1]);
    vec256_t rnd_3_2; rnd_3_2.v[0] = vtrn1q_u32(rnd_2_2.v[0], rnd_2_6.v[0]); rnd_3_2.v[1] = vtrn1q_u32(rnd_2_2.v[1], rnd_2_6.v[1]);
    vec256_t rnd_3_6; rnd_3_6.v[0] = vtrn2q_u32(rnd_2_2.v[0], rnd_2_6.v[0]); rnd_3_6.v[1] = vtrn2q_u32(rnd_2_2.v[1], rnd_2_6.v[1]);
    vec256_t rnd_3_3; rnd_3_3.v[0] = vtrn1q_u32(rnd_2_3.v[0], rnd_2_7.v[0]); rnd_3_3.v[1] = vtrn1q_u32(rnd_2_3.v[1], rnd_2_7.v[1]);
    vec256_t rnd_3_7; rnd_3_7.v[0] = vtrn2q_u32(rnd_2_3.v[0], rnd_2_7.v[0]); rnd_3_7.v[1] = vtrn2q_u32(rnd_2_3.v[1], rnd_2_7.v[1]);
    vec256_t rnd_3_8 ; rnd_3_8.v[0]  = vtrn1q_u32(rnd_2_8.v[0], rnd_2_12.v[0]); rnd_3_8.v[1]  = vtrn1q_u32(rnd_2_8.v[1], rnd_2_12.v[1]);
    vec256_t rnd_3_12; rnd_3_12.v[0] = vtrn2q_u32(rnd_2_8.v[0], rnd_2_12.v[0]); rnd_3_12.v[1] = vtrn2q_u32(rnd_2_8.v[1], rnd_2_12.v[1]);
    vec256_t rnd_3_9;  rnd_3_9.v[0] =  vtrn1q_u32(rnd_2_9.v[0], rnd_2_13.v[0]); rnd_3_9.v[1] =  vtrn1q_u32(rnd_2_9.v[1], rnd_2_13.v[1]);
    vec256_t rnd_3_13; rnd_3_13.v[0] = vtrn2q_u32(rnd_2_9.v[0], rnd_2_13.v[0]); rnd_3_13.v[1] = vtrn2q_u32(rnd_2_9.v[1], rnd_2_13.v[1]);
    vec256_t rnd_3_10; rnd_3_10.v[0] = vtrn1q_u32(rnd_2_10.v[0], rnd_2_14.v[0]); rnd_3_10.v[1] = vtrn1q_u32(rnd_2_10.v[1], rnd_2_14.v[1]);
    vec256_t rnd_3_14; rnd_3_14.v[0] = vtrn2q_u32(rnd_2_10.v[0], rnd_2_14.v[0]); rnd_3_14.v[1] = vtrn2q_u32(rnd_2_10.v[1], rnd_2_14.v[1]);
    vec256_t rnd_3_11; rnd_3_11.v[0] = vtrn1q_u32(rnd_2_11.v[0], rnd_2_15.v[0]); rnd_3_11.v[1] = vtrn1q_u32(rnd_2_11.v[1], rnd_2_15.v[1]);
    vec256_t rnd_3_15; rnd_3_15.v[0] = vtrn2q_u32(rnd_2_11.v[0], rnd_2_15.v[0]); rnd_3_15.v[1] = vtrn2q_u32(rnd_2_11.v[1], rnd_2_15.v[1]);
    vec256_t rnd_3_16; rnd_3_16.v[0] = vtrn1q_u32(rnd_2_16.v[0], rnd_2_20.v[0]); rnd_3_16.v[1] = vtrn1q_u32(rnd_2_16.v[1], rnd_2_20.v[1]);
    vec256_t rnd_3_20; rnd_3_20.v[0] = vtrn2q_u32(rnd_2_16.v[0], rnd_2_20.v[0]); rnd_3_20.v[1] = vtrn2q_u32(rnd_2_16.v[1], rnd_2_20.v[1]);
    vec256_t rnd_3_17; rnd_3_17.v[0] = vtrn1q_u32(rnd_2_17.v[0], rnd_2_21.v[0]); rnd_3_17.v[1] = vtrn1q_u32(rnd_2_17.v[1], rnd_2_21.v[1]);
    vec256_t rnd_3_21; rnd_3_21.v[0] = vtrn2q_u32(rnd_2_17.v[0], rnd_2_21.v[0]); rnd_3_21.v[1] = vtrn2q_u32(rnd_2_17.v[1], rnd_2_21.v[1]);
    vec256_t rnd_3_18; rnd_3_18.v[0] = vtrn1q_u32(rnd_2_18.v[0], rnd_2_22.v[0]); rnd_3_18.v[1] = vtrn1q_u32(rnd_2_18.v[1], rnd_2_22.v[1]);
    vec256_t rnd_3_22; rnd_3_22.v[0] = vtrn2q_u32(rnd_2_18.v[0], rnd_2_22.v[0]); rnd_3_22.v[1] = vtrn2q_u32(rnd_2_18.v[1], rnd_2_22.v[1]);
    vec256_t rnd_3_19; rnd_3_19.v[0] = vtrn1q_u32(rnd_2_19.v[0], rnd_2_23.v[0]); rnd_3_19.v[1] = vtrn1q_u32(rnd_2_19.v[1], rnd_2_23.v[1]);
    vec256_t rnd_3_23; rnd_3_23.v[0] = vtrn2q_u32(rnd_2_19.v[0], rnd_2_23.v[0]); rnd_3_23.v[1] = vtrn2q_u32(rnd_2_19.v[1], rnd_2_23.v[1]);
    vec256_t rnd_3_24; rnd_3_24.v[0] = vtrn1q_u32(rnd_2_24.v[0], rnd_2_28.v[0]); rnd_3_24.v[1] = vtrn1q_u32(rnd_2_24.v[1], rnd_2_28.v[1]);
    vec256_t rnd_3_28; rnd_3_28.v[0] = vtrn2q_u32(rnd_2_24.v[0], rnd_2_28.v[0]); rnd_3_28.v[1] = vtrn2q_u32(rnd_2_24.v[1], rnd_2_28.v[1]);
    vec256_t rnd_3_25; rnd_3_25.v[0] = vtrn1q_u32(rnd_2_25.v[0], rnd_2_29.v[0]); rnd_3_25.v[1] = vtrn1q_u32(rnd_2_25.v[1], rnd_2_29.v[1]);
    vec256_t rnd_3_29; rnd_3_29.v[0] = vtrn2q_u32(rnd_2_25.v[0], rnd_2_29.v[0]); rnd_3_29.v[1] = vtrn2q_u32(rnd_2_25.v[1], rnd_2_29.v[1]);
    vec256_t rnd_3_26; rnd_3_26.v[0] = vtrn1q_u32(rnd_2_26.v[0], rnd_2_30.v[0]); rnd_3_26.v[1] = vtrn1q_u32(rnd_2_26.v[1], rnd_2_30.v[1]);
    vec256_t rnd_3_30; rnd_3_30.v[0] = vtrn2q_u32(rnd_2_26.v[0], rnd_2_30.v[0]); rnd_3_30.v[1] = vtrn2q_u32(rnd_2_26.v[1], rnd_2_30.v[1]);
    vec256_t rnd_3_27; rnd_3_27.v[0] = vtrn1q_u32(rnd_2_27.v[0], rnd_2_31.v[0]); rnd_3_27.v[1] = vtrn1q_u32(rnd_2_27.v[1], rnd_2_31.v[1]);
    vec256_t rnd_3_31; rnd_3_31.v[0] = vtrn2q_u32(rnd_2_27.v[0], rnd_2_31.v[0]); rnd_3_31.v[1] = vtrn2q_u32(rnd_2_27.v[1], rnd_2_31.v[1]);

    vec256_t rnd_4_0; rnd_4_0.v[0] = vtrn1q_u64(rnd_3_0.v[0], rnd_3_8.v[0]); rnd_4_0.v[1] = vtrn1q_u64(rnd_3_0.v[1], rnd_3_8.v[1]);
    vec256_t rnd_4_8; rnd_4_8.v[0] = vtrn2q_u64(rnd_3_0.v[0], rnd_3_8.v[0]); rnd_4_8.v[1] = vtrn2q_u64(rnd_3_0.v[1], rnd_3_8.v[1]);
    vec256_t rnd_4_1; rnd_4_1.v[0] = vtrn1q_u64(rnd_3_1.v[0], rnd_3_9.v[0]); rnd_4_1.v[1] = vtrn1q_u64(rnd_3_1.v[1], rnd_3_9.v[1]);
    vec256_t rnd_4_9; rnd_4_9.v[0] = vtrn2q_u64(rnd_3_1.v[0], rnd_3_9.v[0]); rnd_4_9.v[1] = vtrn2q_u64(rnd_3_1.v[1], rnd_3_9.v[1]);
    vec256_t rnd_4_2 ; rnd_4_2.v[0]  = vtrn1q_u64(rnd_3_2.v[0], rnd_3_10.v[0]); rnd_4_2.v[1]  = vtrn1q_u64(rnd_3_2.v[1], rnd_3_10.v[1]);
    vec256_t rnd_4_10; rnd_4_10.v[0] = vtrn2q_u64(rnd_3_2.v[0], rnd_3_10.v[0]); rnd_4_10.v[1] = vtrn2q_u64(rnd_3_2.v[1], rnd_3_10.v[1]);
    vec256_t rnd_4_3;  rnd_4_3.v[0]  = vtrn1q_u64(rnd_3_3.v[0], rnd_3_11.v[0]); rnd_4_3.v[1]  = vtrn1q_u64(rnd_3_3.v[1], rnd_3_11.v[1]);
    vec256_t rnd_4_11; rnd_4_11.v[0] = vtrn2q_u64(rnd_3_3.v[0], rnd_3_11.v[0]); rnd_4_11.v[1] = vtrn2q_u64(rnd_3_3.v[1], rnd_3_11.v[1]);
    vec256_t rnd_4_4;  rnd_4_4.v[0]  = vtrn1q_u64(rnd_3_4.v[0], rnd_3_12.v[0]); rnd_4_4.v[1]  = vtrn1q_u64(rnd_3_4.v[1], rnd_3_12.v[1]);
    vec256_t rnd_4_12; rnd_4_12.v[0] = vtrn2q_u64(rnd_3_4.v[0], rnd_3_12.v[0]); rnd_4_12.v[1] = vtrn2q_u64(rnd_3_4.v[1], rnd_3_12.v[1]);
    vec256_t rnd_4_5;  rnd_4_5.v[0]  = vtrn1q_u64(rnd_3_5.v[0], rnd_3_13.v[0]); rnd_4_5.v[1]  = vtrn1q_u64(rnd_3_5.v[1], rnd_3_13.v[1]);
    vec256_t rnd_4_13; rnd_4_13.v[0] = vtrn2q_u64(rnd_3_5.v[0], rnd_3_13.v[0]); rnd_4_13.v[1] = vtrn2q_u64(rnd_3_5.v[1], rnd_3_13.v[1]);
    vec256_t rnd_4_6;  rnd_4_6.v[0]  = vtrn1q_u64(rnd_3_6.v[0], rnd_3_14.v[0]); rnd_4_6.v[1]  = vtrn1q_u64(rnd_3_6.v[1], rnd_3_14.v[1]);
    vec256_t rnd_4_14; rnd_4_14.v[0] = vtrn2q_u64(rnd_3_6.v[0], rnd_3_14.v[0]); rnd_4_14.v[1] = vtrn2q_u64(rnd_3_6.v[1], rnd_3_14.v[1]);
    vec256_t rnd_4_7;  rnd_4_7.v[0]  = vtrn1q_u64(rnd_3_7.v[0], rnd_3_15.v[0]); rnd_4_7.v[1]  = vtrn1q_u64(rnd_3_7.v[1], rnd_3_15.v[1]);
    vec256_t rnd_4_15; rnd_4_15.v[0] = vtrn2q_u64(rnd_3_7.v[0], rnd_3_15.v[0]); rnd_4_15.v[1] = vtrn2q_u64(rnd_3_7.v[1], rnd_3_15.v[1]);
    vec256_t rnd_4_16; rnd_4_16.v[0] = vtrn1q_u64(rnd_3_16.v[0], rnd_3_24.v[0]); rnd_4_16.v[1] = vtrn1q_u64(rnd_3_16.v[1], rnd_3_24.v[1]);
    vec256_t rnd_4_24; rnd_4_24.v[0] = vtrn2q_u64(rnd_3_16.v[0], rnd_3_24.v[0]); rnd_4_24.v[1] = vtrn2q_u64(rnd_3_16.v[1], rnd_3_24.v[1]);
    vec256_t rnd_4_17; rnd_4_17.v[0] = vtrn1q_u64(rnd_3_17.v[0], rnd_3_25.v[0]); rnd_4_17.v[1] = vtrn1q_u64(rnd_3_17.v[1], rnd_3_25.v[1]);
    vec256_t rnd_4_25; rnd_4_25.v[0] = vtrn2q_u64(rnd_3_17.v[0], rnd_3_25.v[0]); rnd_4_25.v[1] = vtrn2q_u64(rnd_3_17.v[1], rnd_3_25.v[1]);
    vec256_t rnd_4_18; rnd_4_18.v[0] = vtrn1q_u64(rnd_3_18.v[0], rnd_3_26.v[0]); rnd_4_18.v[1] = vtrn1q_u64(rnd_3_18.v[1], rnd_3_26.v[1]);
    vec256_t rnd_4_26; rnd_4_26.v[0] = vtrn2q_u64(rnd_3_18.v[0], rnd_3_26.v[0]); rnd_4_26.v[1] = vtrn2q_u64(rnd_3_18.v[1], rnd_3_26.v[1]);
    vec256_t rnd_4_19; rnd_4_19.v[0] = vtrn1q_u64(rnd_3_19.v[0], rnd_3_27.v[0]); rnd_4_19.v[1] = vtrn1q_u64(rnd_3_19.v[1], rnd_3_27.v[1]);
    vec256_t rnd_4_27; rnd_4_27.v[0] = vtrn2q_u64(rnd_3_19.v[0], rnd_3_27.v[0]); rnd_4_27.v[1] = vtrn2q_u64(rnd_3_19.v[1], rnd_3_27.v[1]);
    vec256_t rnd_4_20; rnd_4_20.v[0] = vtrn1q_u64(rnd_3_20.v[0], rnd_3_28.v[0]); rnd_4_20.v[1] = vtrn1q_u64(rnd_3_20.v[1], rnd_3_28.v[1]);
    vec256_t rnd_4_28; rnd_4_28.v[0] = vtrn2q_u64(rnd_3_20.v[0], rnd_3_28.v[0]); rnd_4_28.v[1] = vtrn2q_u64(rnd_3_20.v[1], rnd_3_28.v[1]);
    vec256_t rnd_4_21; rnd_4_21.v[0] = vtrn1q_u64(rnd_3_21.v[0], rnd_3_29.v[0]); rnd_4_21.v[1] = vtrn1q_u64(rnd_3_21.v[1], rnd_3_29.v[1]);
    vec256_t rnd_4_29; rnd_4_29.v[0] = vtrn2q_u64(rnd_3_21.v[0], rnd_3_29.v[0]); rnd_4_29.v[1] = vtrn2q_u64(rnd_3_21.v[1], rnd_3_29.v[1]);
    vec256_t rnd_4_22; rnd_4_22.v[0] = vtrn1q_u64(rnd_3_22.v[0], rnd_3_30.v[0]); rnd_4_22.v[1] = vtrn1q_u64(rnd_3_22.v[1], rnd_3_30.v[1]);
    vec256_t rnd_4_30; rnd_4_30.v[0] = vtrn2q_u64(rnd_3_22.v[0], rnd_3_30.v[0]); rnd_4_30.v[1] = vtrn2q_u64(rnd_3_22.v[1], rnd_3_30.v[1]);
    vec256_t rnd_4_23; rnd_4_23.v[0] = vtrn1q_u64(rnd_3_23.v[0], rnd_3_31.v[0]); rnd_4_23.v[1] = vtrn1q_u64(rnd_3_23.v[1], rnd_3_31.v[1]);
    vec256_t rnd_4_31; rnd_4_31.v[0] = vtrn2q_u64(rnd_3_23.v[0], rnd_3_31.v[0]); rnd_4_31.v[1] = vtrn2q_u64(rnd_3_23.v[1], rnd_3_31.v[1]);

    // NOTE: maybe It's useful to reoptimize the following code. As we are just shuffling some register we could do
    // already do this in round 4, just above.
    vec256_t rnd_5_0;  rnd_5_0.v[0]  = rnd_4_0.v[0]; rnd_5_0.v[1]  = rnd_4_16.v[0];
    vec256_t rnd_5_16; rnd_5_16.v[0] = rnd_4_0.v[1]; rnd_5_16.v[1] = rnd_4_16.v[1];
    *((vec256_t *)(dst_origin +  0*dst_stride)) = rnd_5_0;
    *((vec256_t *)(dst_origin + 16*dst_stride)) = rnd_5_16;
    vec256_t rnd_5_1;  rnd_5_1.v[0]  = rnd_4_1.v[0]; rnd_5_1.v[1]  = rnd_4_17.v[0];
    vec256_t rnd_5_17; rnd_5_17.v[0] = rnd_4_1.v[1]; rnd_5_17.v[1] = rnd_4_17.v[1];
    *((vec256_t *)(dst_origin +  1*dst_stride)) = rnd_5_1;
    *((vec256_t *)(dst_origin + 17*dst_stride)) = rnd_5_17;
    vec256_t rnd_5_2;  rnd_5_2.v[0]  = rnd_4_2.v[0];  rnd_5_2.v[1] = rnd_4_18.v[0];
    vec256_t rnd_5_18; rnd_5_18.v[0] = rnd_4_2.v[1]; rnd_5_18.v[1] = rnd_4_18.v[1];
    *((vec256_t *)(dst_origin +  2*dst_stride)) = rnd_5_2;
    *((vec256_t *)(dst_origin + 18*dst_stride)) = rnd_5_18;
    vec256_t rnd_5_3;  rnd_5_3.v[0]  = rnd_4_3.v[0];  rnd_5_3.v[1] = rnd_4_19.v[0];
    vec256_t rnd_5_19; rnd_5_19.v[0] = rnd_4_3.v[1]; rnd_5_19.v[1] = rnd_4_19.v[1];
    *((vec256_t *)(dst_origin +  3*dst_stride)) = rnd_5_3;
    *((vec256_t *)(dst_origin + 19*dst_stride)) = rnd_5_19;
    vec256_t rnd_5_4;  rnd_5_4.v[0]  = rnd_4_4.v[0];  rnd_5_4.v[1] = rnd_4_20.v[0];
    vec256_t rnd_5_20; rnd_5_20.v[0] = rnd_4_4.v[1]; rnd_5_20.v[1] = rnd_4_20.v[1];
    *((vec256_t *)(dst_origin +  4*dst_stride)) = rnd_5_4;
    *((vec256_t *)(dst_origin + 20*dst_stride)) = rnd_5_20;
    vec256_t rnd_5_5;  rnd_5_5.v[0]  = rnd_4_5.v[0];  rnd_5_5.v[1] = rnd_4_21.v[0];
    vec256_t rnd_5_21; rnd_5_21.v[0] = rnd_4_5.v[1]; rnd_5_21.v[1] = rnd_4_21.v[1];
    *((vec256_t *)(dst_origin +  5*dst_stride)) = rnd_5_5;
    *((vec256_t *)(dst_origin + 21*dst_stride)) = rnd_5_21;
    vec256_t rnd_5_6;  rnd_5_6.v[0]  = rnd_4_6.v[0];  rnd_5_6.v[1] = rnd_4_22.v[0];
    vec256_t rnd_5_22; rnd_5_22.v[0] = rnd_4_6.v[1]; rnd_5_22.v[1] = rnd_4_22.v[1];
    *((vec256_t *)(dst_origin +  6*dst_stride)) = rnd_5_6;
    *((vec256_t *)(dst_origin + 22*dst_stride)) = rnd_5_22;
    vec256_t rnd_5_7;  rnd_5_7.v[0]  = rnd_4_7.v[0];  rnd_5_7.v[1] = rnd_4_23.v[0];
    vec256_t rnd_5_23; rnd_5_23.v[0] = rnd_4_7.v[1]; rnd_5_23.v[1] = rnd_4_23.v[1];
    *((vec256_t *)(dst_origin +  7*dst_stride)) = rnd_5_7;
    *((vec256_t *)(dst_origin + 23*dst_stride)) = rnd_5_23;
    vec256_t rnd_5_8;  rnd_5_8.v[0]  = rnd_4_8.v[0];  rnd_5_8.v[1] = rnd_4_24.v[0];
    vec256_t rnd_5_24; rnd_5_24.v[0] = rnd_4_8.v[1]; rnd_5_24.v[1] = rnd_4_24.v[1];
    *((vec256_t *)(dst_origin +  8*dst_stride)) = rnd_5_8;
    *((vec256_t *)(dst_origin + 24*dst_stride)) = rnd_5_24;
    vec256_t rnd_5_9;  rnd_5_9.v[0]  = rnd_4_9.v[0];  rnd_5_9.v[1] = rnd_4_25.v[0];
    vec256_t rnd_5_25; rnd_5_25.v[0] = rnd_4_9.v[1]; rnd_5_25.v[1] = rnd_4_25.v[1];
    *((vec256_t *)(dst_origin +  9*dst_stride)) = rnd_5_9;
    *((vec256_t *)(dst_origin + 25*dst_stride)) = rnd_5_25;
    vec256_t rnd_5_10; rnd_5_10.v[0] = rnd_4_10.v[0]; rnd_5_10.v[1] = rnd_4_26.v[0];
    vec256_t rnd_5_26; rnd_5_26.v[0] = rnd_4_10.v[1]; rnd_5_26.v[1] = rnd_4_26.v[1];
    *((vec256_t *)(dst_origin + 10*dst_stride)) = rnd_5_10;
    *((vec256_t *)(dst_origin + 26*dst_stride)) = rnd_5_26;
    vec256_t rnd_5_11; rnd_5_11.v[0] = rnd_4_11.v[0]; rnd_5_11.v[1] = rnd_4_27.v[0];
    vec256_t rnd_5_27; rnd_5_27.v[0] = rnd_4_11.v[1]; rnd_5_27.v[1] = rnd_4_27.v[1];
    *((vec256_t *)(dst_origin + 11*dst_stride)) = rnd_5_11;
    *((vec256_t *)(dst_origin + 27*dst_stride)) = rnd_5_27;
    vec256_t rnd_5_12; rnd_5_12.v[0] = rnd_4_12.v[0]; rnd_5_12.v[1] = rnd_4_28.v[0];
    vec256_t rnd_5_28; rnd_5_28.v[0] = rnd_4_12.v[1]; rnd_5_28.v[1] = rnd_4_28.v[1];
    *((vec256_t *)(dst_origin + 12*dst_stride)) = rnd_5_12;
    *((vec256_t *)(dst_origin + 28*dst_stride)) = rnd_5_28;
    vec256_t rnd_5_13; rnd_5_13.v[0] = rnd_4_13.v[0]; rnd_5_13.v[1] = rnd_4_29.v[0];
    vec256_t rnd_5_29; rnd_5_29.v[0] = rnd_4_13.v[1]; rnd_5_29.v[1] = rnd_4_29.v[1];
    *((vec256_t *)(dst_origin + 13*dst_stride)) = rnd_5_13;
    *((vec256_t *)(dst_origin + 29*dst_stride)) = rnd_5_29;
    vec256_t rnd_5_14; rnd_5_14.v[0] = rnd_4_14.v[0]; rnd_5_14.v[1] = rnd_4_30.v[0];
    vec256_t rnd_5_30; rnd_5_30.v[0] = rnd_4_14.v[1]; rnd_5_30.v[1] = rnd_4_30.v[1];
    *((vec256_t *)(dst_origin + 14*dst_stride)) = rnd_5_14;
    *((vec256_t *)(dst_origin + 30*dst_stride)) = rnd_5_30;
    vec256_t rnd_5_15; rnd_5_15.v[0] = rnd_4_15.v[0]; rnd_5_15.v[1] = rnd_4_31.v[0];
    vec256_t rnd_5_31; rnd_5_31.v[0] = rnd_4_15.v[1]; rnd_5_31.v[1] = rnd_4_31.v[1];
    *((vec256_t *)(dst_origin + 15*dst_stride)) = rnd_5_15;
    *((vec256_t *)(dst_origin + 31*dst_stride)) = rnd_5_31;
}
