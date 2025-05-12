/**
 *
 * Optimized Implementation of LESS.
 *
 * @version 1.2 (May 2023)
 *
 * @author Duc Tri Nguyen <dnguye69@gmu.edu>
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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

#include "parameters.h"
#include "rng.h"

// number of needed to represent q = 7
#define NUM_BITS_Q (BITS_TO_REPRESENT(Q))


#define DEF_RAND_STATE(FUNC_NAME, EL_T, MINV, MAXV) \
static inline void FUNC_NAME(SHAKE_STATE_STRUCT *shake_monomial_state, EL_T *buffer, size_t num_elements) { \
   typedef uint64_t WORD_T; \
   static const EL_T MIN_VALUE = (MINV);\
   static const EL_T MAX_VALUE = (MAXV); \
   static const EL_T SPAN = MAX_VALUE - MIN_VALUE; \
   static const size_t REQ_BITS = BITS_TO_REPRESENT(SPAN); \
   static const EL_T EL_MASK = ((EL_T) 1 << REQ_BITS) - 1; \
   WORD_T word; \
   size_t count = 0; \
   do { \
      csprng_randombytes((unsigned char *) &word, sizeof(WORD_T), shake_monomial_state); \
      for (uint32_t i = 0; i < ((sizeof(WORD_T)*8u) / REQ_BITS); i++) { \
         EL_T rnd_value = word & EL_MASK; \
         if (rnd_value <= SPAN) buffer[count++] = rnd_value + MIN_VALUE; \
         if (count >= num_elements) return; \
         word >>= REQ_BITS; \
      } \
   } while (1); }

#define DEF_RAND(FUNC_NAME, EL_T, MINV, MAXV) \
static inline void FUNC_NAME(EL_T *buffer, size_t num_elements) { \
   typedef uint64_t WORD_T; \
   static const EL_T MIN_VALUE = (MINV); \
   static const EL_T MAX_VALUE = (MAXV); \
   static const EL_T SPAN = MAX_VALUE - MIN_VALUE; \
   static const size_t REQ_BITS = BITS_TO_REPRESENT(SPAN); \
   static const EL_T EL_MASK = ((EL_T) 1 << REQ_BITS) - 1; \
   WORD_T word; \
   size_t count = 0; \
   do { \
      randombytes((unsigned char *) &word, sizeof(WORD_T)); \
      for (uint32_t i = 0; i < ((sizeof(WORD_T)*8) / REQ_BITS); i++) { \
         EL_T rnd_value = word & EL_MASK; \
         if (rnd_value <= SPAN) buffer[count++] = rnd_value + MIN_VALUE; \
         if (count >= num_elements) return; \
         word >>= REQ_BITS; \
      } \
   } while (1); }



/// constant time implementation
/// \param x[in]: input number < 2*127
/// \return x mod 127
static inline
FQ_ELEM fq_cond_sub(const FQ_ELEM x) {
    // equivalent to: (x >= Q) ? (x - Q) : x
    // likely to be ~ constant-time (a "smart" compiler might turn this into conditionals though)
    FQ_ELEM sub_q = x - Q;
    FQ_ELEM mask = -(sub_q >> NUM_BITS_Q);
    return (mask & Q) + sub_q;
} /* end fq_cond_sub */

/// constant time implementation
/// \param x[in]: input number < 127**2
/// \return x mod 127
static inline
FQ_ELEM fq_red(const FQ_DOUBLEPREC x) {
    return fq_cond_sub((x >> NUM_BITS_Q) + ((FQ_ELEM) x & Q));
} /* end fq_red */

/// constant time implementation
/// \param x[in]: minuend < 127
/// \param y[in]: subtrahend < 127
/// \return difference x-y mod 127
static inline
FQ_ELEM fq_sub(const FQ_ELEM x,
			   const FQ_ELEM y) {
    return fq_cond_sub(x + Q - y);
} /* end fq_sub */

/// constant time implementation
/// \param x[in]: first summand < 127
/// \param y[in]: second summand < 127
/// \return sum x+y mod 127
static inline
FQ_ELEM fq_mul(const FQ_ELEM x,
			   const FQ_ELEM y) {
    return fq_red(((FQ_DOUBLEPREC)x) *(FQ_DOUBLEPREC)y);
} /* end fq_mul */

/// constant time implementation
/// \param x[in]: first factor < 127
/// \param y[in]: second factor < 127
/// \return product x*y mod 127
static inline
FQ_ELEM fq_add(const FQ_ELEM x,
			   const FQ_ELEM y) {
      return (x + y) % Q;
} /* end fq_add */

/// NOTE: maybe dont use it for sensitive data
static const uint8_t fq_inv_table[128] __attribute__((aligned(64))) = {
   0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0
};

/// NOTE: non constant-time. Don't use for anything important.
/// NOTE: input must be reduced.
/// \param x[in]: input value < 127
/// \return x^{-1} mod 127
static inline
FQ_ELEM fq_inv(FQ_ELEM x) {
   return fq_inv_table[x];
} /* end fq_inv */

/// Sampling functions from the global TRNG state
DEF_RAND(fq_star_rnd_elements, FQ_ELEM, 1, Q-1)
DEF_RAND(rand_range_q_elements, FQ_ELEM, 0, Q-1)

/// Sampling functions from the taking the PRNG state as a parameter
DEF_RAND_STATE(fq_star_rnd_state_elements, FQ_ELEM, 1, Q-1)
DEF_RAND_STATE(rand_range_q_state_elements, FQ_ELEM, 0, Q-1)

// load neon macros
#include "macro.h"

  
/// \param a[in]: [a_0 mod 127, ..., a_15 mod 127]
/// \param b[in]: [b_0 mod 127, ..., b_15 mod 127]
/// \return [a_0 * b_0 mod 127, ..., a_15 * b_15 mod 127]
static inline
uint8x16_t gf127v_mul_u128(const uint8x16_t a,
                           const uint8x16_t b) {
    union {
        uint8x16_t  x1;
        uint8x8x2_t x2;
    } tmp;
    const uint8x16_t q  = vdupq_n_u8(0x7f);
    const uint16x8_t q2 = vdupq_n_u16(0x007f);

    const uint8x8_t la = vget_low_u8(a);
    const uint8x8_t lb = vget_low_u8(b);
    const uint8x8_t ha = vget_high_u8(a);
    const uint8x8_t hb = vget_high_u8(b);

    const uint16x8_t lc = vmull_u8(la, lb);
    const uint16x8_t hc = vmull_u8(ha, hb);

    const uint16x8_t lt1 = vshrq_n_u16(lc, 7);
    const uint16x8_t ht1 = vshrq_n_u16(hc, 7);
    const uint16x8_t lt2 = vaddq_u16(lt1, vandq_u16(lc, q2));
    const uint16x8_t ht2 = vaddq_u16(ht1, vandq_u16(hc, q2));
    const uint8x8_t tl = vmovn_u16(lt2);
    const uint8x8_t th = vmovn_u16(ht2);
    tmp.x2.val[0] = tl;
    tmp.x2.val[1] = th;
    const uint8x16_t c = tmp.x1;

    const uint8x16_t t   = vsubq_u8(c, q);
    const uint8x16_t m   = vshrq_n_s8(c, 7);
    const uint8x16_t r   = vbslq_u8(m, t, c);
    return r;
} /* end gf127v_mul_u128 */

/// \param ret[out]: pointer to two uint8x16x4_t which will contain the table
/// \param a[in]: scalar to which the table should be loaded
static inline void gf127v_scalar_compute_table(uint8x16x4_t *ret,
											   const uint8_t a) {
	ret[0] = vld1q_u8_x4(__gf127_lookuptable + a*128 +  0);
	ret[1] = vld1q_u8_x4(__gf127_lookuptable + a*128 + 64);
} /* end gf127v_scalar_compute_table */

/// \param b[in]: vector register container 16 Fq elements: [b_0, ..., b_15]
/// \param table[in]: output of `gf127v_scalar_compute_table` 
/// \return [b_0 * a, ..., b_15 * a], where `a` is the input to `gf127v_scalar_compute_table`
static inline uint8x16_t gf127v_scalar_table(const uint8x16_t b,
											 const uint8x16x4_t *table) {
	const uint8x16_t zero = vdupq_n_u8(0);
	const uint8x16_t mask = vdupq_n_u8(64);
	const uint8x16_t c = vsubq_u8(b, mask);
	const uint8x16_t ret1 = vqtbx4q_u8(zero, table[0], b);
	const uint8x16_t ret2 = vqtbx4q_u8(ret1, table[1], c);
	return ret2;
} /* end gf127v_scalar_table */

/// \param a[in]: [a_0 mod 127, ..., a_31 mod 127]
/// \param b[in]: [b_0 mod 127, ..., b_31 mod 127]
/// \return [a_0 * b_0 mod 127, ..., a_31 * b_31 mod 127]
static inline
vec256_t gf127v_mul_u256(const vec256_t a,
						 const vec256_t b) {
	vec256_t ret;
	ret.v[0] = gf127v_mul_u128(a.v[0], b.v[0]);
	ret.v[1] = gf127v_mul_u128(a.v[1], b.v[1]);
	return ret;
} /* end gf127v_mul_u256 */

/// \param b[in]: vector register container 32 Fq elements: [b_0, ..., b_31]
/// \param table[in]: output of `gf127v_scalar_compute_table` 
/// \return [b_0 * a, ..., b_31 * a], where `a` is the input to `gf127v_scalar_compute_table`
static inline vec256_t gf127v_scalar_table_full(const vec256_t b,
												const uint8x16x4_t *table) {
	vec256_t ret;
	ret.v[0] = gf127v_scalar_table(b.v[0], table);
	ret.v[1] = gf127v_scalar_table(b.v[1], table);
	return ret;
} /* end gf127v_scalar_table_full */

/// \param t[in]: two neon registers: [t_0, ..., t_31]
/// \return t_0 + t_1 + ... + t_3§ mod 127
static inline uint8_t vhadd8(const vec256_t t) {
	vec128_t v;
	v.v = vaddq_u8(t.v[0], t.v[1]);
	uint8_t r = fq_red(v.v8[0]);
	for (uint32_t i = 1; i < 16; i++) {
		r = fq_add(r, fq_red(v.v8[i]));
	}

	return r;
} /* end vhadd8 */

/// NOTE: these functions are outsourced to this file, to make the
/// optimizied implementation as easy as possible.
/// accumulates a row
/// \param row[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \return sum(d[i]) for i in range(N-K)
static inline
FQ_ELEM row_acc(const FQ_ELEM *row) {
    vec256_t s, t;
	uint8x16_t c7f = vdupq_n_u8(0x7F);
    vset8(s, 0);

    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(t, (const vec256_t *)(row + col));
        vadd8(s, s, t);
        vred8(s, t, c7f);
	 }

    uint32_t k = vhadd8(s);
    return fq_red(k);
} /* end row_acc */

/// accumulates the inverse to a row of length N_pad
/// \param row[in]: pointer to a row
/// \return sum(d[i]**-1) for i in range(N-K)
static inline
FQ_ELEM row_acc_inv(const FQ_ELEM *row) {
    vec256_t a, b, c;
    uint8x16x4_t table[2];
	table[0] = vld1q_u8_x4(fq_inv_table +  0);
	table[1] = vld1q_u8_x4(fq_inv_table + 64);

	const uint8x16_t zero = vdupq_n_u8(0);
	const uint8x16_t mask = vdupq_n_u8(64);

    // NOTE: actually only the last pos need to be 0
    static FQ_ELEM inv_data[N_K_pad] = {0}; 
    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(row + col));

	    b.v[0] = vsubq_u8(a.v[0], mask);
	    b.v[1] = vsubq_u8(a.v[1], mask);
	    const uint8x16_t ret11 = vqtbx4q_u8(zero, table[0], a.v[0]);
	    const uint8x16_t ret12 = vqtbx4q_u8(zero, table[1], b.v[0]);
	    c.v[0] = ret11 ^ ret12;
	    const uint8x16_t ret21 = vqtbx4q_u8(zero, table[0], a.v[1]);
	    const uint8x16_t ret22 = vqtbx4q_u8(zero, table[1], b.v[1]);
	    c.v[1] = ret21 ^ ret22;
        vstore256((vec256_t *)(inv_data + col), c);
	}

    return row_acc(inv_data);
} /* end row_acc_inv */

/// scalar multiplication of a row
/// \param row[in/out] *= s for _ in range(N-K), pointer to a row of length ROUND_UP(N-K, 32)
/// \param s[in]: scalar value
static inline
void row_mul(FQ_ELEM *row,
             const FQ_ELEM s) {
    vec256_t a;
    uint8x16x4_t table[2];
    gf127v_scalar_compute_table(table, s);

    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(row + col));
        a = gf127v_scalar_table_full(a, table);
        vstore256((vec256_t *)(row + col), a);
    }
} /* end row_mul */

/// scalar multiplication of a row
/// \param out[out]: = s*in[i] for i in range(N-K)
/// \param in[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \param s[in]: scalar value
static inline
void row_mul2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in, const FQ_ELEM s) {
    vec256_t a;
    uint8x16x4_t table[2];
    gf127v_scalar_compute_table(table, s);

    for (uint32_t col = 0; (col+32) <= N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in + col));
        a = gf127v_scalar_table_full(a, table);
        vstore256((vec256_t *)(out + col), a);
    }
} /* end row_mul2 */

/// full element wise multiplication of two rows
/// \param out[out]: = in1[i]*in2[i] for i in range(N-K)
/// \param in1[in]: pointer to a row of length ROUND_UP(N-K, 32)
/// \param in2[in]: pointer to a row of length ROUND_UP(N-K, 32)
static inline
void row_mul3(FQ_ELEM *__restrict__ out,
			  const FQ_ELEM *__restrict__ in1,
			  const FQ_ELEM *__restrict__ in2) {
    for (uint32_t col = 0; (col+16) <= N_K_pad; col+=16) {
        const uint8x16_t a = vld1q_u8((uint8_t *)(in1 + col));
        const uint8x16_t b = vld1q_u8((uint8_t *)(in2 + col));
        const uint8x16_t c = gf127v_mul_u128(a, b);
        vst1q_u8((uint8_t *)(out + col), c);
    }
} /* end row_mul3 */

/// invert a row
/// \param out[out]: in[i]**-1 for i in range(N-K)
/// \param in [in]: pointer to a row of length N-K
static inline
void row_inv2(FQ_ELEM *__restrict__ out,
              const FQ_ELEM *__restrict__ in) {
    vec256_t a, b, c;
    uint8x16x4_t table[2];
	table[0] = vld1q_u8_x4(fq_inv_table +  0);
	table[1] = vld1q_u8_x4(fq_inv_table + 64);

	const uint8x16_t zero = vdupq_n_u8(0);
	const uint8x16_t mask = vdupq_n_u8(64);

    for (uint32_t col = 0; col < N_K_pad; col+=32) {
        vload256(a, (vec256_t *)(in + col));

	    b.v[0] = vsubq_u8(a.v[0], mask);
	    b.v[1] = vsubq_u8(a.v[1], mask);
	    const uint8x16_t ret11 = vqtbx4q_u8(zero, table[0], a.v[0]);
	    const uint8x16_t ret12 = vqtbx4q_u8(zero, table[1], b.v[0]);
	    c.v[0] = ret11 ^ ret12;
	    const uint8x16_t ret21 = vqtbx4q_u8(zero, table[0], a.v[1]);
	    const uint8x16_t ret22 = vqtbx4q_u8(zero, table[1], b.v[1]);
	    c.v[1] = ret21 ^ ret22;
        vstore256((vec256_t *)(out + col), c);
	}
} /* end row_inv2 */

/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t row_all_same(const FQ_ELEM *in) {
    vec256_t t1, t2, acc;
    vset8(acc, -1);
    vset8(t2, in[0]);

    uint32_t col = 0;
    for (; col < N_K_pad-32; col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, t2);
        vand(acc, acc, t1);
    }

    const uint32_t t3 = (uint32_t)vmovemask8(acc);
    for (;col < N-K; col++) {
        if (in[col-1] != in[col]) {
            return 0;
        }
    }

    return t3 == -1u;
} /* end row_all_same */

/// NOTE: not CT, but close.
/// \param in[in]: row
/// \return 1 if a zero was found
///         0 else
static inline
uint32_t row_contains_zero(const FQ_ELEM *in) {
    vec256_t t1, t2, acc;
    vset8(t2, 0);
    vset8(acc, 0);
    uint32_t col = 0;
    for (; col < (N_K_pad-32); col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, t2);
        vor(acc, acc, t1);
    }
    
    const uint32_t t3 = (uint32_t)vmovemask8(acc);
    if (t3 != 0ul) { return 1; }

    for (;col < N-K; col++) {
        if (in[col] == 0) {
            return 1;
        }
    }
    return 0;
} /* end row_contains_zero */

/// NOTE: even though the input vectors are of length ROUND_UP(N-K, 32), only
/// N-K bytes are read. As the remaining alignment bytes are zero, which would  
/// alter the result.
/// \param in[in]: vector of length N-K
/// \return: the number of zeros in `in`
static inline
uint32_t row_count_zero(const FQ_ELEM *in) {
    vec256_t t1, zero, acc, mask;
    vset8(zero, 0);
    vset8(acc, 0);
    vset8(mask, 1);
    uint32_t col = 0;
    for (; col < (N_K_pad-32); col += 32) {
        vload256(t1, (vec256_t *)(in + col));
        vcmp8(t1, t1, zero);
        vand(t1, t1, mask);
        vadd8(acc, acc, t1);
    }
  
    uint32_t a = vhadd8(acc);
    for (;col < N-K; col++) {
        a += (in[col] == 0);
    }
    return a;
} /* end row_count_zero */
