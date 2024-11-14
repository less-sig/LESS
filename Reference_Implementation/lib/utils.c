/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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

#include "utils.h"
#include <stdlib.h>

/// swaps a and b if
void cswap(uintptr_t *a,
           uintptr_t *b,
           const uintptr_t mask) {
    *a ^= (mask & *b);
    *b ^= (mask & *a);
    *a ^= (mask & *b);
}

/// swaps a and b if f==1, if f==0, nothing will happen
void cswap_bit(uintptr_t *a,
               uintptr_t *b,
               const uintptr_t f) {
	const uint64_t mask = -f;
	cswap(a, b, mask);
}

/// swaps a and b
void cswap_array(uintptr_t *a,
                 uintptr_t *b,
                 const uintptr_t mask,
                 const uint32_t n) {
    for (uint32_t i = 0; i < n; i++) {
        MASKED_SWAP(a[i], b[i], mask);
    }
}

/// Description: Compare two arrays for equality in constant time.
///
/// Arguments:   const uint8_t *a: pointer to first byte array
///              const uint8_t *b: pointer to second byte array
///              size_t len:       length of the byte arrays
///
/// Returns 0 if the byte arrays are equal, 1 otherwise
int verify(const uint8_t *a,
           const uint8_t *b,
           const size_t len) {
    uint8_t r = 0;

    for(size_t i=0;i<len;i++) {
        r |= a[i] ^ b[i];
    }

    return (-(uint64_t)r) >> 63;
}

/// Description: Copy len bytes from x to r if b is 1;
///              don't modify x if b is 0. Requires b to be in {0,1};
///              assumes two's complement representation of negative integers.
///              Runs in constant time.
///
/// Arguments:   uint8_t *r:       pointer to output byte array
///              const uint8_t *x: pointer to input byte array
///              size_t len:       Amount of bytes to be copied
///              uint8_t b:        Condition bit; has to be in {0,1}
void cmov(uint8_t *r, const uint8_t *x, size_t len, uint8_t b) {
#if defined(__GNUC__) || defined(__clang__)
    // Prevent the compiler from
    //    1) inferring that b is 0/1-valued, and
    //    2) handling the two cases with a branch.
    __asm__("" : "+r"(b) : /* no inputs */);
#endif

    b = -b;
    for(size_t i=0;i<len;i++) {
        r[i] ^= b & (r[i] ^ x[i]);
    }
}

// TODO move to optimized
#ifdef USE_AVX
#include <immintrin.h>

void cswap_array(uintptr_t *a,
                 uintptr_t *b,
                 const uintptr_t mask,
                 const uint32_t n) {
    __m256i m = __m256_set1_epi64x(mask);
    __m256i *aa = (__m256i *)a;
    __m256i *bb = (__m256i *)b;
    for (uint32_t i = 0; i < (n/4); i++) {
        MASKED_SWAP(ab[i], bb[i], m);
    }
}

/// taken from kyber
/// Description: Compare two arrays for equality in constant time.
/// Arguments:   const uint8_t *a: pointer to first byte array
///              const uint8_t *b: pointer to second byte array
///              size_t len: length of the byte arrays
///
/// Returns 0 if the byte arrays are equal, 1 otherwise
int verify(const uint8_t *a, const uint8_t *b, size_t len) {
    size_t i;
    uint64_t r;
    __m256i f, g, h;

    h = _mm256_setzero_si256();
    for(i=0;i<len/32;i++) {
        f = _mm256_loadu_si256((__m256i *)&a[32*i]);
        g = _mm256_loadu_si256((__m256i *)&b[32*i]);
        f = _mm256_xor_si256(f,g);
        h = _mm256_or_si256(h,f);
    }
    r = 1 - _mm256_testz_si256(h,h);

    a += 32*i;
    b += 32*i;
    len -= 32*i;

    for(i=0;i<len;i++) {
        r |= a[i] ^ b[i];
    }

    r = (-r) >> 63;
    return r;
}


/// Description: Copy len bytes from x to r if b is 1;
///              don't modify x if b is 0. Requires b to be in {0,1};
///              assumes two's complement representation of negative integers.
///              Runs in constant time.
///
/// Arguments:   uint8_t *r: pointer to output byte array
///              const uint8_t *x: pointer to input byte array
///              size_t len: Amount of bytes to be copied
///              uint8_t b: Condition bit; has to be in {0,1}
void cmov(uint8_t * restrict r, const uint8_t *x, size_t len, uint8_t b){
      size_t i;
      __m256i xvec, rvec, bvec;

#if defined(__GNUC__) || defined(__clang__)
    // Prevent the compiler from
    //    1) inferring that b is 0/1-valued, and
    //    2) handling the two cases with a branch.
    // This is not necessary when verify.c and kem.c are separate translation
    // units, but we expect that downstream consumers will copy this code and/or
    // change how it is built.
    __asm__("" : "+r"(b) : /* no inputs */);
#endif

    bvec = _mm256_set1_epi64x(-(uint64_t)b);
    for(i=0;i<len/32;i++) {
        rvec = _mm256_loadu_si256((__m256i *)&r[32*i]);
        xvec = _mm256_loadu_si256((__m256i *)&x[32*i]);
        rvec = _mm256_blendv_epi8(rvec,xvec,bvec);
        _mm256_storeu_si256((__m256i *)&r[32*i],rvec);
    }

    r += 32*i;
    x += 32*i;
    len -= 32*i;
    for(i=0;i<len;i++) {
        r[i] ^= -b & (x[i] ^ r[i]);
    }
}
#endif

#define MAX_KEYPAIR_INDEX (NUM_KEYPAIRS-1)
#define KEYPAIR_INDEX_MASK ( ((uint16_t)1 << BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX)) -1 )
/* bitmask for rejection sampling of the position */
#define  POSITION_MASK (( (uint16_t)1 << BITS_TO_REPRESENT(T-1))-1)

/* Expands a digest expanding it into a fixed weight string with elements in
 * Z_{NUM_KEYPAIRS}. */
void expand_digest_to_fixed_weight( uint8_t fixed_weight_string[T],
                                    const uint8_t digest[HASH_DIGEST_LENGTH]){
    SHAKE_STATE_STRUCT shake_state;
    initialize_csprng(&shake_state,
                      (const unsigned char *) digest,
                      HASH_DIGEST_LENGTH);

    uint16_t rnd_buf;
    int placed_elements = 0;
    while (placed_elements < W) {
        uint8_t value;
        POSITION_T pos;
        do {
            csprng_randombytes((unsigned char *) &rnd_buf,
                               sizeof(uint16_t),
                               &shake_state);

            value = rnd_buf & (KEYPAIR_INDEX_MASK);
            pos   = rnd_buf >> BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX) ;
            pos   = pos & POSITION_MASK;
        } while ( (value >= NUM_KEYPAIRS) || /* for non-power-of-two keypair numbers */
                  (  pos >= T) ||             /* rejection sampling */
                  (fixed_weight_string[pos] != 0) ); /* skip elements already placed */
        fixed_weight_string[pos] = value;
        placed_elements += (value != 0);
    }
} /* end parse_digest */
