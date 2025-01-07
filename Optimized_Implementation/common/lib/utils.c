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
#include <string.h>
#include <stdlib.h>


/// swaps a and b if
void cswap(uintptr_t *a,
           uintptr_t *b,
           const uintptr_t mask) {
    *a ^= (mask & *b);
    *b ^= (mask & *a);
    *a ^= (mask & *b);
}

#include <immintrin.h>

/// taken from kyber
/// Description: Compare two arrays for equality in constant time.
/// Arguments:   const uint8_t *a: pointer to first byte array
///              const uint8_t *b: pointer to second byte array
///              size_t len: length of the byte arrays
///
/// Returns 0 if the byte arrays are equal, 1 otherwise
int verify(const uint8_t *a,
           const uint8_t *b,
           size_t len) {
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
    r = 1u - (uint32_t)_mm256_testz_si256(h,h);

    a += 32*i;
    b += 32*i;
    len -= 32*i;

    for(i=0;i<len;i++) {
        r |= a[i] ^ b[i];
    }

    r = (-r) >> 63;
    return r;
}

#define MAX_KEYPAIR_INDEX (NUM_KEYPAIRS-1)
#define KEYPAIR_INDEX_MASK ( ((uint16_t)1 << BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX)) -1 )
/* bitmask for rejection sampling of the position */
#define  POSITION_MASK (( (uint16_t)1 << BITS_TO_REPRESENT(T-1))-1)

/* Expands a digest expanding it into a fixed weight string with elements in
 * Z_{NUM_KEYPAIRS}. */
void DigestToFixedWeight(uint8_t fixed_weight_string[T],
                         const uint8_t digest[HASH_DIGEST_LENGTH]){
   SHAKE_STATE_STRUCT shake_state;
   initialize_csprng(&shake_state,
                     (const unsigned char *) digest,
                     HASH_DIGEST_LENGTH);

   uint16_t rnd_buf;
   for (int i = 0; i < T-W; i++) 
      fixed_weight_string[i] = 0;

    if (NUM_KEYPAIRS != 2) {
        for (int i = T-W; i < T; i++) {
            uint8_t value;
            do {
                csprng_randombytes((unsigned char *) &rnd_buf,
                                 sizeof(uint8_t),
                                 &shake_state);

                value = rnd_buf & (KEYPAIR_INDEX_MASK);
          } while (value >= NUM_KEYPAIRS-1);
          fixed_weight_string[i] = value + 1;
       }
    } else {
        for (int i = T-W; i < T; i++) {
            fixed_weight_string[i] = 1;
        }
   }

   for (int p = T-W; p < T; p++) {
      POSITION_T pos;
      uint8_t tmp;
      do {
         csprng_randombytes((unsigned char *) &rnd_buf,
                             sizeof(POSITION_T),
                             &shake_state);

         pos = rnd_buf & (POSITION_MASK);
      } while (pos > p);
      tmp = fixed_weight_string[p];
      fixed_weight_string[p] = fixed_weight_string[pos];
      fixed_weight_string[pos] = tmp;
   }      
} /* end parse_digest */
