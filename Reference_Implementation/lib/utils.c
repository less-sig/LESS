/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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

#include "utils.h"
#include <stdlib.h>

/// swaps a and b if mask == -1ull
void cswap(uintptr_t *a,
           uintptr_t *b,
           const uintptr_t mask) {
    *a ^= (mask & *b);
    *b ^= (mask & *a);
    *a ^= (mask & *b);
}

/// taken from the kyber impl.
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

///
#define MAX_KEYPAIR_INDEX (NUM_KEYPAIRS-1)
///
#define KEYPAIR_INDEX_MASK (((uint16_t)1u << BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX)) - 1u)
/* bitmask for rejection sampling of the position */
#define  POSITION_MASK (( (uint16_t)1 << BITS_TO_REPRESENT(T-1))-1)

/* Expands a digest expanding it into a fixed weight string with elements in
 * Z_{NUM_KEYPAIRS}. */
void SampleChallenge(uint8_t fixed_weight_string[T],
                     const uint8_t digest[HASH_DIGEST_LENGTH]) {
    SHAKE_STATE_STRUCT shake_state;
    initialize_csprng(&shake_state,
                      (const unsigned char *) digest,
                      HASH_DIGEST_LENGTH);

    uint64_t rnd_buf;
    uint32_t c = 0;
    for (uint32_t i = 0; i < T-W; i++) {
        fixed_weight_string[i] = 0;
    }

    if (NUM_KEYPAIRS != 2) {
        for (uint32_t i = T-W; i < T; i++) {
            uint8_t value;
            do {
                if (c == 0) {
                    csprng_randombytes((unsigned char *) &rnd_buf,
                                     sizeof(uint64_t),
                                     &shake_state);
                    c = 64u / BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX);
                }

                value = rnd_buf & (KEYPAIR_INDEX_MASK);
                rnd_buf >>= BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX);
                c -= 1;
          } while (value >= (NUM_KEYPAIRS-1));
          fixed_weight_string[i] = value + 1;
       }
    } else {
        for (uint32_t i = T-W; i < T; i++) {
            fixed_weight_string[i] = 1;
        }
    }

    for (uint32_t p = T - W; p < T; p++) {
        POSITION_T pos;
        do {
            if (c == 0) {
                csprng_randombytes((unsigned char *) &rnd_buf,
                                   sizeof(uint64_t),
                                   &shake_state);
                c = 64u / BITS_TO_REPRESENT(T-1);
            }
            pos = rnd_buf & (POSITION_MASK);
            rnd_buf >>= BITS_TO_REPRESENT(T-1);
            c -= 1;
        } while (pos > p);
        const uint8_t tmp = fixed_weight_string[p];
        fixed_weight_string[p] = fixed_weight_string[pos];
        fixed_weight_string[pos] = tmp;
    }
}
