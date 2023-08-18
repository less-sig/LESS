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

void compute_digest( uint8_t digest[HASH_DIGEST_LENGTH],
                     const char *const m,
                     const uint64_t mlen,
                     normalized_IS_t V_array[T],
                     const unsigned char tree_salt[HASH_DIGEST_LENGTH]) {
   unsigned char* message_buffer = malloc( T*sizeof(normalized_IS_t)+ mlen + HASH_DIGEST_LENGTH);

   memcpy(message_buffer, V_array, T*sizeof(normalized_IS_t) );
   memcpy(message_buffer + T*sizeof(normalized_IS_t), m, mlen);
   memcpy(message_buffer + (T*sizeof(normalized_IS_t) + mlen),
          tree_salt,
          HASH_DIGEST_LENGTH);

#if (HASH_DIGEST_LENGTH*8 == 256)
   sha3_256(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#elif (HASH_DIGEST_LENGTH*8 == 384)
   sha3_384(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#elif (HASH_DIGEST_LENGTH*8 == 512)
   sha3_512(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#else
#error digest length unsupported by SHA-3
#endif
   free(message_buffer);
} /* end compute_digest */

void compute_digest_array( uint8_t digest[HASH_DIGEST_LENGTH],
                     const char *const m,
                     const uint64_t mlen,
                     uint8_t V_array_compressed[T*RREF_IS_COLUMNS_PACKEDBYTES],
                     const unsigned char tree_salt[HASH_DIGEST_LENGTH]) {
   unsigned char* message_buffer = malloc( T*sizeof(normalized_IS_t)+ mlen + HASH_DIGEST_LENGTH);
  
   memcpy(message_buffer, V_array_compressed, T*RREF_IS_COLUMNS_PACKEDBYTES );
   memcpy(message_buffer + T*sizeof(normalized_IS_t), m, mlen);
   memcpy(message_buffer + (T*sizeof(normalized_IS_t) + mlen), 
          tree_salt, 
          HASH_DIGEST_LENGTH);

#if (HASH_DIGEST_LENGTH*8 == 256)
   sha3_256(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#elif (HASH_DIGEST_LENGTH*8 == 384)
   sha3_384(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#elif (HASH_DIGEST_LENGTH*8 == 512)
   sha3_512(digest,message_buffer, mlen+T*sizeof(normalized_IS_t));
#else
#error digest length unsupported by SHA-3
#endif
   free(message_buffer);
} /* end compute_digest */


void sha3_hash(uint8_t digest[HASH_DIGEST_LENGTH],
               const unsigned char *const m,
               const uint64_t mlen)
{
#if (HASH_DIGEST_LENGTH*8 == 256)
   sha3_256(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 384)
   sha3_384(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 512)
   sha3_512(digest,m, mlen);
#else
#error digest length unsupported by SHA-3
#endif
} /* end compute_digest */


#define MAX_KEYPAIR_INDEX (NUM_KEYPAIRS-1)
#define KEYPAIR_INDEX_MASK ( ((uint16_t)1 << BITS_TO_REPRESENT(MAX_KEYPAIR_INDEX)) -1 )
/* bitmask for rejection sampling of the position */
#define  POSITION_MASK (( (uint16_t)1 << BITS_TO_REPRESENT(T-1))-1)

/* Expands a digest expanding it into a fixed weight string with elements in
 * Z_{NUM_KEYPAIRS}. */
#include <stdio.h>
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
