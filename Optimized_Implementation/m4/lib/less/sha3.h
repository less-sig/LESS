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

#pragma once

#include "fips202.h"
#include <stddef.h>

/* standalone SHA-3 implementation has no visible state for single-call SHA-3 */
// #define SHA3_STATE_STRUCT shake256ctx
/* and has different states for SHAKE depending on security level*/
#if CATEGORY == 252
#define SHAKE_STATE_STRUCT shake128incctx
#else
#define SHAKE_STATE_STRUCT shake256incctx
#endif
// %%%%%%%%%%%%%%%%%% Self-contained SHAKE Wrappers %%%%%%%%%%%%%%%%%%%%%%%%%%%%

static inline
void xof_shake_init(SHAKE_STATE_STRUCT *state, int val) {
	(void)val;
#if CATEGORY == 252
   shake128_inc_init(state);
#else
   shake256_inc_init(state);
#endif
}

static inline
void xof_shake_update(SHAKE_STATE_STRUCT *state,
                      const unsigned char *input,
                      size_t inputByteLen) {
#if CATEGORY == 252
   shake128_inc_absorb(state,
                       (const uint8_t *)input,
                       inputByteLen);
#else
   shake256_inc_absorb(state,
                       (const uint8_t *)input,
                       inputByteLen);
#endif
}

static inline
void xof_shake_final(SHAKE_STATE_STRUCT *state) {
#if CATEGORY == 252
   shake128_inc_finalize(state);
#else
   shake256_inc_finalize(state);
#endif
}

static inline
void xof_shake_extract(SHAKE_STATE_STRUCT *state,
                       unsigned char *output,
                       unsigned int outputByteLen) {
#if CATEGORY == 252
   shake128_inc_squeeze(output, outputByteLen, state);
#else
   shake256_inc_squeeze(output, outputByteLen, state);
#endif
}


// This abstract away the SHA3 interface.
#if (HASH_DIGEST_LENGTH*8 == 256)
#define LESS_SHA3_INC_CTX                     sha3_256incctx
#define LESS_SHA3_INC_INIT(state)             sha3_256_inc_init(state)
#define LESS_SHA3_INC_ABSORB(state, ptr, len) sha3_256_inc_absorb(state, ptr, len)
#define LESS_SHA3_INC_FINALIZE(output, state) sha3_256_inc_finalize(output, state)
#elif (HASH_DIGEST_LENGTH*8 == 384)
#define LESS_SHA3_INC_CTX                     sha3_384incctx
#define LESS_SHA3_INC_INIT(state)             sha3_384_inc_init(state)
#define LESS_SHA3_INC_ABSORB(state, ptr, len) sha3_384_inc_absorb(state, ptr, len)
#define LESS_SHA3_INC_FINALIZE(output, state) sha3_384_inc_finalize(output, state)
#elif (HASH_DIGEST_LENGTH*8 == 512)
#define LESS_SHA3_INC_CTX                     sha3_512incctx
#define LESS_SHA3_INC_INIT(state)             sha3_512_inc_init(state)
#define LESS_SHA3_INC_ABSORB(state, ptr, len) sha3_512_inc_absorb(state, ptr, len)
#define LESS_SHA3_INC_FINALIZE(output, state) sha3_512_inc_finalize(output, state)
#else
#error digest length unsupported by SHA-3
#endif

