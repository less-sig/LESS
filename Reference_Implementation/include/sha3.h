/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Floyd Zweydinger <floyd.zweydinger+github@rub.de>
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
#include <stdint.h>

#if CATEGORY == 252
#define SHAKE_STATE_STRUCT shake128incctx
#else
#define SHAKE_STATE_STRUCT shake256incctx
#endif

/// initialized the shake structure
/// \param state[out]: pointer to a uninitialized state 
static inline
void xof_shake_init(SHAKE_STATE_STRUCT *state) {
#if CATEGORY == 252
   shake128_inc_init(state);
#else
   shake256_inc_init(state);
#endif
}

/// absorbs some bytes into the shake struct
/// \param state[in/out]: pointer to an initialized state 
/// \param input[in]: data to fed into shake 
/// \param inputByteLen[in]: length of the input in bytes
static inline
void xof_shake_update(SHAKE_STATE_STRUCT *state,
                      const uint8_t *input,
                      const size_t inputByteLen) {
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

/// finalizes the absorb phase. After this function call, one can extract 
/// random bytes from the shake state.
/// \param state[in]: pointer to an initialized state 
static inline
void xof_shake_final(SHAKE_STATE_STRUCT *state) {
#if CATEGORY == 252
   shake128_inc_finalize(state);
#else
   shake256_inc_finalize(state);
#endif
}

/// \param state[in/out]: pointer to an initialized state 
/// \param output[out]: pointer to the output byte buffer 
/// \param outputByteLen[in]: number of bytes to extract from the state
static inline
void xof_shake_extract(SHAKE_STATE_STRUCT *state,
                       uint8_t *output,
                       const size_t outputByteLen) {
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

