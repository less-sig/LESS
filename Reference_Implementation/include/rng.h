/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
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

#include "parameters.h"
#include "sha3.h"
#include <stddef.h>

/// initializes a CSPRNG, given the seed and a state pointer
/// \param shake_state[out]: pointer to a uninitialized state
/// \param seed[in]: seed passed to the prng
/// \param seed_len_bytes[in]: number of bytes of the seed
void initialize_csprng(SHAKE_STATE_STRUCT *shake_state,
                       const unsigned char *seed,
                       uint32_t seed_len_bytes);

/// initializes a CSPRNG, given the seed, a state pointer and a domain separation
/// constant
/// extracts xlen bytes from the CSPRNG, given the state
/// \param shake_state[out] uninitialized shake structure
/// \param seed[in]: seed which is feed into the prng
/// \param seed_len_bytes[in]: length of the seed
/// \param domain_sep_constant[in]: domain seperator as defined in the spec,
///         which is fed into the prng after the seed
void initialize_csprng_ds(SHAKE_STATE_STRUCT *shake_state,
                          const unsigned char *seed,
                          uint32_t seed_len_bytes,
                          uint16_t domain_sep_constant );


/// extracts xlen bytes from the CSPRNG, given the state
/// \param x[out]: output buffer of size
/// \param xlen[in]: number of bytes to extract from the PRNG
/// \param shake_state[out] uninitialized shake structure
static inline
void csprng_randombytes(unsigned char *x,
                        unsigned long long xlen,
                        SHAKE_STATE_STRUCT *shake_state) {
    xof_shake_extract(shake_state, x, xlen);
}

/// global csprng state employed to have a deterministic randombytes for testing
extern SHAKE_STATE_STRUCT platform_csprng_state;

/// extracts `xlen` bytes from the global CSPRNG
/// \param x[out]: output buffer 
/// \param xlen[in]: length of the bytes to extract
static inline
void randombytes(uint8_t *x,
                 const size_t xlen) {
    xof_shake_extract(&platform_csprng_state, x, xlen);
}

/// maybe unused
/// \param seed[in]: seed passed to the global prng
/// \param seed_len_bytes[in]: number of bytes of the seed
__attribute__((unused))
static inline
void init_randombytes(const uint8_t *seed,
                      const size_t seed_len_bytes) {
    initialize_csprng(&platform_csprng_state, seed, seed_len_bytes);
}

void merge_exchange(uint16_t *P,
                    const uint32_t n,
                    SHAKE_STATE_STRUCT *state);
