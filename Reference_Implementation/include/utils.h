/**
 *
 * Reference ISO-C11 Implementation of LESS-CF.
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
#include <stddef.h>

#include "parameters.h"

/// NOTE: this implementation does not use additional space
/// \param a[in/out]: will be set to b
/// \param b[in/out]: will be set to a
/// \return nothing
#define SWAP(a, b) { (a)^=(b); (b)^=(a); (a)^=(b); }

/// \param i[in] first index
/// \param j[in] second index
/// \return -1ull if i == j
///         0 else
#define COMPUTE_CT_MASK(i, j) (-(i == j))

#ifdef GCC_SUPPORT_UNROLL_PRAGMA
#  define LOOP_UNROLL_2 _Pragma("GCC unroll 2")
#  define LOOP_UNROLL_4 _Pragma("GCC unroll 4")
#  define LOOP_UNROLL_8 _Pragma("GCC unroll 8")
#elif defined(__clang__)
#  define LOOP_UNROLL_2 _Pragma("clang loop unroll_count(2)")
#  define LOOP_UNROLL_4 _Pragma("clang loop unroll_count(4)")
#  define LOOP_UNROLL_8 _Pragma("clang loop unroll_count(8)")
#else
#  define LOOP_UNROLL_2
#  define LOOP_UNROLL_4
#  define LOOP_UNROLL_8
#endif


/// Compare two arrays for equality in constant time.
/// \param a[in]: pointer to the first byte array
/// \param b[in]: pointer to the second byte array
/// \param len[in]: length of the byte array
/// \returns 0 if the byte arrays are equal, 1 otherwise
int verify(const uint8_t *a,
           const uint8_t *b,
           size_t len);

/// \param fixed_weight_string[out]: array of length T, with hamming weight W,
///     where the values > 0 are between [1, MAX_KEYPAIR_INDEX)
/// \param digest[in]: commitment hash
void SampleChallenge(uint8_t fixed_weight_string[T],
                     const uint8_t digest[HASH_DIGEST_LENGTH]);

