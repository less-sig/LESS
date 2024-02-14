/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Floyd Zweydinger
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * 
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
#include "monomial_mat.h"
#include "codes.h"
#include <stdint.h>



/* Public key: the first gen. matrix is shrunk to just a seed, all the
 * others are stored in RREF form  */
// __attribute__((packed))
typedef struct {
   unsigned char G_0_seed[SEED_LENGTH_BYTES];
   uint8_t SF_G [NUM_KEYPAIRS-1][RREF_MAT_PACKEDBYTES];
} pubkey_t;

/* Private key: it contains both a single seed generating all private *
 * (inverse) monomials and the seed to geneate the public code */
// __attribute__((packed))
typedef struct {
   /*the private key is compressible down to a single seed*/
   unsigned char compressed_sk[PRIVATE_KEY_SEED_LENGTH_BYTES];
   unsigned char G_0_seed[SEED_LENGTH_BYTES];
} prikey_t;

// __attribute__((packed))
typedef struct sig_t {
    unsigned char tree_salt[HASH_DIGEST_LENGTH];
    uint8_t monom_actions[W][MONO_ACTION_PACKEDBYTES];
    uint8_t digest[HASH_DIGEST_LENGTH];
    /// we need an additional byte to store the number of published seeds
    unsigned char seed_storage[SEED_TREE_MAX_PUBLISHED_BYTES + 1u];
} sign_t;

/* keygen cannot fail */
void LESS_keygen(prikey_t *SK,
                 pubkey_t *PK);

/* sign cannot fail */
size_t LESS_sign(const prikey_t *SK,
               const char *const m,
               const uint64_t mlen,
               sign_t *sig);

/* verify returns 1 if signature is ok, 0 otherwise */
int LESS_verify(const pubkey_t *const PK,
                const char *const m,
                const uint64_t mlen,
                const sign_t *const sig);
