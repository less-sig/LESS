/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
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



#include "api.h"
#include <stddef.h>  //  size_t
#include <string.h>  // memcpy
/*----------------------------------------------------------------------------*/

int crypto_sign_keypair(unsigned char *pk,
                        unsigned char *sk) {
    /* keygen cannot fail */
    LESS_keygen((prikey_t *) sk,
                (pubkey_t *) pk);

    return 0; // NIST convention: 0 == zero errors
} // end crypto_sign_keypair

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*... generating a signed message sm[0],sm[1],...,sm[*smlen-1]                */
/*... from original message m[0],m[1],...,m[mlen-1]                           */
/*... under secret key sk[0],sk[1],...                                        */
int crypto_sign(unsigned char *sm,
                unsigned long long *smlen,                        // out parameter
                const unsigned char *m, unsigned long long mlen,  // in parameter
                const unsigned char *sk)                          // in parameter
{
    /* sign cannot fail */
    memcpy((unsigned char *) sm, (const unsigned char *) m, (size_t) mlen);

    /// TODO: we need to fix this somehow: But how do we do this?
    /// The problem, iirc, is that always a different amount of seeds are opened.
    /// So we need to keep track of them
    /// TODO second problem: the `without tree` version does not return seeds
    const uint32_t num_seeds_published = LESS_sign(
            (const prikey_t *)
                    sk,                                             // in parameter
            (const char *const) m, (const uint64_t) mlen,           // in parameter
            (sign_t *) (sm + mlen));                                // out parameter
#if defined(SEED_TREE)
    const uint32_t sig_len = LESS_CRYPTO_BYTES(num_seeds_published);
#else
    (void)num_seeds_published;
    const uint32_t sig_len = sizeof(sign_t);//LESS_CRYPTO_BYTES;
#endif

    *smlen = mlen + sig_len;
    return 0;  // NIST convention: 0 == zero errors
} // end crypto_sign

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*.  ... verifying a signed message sm[0],sm[1],...,sm[smlen-1]               */
/*.  ... under public key pk[0],pk[1],...                                     */
/*.  ... and producing original message m[0],m[1],...,m[*mlen-1]              */
int crypto_sign_open(unsigned char *m,
                     unsigned long long *mlen,        // out parameter
                     const unsigned char *sm, unsigned long long smlen, // in parameter
                     const unsigned char *pk)                           // in parameter
{
#if defined(SEED_TREE)
    /// TODO, see sign for explanation
    const uint8_t num_seeds_published = sm[smlen - 1u];
    if (num_seeds_published >= MAX_PUBLISHED_SEEDS) {
        return -1;
    }

    // the size of the signature in bytes
    const uint32_t sig_len = LESS_CRYPTO_BYTES((uint32_t)num_seeds_published);
#else
    const uint32_t sig_len = sizeof(sign_t);//LESS_CRYPTO_BYTES;
#endif

    *mlen = smlen - (unsigned long long)sig_len;
    memcpy((unsigned char *) m, (const unsigned char *) sm, (size_t) *mlen);
    /* verify returns 1 if signature is ok, 0 otherwise */
    int ok = LESS_verify((const pubkey_t *const)
                                 pk,                                    // in parameter
                         (const char *const) m, (const uint64_t) *mlen, // in parameter
                         (const sign_t *const) (sm + *mlen));           // in parameter

    // NIST convention: 0 == zero errors, -1 == error condition
    return ok - 1;
} // end crypto_sign_open

/*----------------------------------------------------------------------------*/
