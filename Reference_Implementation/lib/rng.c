/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Floyd Zweydinger <zweydfg+github@rub.de>
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


#include "rng.h"

/* Initializes a CSPRNG from either an input seed or the output of
 * clock_gettime. Input seed assumed to be a C convention string */
SHAKE_STATE_STRUCT platform_csprng_state;

/// \param shake_state[in/out]
/// \param seed[in]
/// \param seed_len_bytes[in]
void initialize_csprng(SHAKE_STATE_STRUCT *shake_state, const unsigned char *seed, const uint32_t seed_len_bytes) {
    // the second parameter is the security level of the SHAKE instance
    xof_shake_init(shake_state, SEED_LENGTH_BYTES * 8);
    xof_shake_update(shake_state, seed, seed_len_bytes);
    xof_shake_final(shake_state);
} /* end initialize_csprng */

void initialize_csprng_ds(SHAKE_STATE_STRUCT *shake_state,
                       const unsigned char *seed,
                       const uint32_t seed_len_bytes,
                       const uint16_t domain_sep_constant ){
    // the second parameter is the security level of the SHAKE instance
    xof_shake_init(shake_state, SEED_LENGTH_BYTES * 8);
    xof_shake_update(shake_state, seed, seed_len_bytes);
#if (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    unsigned char domain_sep[] = {domain_sep_constant >> 8, domain_sep_constant & 0xff};
#else
    unsigned char* domain_sep = (unsigned char *) &domain_sep_constant;
#endif
    xof_shake_update(shake_state, domain_sep, sizeof(uint16_t));
    xof_shake_final(shake_state);
} /* end initialize_csprng_ds */

