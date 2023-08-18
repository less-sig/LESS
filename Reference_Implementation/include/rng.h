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

#include "parameters.h"
#include "sha3.h"

/* initializes a CSPRNG, given the seed and a state pointer */
void initialize_csprng(SHAKE_STATE_STRUCT *shake_state,
                       const unsigned char *seed,
                       const uint32_t seed_len_bytes);

/* extracts xlen bytes from the CSPRNG, given the state */
static inline
void csprng_randombytes(unsigned char *x,
                        unsigned long long xlen,
                        SHAKE_STATE_STRUCT *shake_state)
{
   xof_shake_extract(shake_state,x,xlen);
}

/* global csprng state employed to have a deterministic randombytes for testing */
extern SHAKE_STATE_STRUCT platform_csprng_state;
/* extracts xlen bytes from the global CSPRNG */
static inline
void randombytes(unsigned char *x,
                 unsigned long long xlen)
{
   xof_shake_extract(&platform_csprng_state,x,xlen);
}
