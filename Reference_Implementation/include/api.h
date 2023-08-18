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
#include "LESS.h"

#define CRYPTO_ALGNAME "LESS"

/*  no. of bytes of the secret key */
#define CRYPTO_SECRETKEYBYTES (sizeof(prikey_t))

/*  no. of bytes of the public key */
#define CRYPTO_PUBLICKEYBYTES (sizeof(pubkey_t))

/* no. of bytes of overhead in a signed message */
#define CRYPTO_BYTES (sizeof(sig_t))

/* required bytes of input randomness */
#define  CRYPTO_RANDOMBYTES (SEED_LENGTH_BYTES) //LESS library -  parameters.h


int crypto_sign_keypair(unsigned char *pk,
                        unsigned char *sk
                       );

int crypto_sign(unsigned char *sm,
                unsigned long long *smlen,
                const unsigned char *m,
                unsigned long long mlen,
                const unsigned char *sk
               );

int crypto_sign_open(unsigned char *m,
                     unsigned long long *mlen,
                     const unsigned char *sm,
                     unsigned long long smlen,
                     const unsigned char *pk
                    );
