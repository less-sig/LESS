/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.0 (February 2022)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "timing_and_stat.h"
#include "LESS.h"
#include "rng.h"


#if defined(CATEGORY_5)
#define NUM_TESTS 6
#elif defined(CATEGORY_3)
#define NUM_TESTS 12
#else
#define NUM_TESTS 24
#endif

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif

welford_t timer;
uint64_t cycles;
pubkey_t pk;
prikey_t sk;
sign_t signature;
char message[8] = "Signme!";

void LESS_sign_mem(void){
    welford_init(&timer);
    for(int i = 0; i <NUM_TESTS; i++) {
        cycles = x86_64_rtdsc();
        LESS_sign(&sk,message,8,&signature);
        welford_update(&timer,(x86_64_rtdsc()-cycles)/1000.0);
    }
    welford_print(timer);
    printf("\n");
}

int main(void){
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"0123456789012345",16);

    LESS_keygen(&sk,&pk);
    LESS_sign_mem();
    return 0;
}
