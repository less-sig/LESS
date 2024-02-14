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

void microbench(){
    welford_t timer;
    welford_init(&timer);

    generator_mat_t G;
    generator_rnd(&G);
    uint8_t is_pivot_column[NN];

    uint64_t cycles;
    for(int i = 0; i <NUM_TESTS; i++) {
        cycles = x86_64_rtdsc();
        generator_RREF(&G,is_pivot_column);
        welford_update(&timer,(x86_64_rtdsc()-cycles)/1000.0);
    }
    fprintf(stdout,"Gaussian elimination kCycles (avg,stddev):");
    welford_print(timer);
    printf("\n");

}

void info(){
    fprintf(stdout,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stdout,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stdout,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stdout,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stdout,"Public key %luB\n", sizeof(pubkey_t));
    fprintf(stdout,"Signature: %luB, %f\n", sizeof(sign_t), ((float) sizeof(sign_t))/1024);
}

void LESS_sign_verify_speed(){
    fprintf(stdout,"Computing number of clock cycles as the average of %d runs\n", NUM_TESTS);
    welford_t timer;
    uint64_t cycles;
    pubkey_t pk;
    prikey_t sk;
    sign_t signature;
    char message[8] = "Signme!";
    info();

    printf("Timings (kcycles):\n");
    welford_init(&timer);
    for(int i = 0; i <NUM_TESTS; i++) {
        cycles = x86_64_rtdsc();
        LESS_keygen(&sk,&pk);
        welford_update(&timer,(x86_64_rtdsc()-cycles)/1000.0);
    }
    printf("Key generation kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n");

    welford_init(&timer);
    for(int i = 0; i <NUM_TESTS; i++) {
        cycles = x86_64_rtdsc();
        LESS_sign(&sk,message,8,&signature);
        welford_update(&timer,(x86_64_rtdsc()-cycles)/1000.0);
    }
    printf("Signature kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n");

    int is_signature_ok;
    welford_init(&timer);
    for(int i = 0; i <NUM_TESTS; i++) {
        cycles = x86_64_rtdsc();
        is_signature_ok = LESS_verify(&pk,message,8,&signature);
        welford_update(&timer,(x86_64_rtdsc()-cycles)/1000.0);
    }
    printf("Verification kCycles (avg,stddev):");
    welford_print(timer);
    printf("\n");
    fprintf(stdout,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
}

int iteration = 0;

int main(int argc, char* argv[]){

#ifdef __APPLE__ 
#ifdef __aarch64__
	__m1_setup_rdtsc();
#endif 
#endif
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"0123456789012345",16);
    fprintf(stdout,"LESS reference implementation benchmarking tool\n");
    microbench();
    LESS_sign_verify_speed();
    return 0;
}
