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

#include <math.h>
#include <stdio.h>

#include "LESS.h"
#include "codes.h"
#include "cycles.h"
#include "monomial_mat.h"
#include "rng.h"
#include "test_helpers.h"
#include "api.h"


typedef struct {
    long double mean;
    long double M2;
    long count;
} welford_t;

static inline
void welford_init(welford_t *state) {
    state->mean = 0.0;
    state->M2 = 0.0;
    state->count = 0;
    return;
}

static inline
void welford_update(welford_t *state, long double sample) {
    long double delta, delta2;
    state->count = state->count + 1;
    delta = sample - state->mean;
    state->mean += delta / (long double)(state->count);
    delta2 = sample - state->mean;
    state->M2 += delta * delta2;
}

static inline
void welford_print(const welford_t state) {
    printf("%.2Lf,%.2Lf",
           state.mean,
           sqrtl(state.M2/(long double)(state.count-1)));
}

#if defined(CATEGORY_5)
#define NUM_RUNS 128
#elif defined(CATEGORY_3)
#define NUM_RUNS 128
#else
#define NUM_RUNS 128
#endif

#define NUM_AVG_RUNS (1u << 10u)

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif


/* samples a random generator matrix */
void generator_rnd(generator_mat_t *res) {
   for(uint32_t i = 0; i < K; i++) {
      rand_range_q_elements(res->values[i], N);
   }
} /* end generator_rnd */


void microbench(void){
    welford_t timer;
    welford_init(&timer);

    generator_mat_t G;
    generator_rnd(&G);
    uint8_t is_pivot_column[NN];

    uint64_t cycles;
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        generator_RREF(&G,is_pivot_column);
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    fprintf(stderr,"Gaussian elimination kCycles (avg,stddev):");
    welford_print(timer);
    printf("\n");

}

void info(void){
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB\n", sizeof(pubkey_t));
    fprintf(stderr,"Signature: %luB, %f\n", sizeof(sign_t), ((float) sizeof(sign_t))/1024);
}

int LESS_avg_sign_size(void) {
    uint8_t seed[] ={0x83,0xC6,0x53,0x70,0x8F,0xAF,0x3E,0x5F,0x6F,0xBC,0x9D,0xFB,0xE6,0xFB,0x5E,0x83,0xE5,0x72,0xA7,0x68,0x86,0x45,0xD7,0x5D,0x2C,0x48,0x35,0xB2,0x86,0x95,0xDE,0xA4,0xBD,0x70,0x93,0x74,0x0D,0x0F,0xF4,0x32,0x37,0x35,0x4E,0xAD,0x1C,0x97,0x8B,0xC2};
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)seed,
                      48);

    const uint32_t mlen = 80;
    unsigned long long smlen = 0;
    unsigned char *m = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char *sm = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char pk[CRYPTO_PUBLICKEYBYTES] = {0}, sk[CRYPTO_SECRETKEYBYTES] = {0};

    int ret_val;
    uint64_t size = 0;
    for(size_t i = 0; i < NUM_AVG_RUNS; i++) {
        init_randombytes(m, mlen);
        if ((ret_val = crypto_sign_keypair(pk, sk)) != 0) {
            return -1;
        }
        if ( (ret_val = crypto_sign(sm, &smlen, m, mlen, sk)) != 0) {
            printf("crypto_sign returned <%d>\n", ret_val);
            return -1;
        }

        size += (smlen - mlen);
    }

    double avg = ((double) size)/((double)NUM_AVG_RUNS);
    printf("WORST sig size: %ld\n", CRYPTO_BYTES);
    printf("AVG   sig size: %f\n", avg);

    free(m);
    free(sm);
    return 0;
}

void LESS_sign_verify_speed(void){
    fprintf(stderr,"Computing number of clock cycles as the average of %d runs\n", NUM_RUNS);
    welford_t timer;
    uint64_t cycles;
    pubkey_t pk;
    prikey_t sk;
    sign_t signature;
    char message[8] = "Signme!";
    info();

    printf("Timings (kcycles):\n");
    welford_init(&timer);
    for(size_t i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        LESS_keygen(&sk,&pk);
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Key generation kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n");


    welford_init(&timer);
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        LESS_sign(&sk,message,8,&signature);
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Signature kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n");

    int is_signature_ok = 1;
    welford_init(&timer);
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        is_signature_ok = LESS_verify(&pk,message,8,&signature); // Message never changes
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Verification kCycles (avg,stddev):");
    welford_print(timer);
    printf("\n");
    fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
}

int main(int argc, char* argv[]){
    (void)argc;
    (void)argv;
    setup_cycle_counter();
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"0123456789012345",16);
    fprintf(stderr,"LESS implementation benchmarking tool\n");
    // microbench();
    LESS_sign_verify_speed();
    //LESS_avg_sign_size();
    return 0;
}
