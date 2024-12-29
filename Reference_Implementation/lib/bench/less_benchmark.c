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
#include <math.h>

#include "LESS.h"
#include "monomial_mat.h"
#include "rng.h"
#include "cycles.h"


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
double welch_t_statistic(const welford_t state1,
                         const welford_t state2) {
    long double num, den, var1, var2;
    var1 = state1.M2/(long double)(state1.count-1);
    var2 = state2.M2/(long double)(state2.count-1);

    num = state1.mean - state2.mean;
    den = sqrtl(var1/(long double) state1.count + var2/(long double) state2.count );

    return num/den;
}

static inline
void welford_print(const welford_t state) {
    printf("%.2Lf,%.2Lf",
           state.mean,
           sqrtl(state.M2/(long double)(state.count-1)));
}

static inline
long double welford_stddev(const welford_t state) {
    return sqrtl(state.M2/(long double)(state.count-1));
}

static inline
long double welford_mean(const welford_t state) {
    return state.mean;
}

#if defined(CATEGORY_5)
#define NUM_RUNS 6
#elif defined(CATEGORY_3)
#define NUM_RUNS 12
#else
#define NUM_RUNS 32
// #define NUM_RUNS 256
#endif

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif

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
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        LESS_keygen(&sk,&pk);
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Key generation kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n\n");

    
    welford_init(&timer);
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        LESS_sign(&sk,message,8,&signature);
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Signature kCycles (avg,stddev): ");
    welford_print(timer);
    printf("\n\n");
    
    int is_signature_ok;
    welford_init(&timer);
    for(int i = 0; i <NUM_RUNS; i++) {
        cycles = read_cycle_counter();
        is_signature_ok = LESS_verify(&pk,message,8,&signature); // Message never changes
        welford_update(&timer,(read_cycle_counter()-cycles)/1000.0);
    }
    printf("Verification kCycles (avg,stddev):");
    welford_print(timer);
    printf("\n\n");
    fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
}

int iteration = 0;

int main(int argc, char* argv[]){
    (void)argc;
    (void)argv;
    setup_cycle_counter();
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"0123456789012345",16);
    fprintf(stderr,"LESS reference implementation benchmarking tool\n");
    // microbench();
    // monomial_distribution();
    LESS_sign_verify_speed();
    return 0;
}
