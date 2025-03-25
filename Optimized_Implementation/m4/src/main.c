#ifndef UNIT_TEST

#include <inttypes.h>
#include <stdio.h>

#if !defined(__APPLE__) && !defined (__x86_64__)
#include "stm32f4xx_hal.h"
#include "m4_utils.h"
#else
static long long get_cycles(void) {
	unsigned long long result;
	asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
	             : "=a"(result)::"%rdx");
	return result;
}
#endif


// LESS includes
#include "codes.h"
#include "LESS.h"
#include "rng.h"
#include "api.h"



#define MAX_STACK_SIZE 132768
#define MSG_LEN 80

const int N_BENCH = 20;     /* Number of tests. */
//const int MSG_LEN = 80;       /* Message length. */

#ifdef USE_M4
/// \brief 
/// \param p 
/// \param n 
void randombytes(unsigned char *p, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        p[i] = i;
    }
}
#endif

/* This goes outside of 'main' to avoid stack overflows. */
int bench_less() {
    int i;
    uint32_t cc_keyg[N_BENCH];
    uint32_t cc_sign[N_BENCH];
    uint32_t cc_verf[N_BENCH];
    uint64_t begin, end;
    printf("\n\nBenchmarks...\n\n");
    
    for (i = 0; i < N_BENCH; i++) {
        pubkey_t pk;
        prikey_t sk;
        sign_t signature;
        char msg[MSG_LEN] = {0};

        /* Generate a random message. */
        randombytes((unsigned char *)msg, MSG_LEN);

        begin = get_cycles();
        {
            LESS_keygen(&sk, &pk);
        }
        end = get_cycles();
        cc_keyg[i] = end - begin;
        

        begin = get_cycles();
        {
           LESS_sign(&sk, msg, MSG_LEN, &signature);
        }
        end = get_cycles();
        cc_sign[i] = end - begin;

        begin = get_cycles();
        {
            /* Verify the message */
            if (LESS_verify(&pk, msg, 8, &signature) != 0) {
                printf("Error: Verification failed!\n");
                return -1;
            }
        }
        end = get_cycles();
        cc_verf[i] = end - begin;
    }

    printf("keygen ");
    for (i = 0; i < N_BENCH; i++) {
        printf("%" PRIu32" ", cc_keyg[i]);
    }

    printf("\n\nsign ");
    for (i = 0; i < N_BENCH; i++) {
        printf("%" PRIu32" ", cc_sign[i]);
    }

    printf("\n\nverf ");
    for (i = 0; i < N_BENCH; i++) {
        printf("%" PRIu32" ", cc_verf[i]);
    }
    printf("\n\nDONE!");
    
    return 0;
}

int test() {
    uint8_t seed[48] = {0};
    uint8_t m[48] = {0};
    // init_randombytes(seed, 48);
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)seed,
                      48);

    const uint32_t mlen = sizeof(m);
    unsigned long long smlen = 0, mlen1;
    unsigned char *m1 = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char *sm = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char pk[CRYPTO_PUBLICKEYBYTES] = {0}, sk[CRYPTO_SECRETKEYBYTES] = {0};

    int ret_val;
    if ((ret_val = crypto_sign_keypair(pk, sk)) != 0) {
        printf("crypto_sign_keypair returned <%d>\n", ret_val);
        return -1;
    }
    if ( (ret_val = crypto_sign(sm, &smlen, m, mlen, sk)) != 0) {
        printf("crypto_sign returned <%d>\n", ret_val);
        return -1;
    }
    if ( (ret_val = crypto_sign_open(m1, &mlen1, sm, smlen, pk)) != 0) {
        printf("crypto_sign_open returned <%d>\n", ret_val);
        return -1;
    }

    free(m1);
    free(sm);
    printf("all good\n");
    return 0;
}

int main(void) {
#if !defined(__APPLE__) && !defined (__x86_64__)
    utils_init();
    bm_decls;

    while (1) {  
        //__disable_irq();
        bm_start();

        //LED_toggle();
        bench_less();

        bm_end();
        //__enable_irq();

        HAL_Delay(1000);
    }
#else
    bench_less();
    //test();
#endif

}

#endif  // UNIT_TEST
