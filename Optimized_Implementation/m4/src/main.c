#ifndef UNIT_TEST

#include <inttypes.h>
#include <stdio.h>

#if !defined(__APPLE__) && !defined (__x86_64__)
#include "stm32f4xx_hal.h"
#else 
static long long get_cycles(void) {
	unsigned long long result;
	asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
	             : "=a"(result)::"%rdx");
	return result;
}
#endif

// non LESS includes
#include "m4_utils.h"

// LESS includes
#include "LESS.h"


#define MAX_STACK_SIZE 132768
const int N_BENCH = 1;     /* Number of tests. */
const int MSG_LEN = 80;       /* Message length. */

/// \brief 
/// \param p 
/// \param n 
void randombytes(unsigned char *p, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        p[i] = i;
    }
}

/* This goes outside of 'main' to avoid stack overflows. */
int bench_less() {
    int i;
    uint32_t cc_keyg[N_BENCH];
    uint32_t cc_sign[N_BENCH];
    uint32_t cc_verf[N_BENCH];
    uint32_t begin, end;
    printf("\n\nBenchmarks...\n\n");
    
    for (i = 0; i < N_BENCH; i++) {
        pubkey_t pk;
        prikey_t sk;
        sign_t signature;
        char msg[MSG_LEN];

        /* Generate a random message. */
        randombytes((unsigned char *)msg, MSG_LEN);

        begin = get_cycles();{
            LESS_keygen(&sk,&pk);
        }
        end = get_cycles();
        printf("kek1\n");
        cc_keyg[i] = end - begin;
        

        begin = get_cycles();
        {
           LESS_sign(&sk, msg, MSG_LEN, &signature);
        }
        end = get_cycles();
        printf("kek2\n");
        cc_sign[i] = end - begin;

        begin = get_cycles();
        {
            /* Verify the message */
            //if (LESS_verify(&pk, msg, 8, &signature) != 0) {
            //    printf("Error: Verification failed!\n");
            //    return -1;
            //}
        }
        end = get_cycles();
        printf("kek3\n");
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

int main(void) {
#if !defined(__APPLE__) && !defined (__x86_64__)
    utils_init();
    bm_decls;

    while (1) {  
        __disable_irq();
        bm_start();

        //LED_toggle();
        bench_less();

        bm_end();
        __enable_irq();

        HAL_Delay(1000);
    }
#else
    bench_less();
#endif

}

#endif  // UNIT_TEST
