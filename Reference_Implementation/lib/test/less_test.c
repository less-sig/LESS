#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fq_arith.h"
#include "monomial_mat.h"
#include "codes.h"
#include "LESS.h"
#include "rng.h"
#include "api.h"

#include "test_helpers.c"

#define GRN "\e[0;32m"
#define WHT "\e[0;37m"

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif

/// taken from the nist package.
void fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L) {
    unsigned long long  i;
    fprintf(fp, "%s", S);
    for ( i=0; i<L; i++ )
        fprintf(fp, "%02X", A[i]);
    if ( L == 0 )
        fprintf(fp, "00");
    fprintf(fp, "\n");
}

/* Exhaustive testing of inverses mod Q */
void inverse_mod_tester(void){
    uint32_t value[Q-1];
    uint32_t inverse[Q-1];
    for(uint32_t i=1; i <= Q-1; i++){
        value[i-1] = i;
        inverse[i-1] = fq_inv(i);
    }
    int all_ok = 1;
    for(uint32_t i=1; i <= Q-1; i++){
        if((value[i-1]*inverse[i-1]) % Q !=1){
           printf("%u*%u=%u\n",
                  value[i-1],
                  inverse[i-1],
                  (value[i-1]*inverse[i-1])%Q);
           all_ok = 0;
        }
    }
    if (all_ok){
        puts("All inverses on F_q ok");
    }
}

/*
 *
 */
void rref_gen_byte_compress_tester(void){
     generator_mat_t G = {0}, Gcheck;
     uint8_t G_compressed [RREF_MAT_PACKEDBYTES];
     uint8_t is_pivot_column[NN];

     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     memcpy(&Gcheck,&G, sizeof(G));
     compress_rref(G_compressed,&G,is_pivot_column);
     generator_rnd(&G); /* fill with garbage to elicit faults */
     expand_to_rref(&G,G_compressed,is_pivot_column);

    if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        printf("Generator SF byte compression: ko\n");
       fprintf(stderr," Comp-decomp\n");
       generator_pretty_print_name("G",&G);

       fprintf(stderr,"is_pivot = \n [ ");
       for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
       fprintf(stderr,"]\n");

       fprintf(stderr," \n\n\n\n\n\n\n\n\nReference\n");
       generator_pretty_print_name("Gcheck",&Gcheck);
    } else {
        printf("Generator SF compression: ok\n");
    }
}

/*
 *
 */
void info(void){
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB, %.3f kiB\n", sizeof(pubkey_t), ((float) sizeof(pubkey_t))/1024);
    fprintf(stderr,"Signature: %luB, %.3f kiB\n", sizeof(sign_t), ((float) sizeof(sign_t))/1024);

}


#define NUMBER_OF_TESTS 1
#define MLEN 160

/* returns 1 if the test is successful, 0 otherwise */
int LESS_sign_verify_test_multiple(void){
    unsigned long long smlen = 0, mlen1;
    unsigned char *m  = (unsigned char *)calloc(MLEN, sizeof(unsigned char));
    unsigned char *sm = (unsigned char *)calloc(MLEN+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES] = {0}, sk[CRYPTO_SECRETKEYBYTES] = {0};

    int ret = 0;
    for (size_t i = 0; i < NUMBER_OF_TESTS; ++i) {
        const uint32_t msg_len = (uint32_t)rand() % MLEN;
        randombytes(m, msg_len);

        int ret_val;
        if ((ret_val = crypto_sign_keypair(pk, sk)) != 0) {
            printf("crypto_sign_keypair returned <%d>\n", ret_val);
            return -1;
        }

        fprintBstr(stdout, "pk = ", pk, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(stdout, "sk = ", sk, CRYPTO_SECRETKEYBYTES);

        if ( (ret_val = crypto_sign(sm, &smlen, m, msg_len, sk)) != 0) {
            printf("crypto_sign returned <%d>\n", ret_val);
            return -1;
        }
        fprintBstr(stdout, "sm = ", sm, smlen);
        if ( (ret_val = crypto_sign_open(m, &mlen1, sm, smlen, pk)) != 0) {
            printf("crypto_sign_open returned <%d>\n", ret_val);
            return -1;
        }

        if(mlen1 != msg_len) {
            printf("crypto_sign_open wrong length\n");
            return -1;
        }
        printf("OK: %zu\n", i);

        ret |= ret_val;
    }

    return ret;
}


int LESS_sign_verify_test_KAT(void) {
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
    fprintBstr(stdout, "pk = ", pk, CRYPTO_PUBLICKEYBYTES);
    fprintBstr(stdout, "sk = ", sk, CRYPTO_SECRETKEYBYTES);

    if ( (ret_val = crypto_sign(sm, &smlen, m, mlen, sk)) != 0) {
        printf("crypto_sign returned <%d>\n", ret_val);
        return -1;
    }

    fprintBstr(stdout, "sm = ", sm, smlen);
    if ( (ret_val = crypto_sign_open(m1, &mlen1, sm, smlen, pk)) != 0) {
        printf("crypto_sign_open returned <%d>\n", ret_val);
        return -1;
    }

    free(m1);
    free(sm);
    printf("all good\n");
    return 0;
}


#define NUM_TEST_ITERATIONS 10
int main(int argc, char* argv[]){
    (void)argc;
    (void)argv;
    // LESS_sign_verify_test_multiple();
    return LESS_sign_verify_test_KAT();
}
