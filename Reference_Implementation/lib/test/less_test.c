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
 * TODO explain
 */
void rref_gen_compress_tester(void){
     generator_mat_t G = {0}, Gcheck;
     rref_generator_mat_t SF_G;
     uint8_t is_pivot_column[NN];

     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     memcpy(&Gcheck,&G, sizeof(G));
     generator_rref_compact(&SF_G,&G,is_pivot_column);
     generator_rnd(&G); /* fill with garbage to elicit faults */
     generator_rref_expand(&G,&SF_G);

    if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        printf("Generator SF compression: ko\n");
       fprintf(stderr," Comp-decomp\n");
       generator_pretty_print_name("G",&G);
       fprintf(stderr,"is_pivot = \n [ ");
       for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
       fprintf(stderr,"]\n");

       generator_rref_pretty_print_name("RREF-G",&SF_G);

       fprintf(stderr," Reference\n");
       generator_pretty_print_name("Gcheck",&Gcheck);
    } else {
        printf("Generator SF compression: ok\n");
    }
}

/*
 * TODO explain
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
 * TODO explain
 */
void info(void){
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB, %.3f kiB\n", sizeof(pubkey_t), ((float) sizeof(pubkey_t))/1024);
    fprintf(stderr,"Signature: %luB, %.3f kiB\n", sizeof(sign_t), ((float) sizeof(sign_t))/1024);

}

/* returns 1 if the test is successful, 0 otherwise */
int LESS_sign_verify_test(void){
    pubkey_t pk = {0};
    prikey_t sk = {0};
    sign_t signature;
    char message[8] = "Signme!";
    LESS_keygen(&sk,&pk);
    LESS_sign(&sk,message,8,&signature);
    int is_signature_ok = LESS_verify(&pk,message,8,&signature);
    fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
    return is_signature_ok;
}

#define NUMBER_OF_TESTS 1
#define MLEN 160

/* returns 1 if the test is successful, 0 otherwise */
int LESS_sign_verify_test_multiple(void){
    unsigned long long smlen = 0, mlen1;
    unsigned char *m  = (unsigned char *)calloc(MLEN, sizeof(unsigned char));
    unsigned char *sm = (unsigned char *)calloc(MLEN+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES] = {0}, sk[CRYPTO_SECRETKEYBYTES] = {0};
    ASSERT(m); ASSERT(sm);

    int ret = 0;
    for (size_t i = 0; i < NUMBER_OF_TESTS; ++i) {
        const uint32_t msg_len = rand() % MLEN;
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

        ASSERT(mlen1 == msg_len);
        printf("OK: %zu\n", i);

        ret |= ret_val;
    }

    return ret;
}


int LESS_sign_verify_test_KAT(void) {
    //unsigned char m[] = {0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
	//uint8_t seed[] = {0x06,0x15,0x50,0x23,0x4D,0x15,0x8C,0x5E,0xC9,0x55,0x95,0xFE,0x04,0xEF,0x7A,0x25,0x76,0x7F,0x2E,0x24,0xCC,0x2B,0xC4,0x79,0xD0,0x9D,0x86,0xDC,0x9A,0xBC,0xFD,0xE7,0x05,0x6A,0x8C,0x26,0x6F,0x9E,0xF9,0x7E,0xD0,0x85,0x41,0xDB,0xD2,0xE1,0xFF,0xA1};
    uint8_t seed[] = {0xF1,0x90,0x2A,0x78,0x15,0xF3,0x7B,0xC7,0xF5,0x80,0x2D,0x8C,0xBC,0xE5,0xB4,0x8D,0x82,0xEB,0x85,0x69,0x17,0x18,0x06,0x2B,0xFB,0x84,0xD8,0xC0,0x6A,0xA4,0x1D,0x6E,0x90,0x39,0xB0,0xA1,0x07,0x24,0x5D,0xAF,0xA4,0xEC,0x10,0x9A,0x57,0x33,0x29,0x14};
    uint8_t m[] = {0x1C,0xDF,0x0A,0xE1,0x12,0x47,0x80,0xA8,0xFF,0x00,0x31,0x8F,0x77,0x9A,0x3B,0x86,0xB3,0x50,0x4D,0x05,0x9C,0xA7,0xAB,0x3F,0xE4,0xD6,0xEA,0xE9,0xFD,0x46,0x42,0x8D,0x1D,0xAB,0xB7,0x04,0xC0,0x73,0x5A,0x8F,0xE8,0x70,0x8F,0x40,0x97,0x41,0x01,0x7B,0x72,0x3D,0x9A,0x30,0x4E,0x54,0xFD,0xC5,0x78,0x9A,0x7B,0x07,0x48,0xC2,0x46,0x4B,0x73,0x08,0xAC,0x96,0x65,0x11,0x56,0x44,0xC5,0x69,0xAE,0x25,0x3D,0x52,0x05,0x75,0x13,0x42,0x57,0x4C,0x03,0x34,0x6D,0xDD,0xC1,0x95,0x0A,0x62,0x73,0x54,0x66,0x16,0xB9,0x6D,0x0C,0x5E,0xCE,0x0A,0x04,0x4A,0xF0,0xED,0xEF,0xBE,0x44,0x5F,0x9A,0xE3,0x7D,0xA5,0xAF,0xB8,0xD2,0x2A,0x56,0xD9,0xFD,0x18,0x01,0x42,0x5A,0x0A,0x27,0x6F,0x48,0x43,0x1D,0x7A,0xF0,0x39,0x52,0x1E,0x54,0x95,0x51,0x48,0x13,0x91,0xFE,0x5F,0x4E,0xBF,0xB7,0x64,0x4D,0x9F,0x97,0x82,0xD8,0x3A,0x95,0x13,0x7E,0x84,0xEA,0x3A,0xEB,0x3C,0x2F,0x80,0x99};

    for (uint32_t i=0; i<48; i++) {
        seed[i] = i;
	}

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

    printf("all good\n");
    return 0;
}


#define NUM_TEST_ITERATIONS 10
int main(int argc, char* argv[]){
    (void)argc;
    (void)argv;
    // LESS_sign_verify_test_multiple();
    return LESS_sign_verify_test_KAT();

    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"012345678912345",
                      16);
    fprintf(stderr,"LESS reference implementation functional testbench\n");
    info();
    return LESS_sign_verify_test();

    /// TODO reenable tests
    int tests_ok = 0;
    for (int i = 0; i < NUM_TEST_ITERATIONS; i++) {
        fputc('.',stderr);
      // fprintf(stderr,"test %d: ",i);
       // inverse_mod_tester();
       // gen_by_monom_tester();
       // monomial_tester();
       // rref_gen_compress_tester();
       // gausselim_tester();
       // rref_gen_by_monom_tester();
       // rref_gen_byte_compress_tester();
      tests_ok += LESS_sign_verify_test();
    }
    fprintf(stderr,"%d tests functional out of %d\n",tests_ok,NUM_TEST_ITERATIONS);    
    return 0;
}
