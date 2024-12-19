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

/* Test monomial matrix multiplication and inversion testing if
 * M*M^-1 = I, where I is a random monomial matrix */
void monomial_tester(void){
    monomial_t mat1, mat2, id,idcheck;
    monomial_mat_id(&idcheck);
    monomial_mat_rnd(&mat1);
    monomial_mat_inv(&mat2,&mat1);
    monomial_mat_mul(&id,&mat2,&mat1);
    if( memcmp( &id,&idcheck,sizeof(monomial_t)) !=0 ){
       monomial_mat_pretty_print_name("M1",&mat1);
       monomial_mat_pretty_print_name("M1^-1",&mat2);
       monomial_mat_pretty_print_name("M1^-1 * M1",&id);
    } else {
        printf("Monomial arith test: ok\n");
    }
}

/*Generate a random full-rank G, keep multiplying G by a monomial and RREF'ing */
#define NUMBER_OF_GE_TESTS 10000
void gausselim_tester(void){
    generator_mat_t G,GMul;
    uint8_t is_pivot_column[NN];

    /* randomly generate a non-singular G */
    do {
        generator_rnd(&G);
        memset(is_pivot_column,0,sizeof(is_pivot_column));
    } while ( generator_RREF(&G,is_pivot_column) == 0);

    /* Stress-test GE by repeatedly bringing in SF GQ*/
    monomial_t mat1;
    int full_rank,rref_ok, all_ok = 1;
    for(int i=0; (i<NUMBER_OF_GE_TESTS) && all_ok; i++ ){
        monomial_mat_rnd(&mat1);
        generator_monomial_mul(&GMul,&G,&mat1);
        memcpy(&G,&GMul,sizeof(generator_mat_t));
        memset(is_pivot_column,0,sizeof(is_pivot_column));
        full_rank = generator_RREF(&GMul,is_pivot_column);
        if(!full_rank){
            all_ok = 0;
            fprintf(stderr,"Singular Matrix (iter:%d)\n",i);
            generator_pretty_print_name("Pre-GE",&G);
            generator_pretty_print_name("Post-GE",&GMul);
        }

        /* check if the matrix is in reduced row echelon form i.e,
         *  - there are k pivots
         *  - each pivot appears as the first element of a row
         *  - is_pivot_column indicator array is correct
         */
        rref_ok = 1;
        for(int row_idx = 0; row_idx < K; row_idx++){
            int found_pivot_column = 0;
            while ( (GMul.values[row_idx][found_pivot_column] == 0) &&
                   (found_pivot_column < N) ){
                found_pivot_column++;
            }
            if ( (GMul.values[row_idx][found_pivot_column] != 1) ){
                fprintf(stderr,"row %d Pivot actually equal to %d\n",row_idx, GMul.values[row_idx][found_pivot_column]);
                     rref_ok = 0;
                }
            if ( (found_pivot_column >= N) ){
                fprintf(stderr,"row %d Pivot missing\n",row_idx);
                     rref_ok = 0;
                }
            if ( (is_pivot_column[found_pivot_column] != 1)){
                fprintf(stderr,"row %d indicator array mismatch\n",row_idx);
                     rref_ok = 0;
                }
        }

        if(full_rank && !rref_ok){
            fprintf(stderr,"RREF incorrect (iter:%d)\n",i);
            fprintf(stderr,"Pre-RREF\n");
            generator_pretty_print_name("Pre-RREF",&G);
            fprintf(stderr,"is_pivot = \n [ ");
            for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
            fprintf(stderr,"]\n");
            fprintf(stderr,"Post-RREF\n");
            generator_pretty_print_name("Post-RREF",&GMul);
            all_ok = 0;
        }

    }
    if(all_ok) {
        printf("GE test: ok\n");
    }
}

/* tests if G*M1*(M1^-1) = G*/
void gen_by_monom_tester(void){
     generator_mat_t G = {0}, G2, Gcheck;
     uint8_t is_pivot_column[NN];
     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);
     monomial_t mat1, mat2;
     monomial_mat_rnd(&mat1);
     monomial_mat_inv(&mat2,&mat1);
     generator_monomial_mul(&G2,&G,&mat1);
     generator_monomial_mul(&Gcheck,&G2,&mat2);
     if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        generator_pretty_print_name("G",&G);
        generator_pretty_print_name("G*Q",&G2);
        generator_pretty_print_name("G*Q*Q^-1",&Gcheck);
     } else {
         printf("Generator-monomial multiplication: ok\n");
     }
}

/* draw random full rank G and pack it
 * compute mu = Q_a^-1 Q_b
 * compute G2 = G Q_a
 * compute G3 = G Q_b
 * compute Gcheck = G2 mu
 */
void rref_gen_by_monom_tester(void){
     generator_mat_t G, G2, G3, Gcheck;

     monomial_t Q_a,Q_a_inv,Q_b,mu,Qcheck;
     monomial_mat_rnd(&Q_a);
     monomial_mat_inv(&Q_a_inv,&Q_a);
     monomial_mat_rnd(&Q_b);

     uint8_t is_pivot_column[NN];
     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     generator_monomial_mul(&G2,&G,&Q_a);
     if (generator_RREF(&G2,is_pivot_column) != 1){
         printf("G2=G Q_a: singular\n");

     };
     generator_monomial_mul(&G3,&G,&Q_b);
     if (generator_RREF(&G3,is_pivot_column) != 1){
         printf("G3=G Q_b: singular\n");

     };
     monomial_mat_mul(&mu,&Q_a_inv,&Q_b);


    monomial_mat_mul(&Qcheck,&Q_a,&mu);
    if( memcmp( &Q_b,&Qcheck,sizeof(monomial_t)) !=0 ){
        monomial_mat_pretty_print_name("mu",&mu);
        monomial_mat_pretty_print_name("Q_a",&Q_a);
        monomial_mat_pretty_print_name("Qcheck",&Qcheck);
        monomial_mat_pretty_print_name("Q_b",&Q_b);
        fprintf(stderr,"Q_a mu != Q_b\n");
    }


     generator_monomial_mul(&Gcheck,&G2,&mu);
     generator_RREF(&Gcheck,is_pivot_column);
     if (generator_RREF(&Gcheck,is_pivot_column) != 1){
         printf("Gcheck=G2 mu: singular\n");

     };

     if( memcmp( &Gcheck,&G3,sizeof(generator_mat_t)) !=0 ){
         printf("SF Generator-monomial multiplication: not ok\n");
        generator_pretty_print_name("G",&G);
        generator_pretty_print_name("G2",&G2);
        generator_pretty_print_name("G3",&G3);
        generator_pretty_print_name("Gcheck",&Gcheck);
        monomial_mat_print_exp_name("Q_a",&Q_a);
        monomial_mat_print_exp_name("Q_a_inv",&Q_a_inv);
        monomial_mat_print_exp_name("Q_b",&Q_b);
        monomial_mat_print_exp_name("mu",&mu);
        monomial_mat_print_exp_name("Qcheck",&Qcheck);
     } else {
         printf("SF Generator-monomial multiplication: ok\n");
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
    int is_signature_ok;
    is_signature_ok = LESS_verify(&pk,message,8,&signature);
    // fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
    return is_signature_ok;
}

#define NUMBER_OF_TESTS 100
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
        // fprintBstr(stdout, "pk = ", pk, CRYPTO_PUBLICKEYBYTES);
        // fprintBstr(stdout, "sk = ", sk, CRYPTO_SECRETKEYBYTES);

        if ( (ret_val = crypto_sign(sm, &smlen, m, msg_len, sk)) != 0) {
            printf("crypto_sign returned <%d>\n", ret_val);
            return -1;
        }

        //fprintBstr(stdout, "sm = ", sm, smlen);
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
    unsigned char m[] = {0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
	uint8_t seed[] = {0x06,0x15,0x50,0x23,0x4D,0x15,0x8C,0x5E,0xC9,0x55,0x95,0xFE,0x04,0xEF,0x7A,0x25,0x76,0x7F,0x2E,0x24,0xCC,0x2B,0xC4,0x79,0xD0,0x9D,0x86,0xDC,0x9A,0xBC,0xFD,0xE7,0x05,0x6A,0x8C,0x26,0x6F,0x9E,0xF9,0x7E,0xD0,0x85,0x41,0xDB,0xD2,0xE1,0xFF,0xA1};
    
	for (uint32_t i=0; i<48; i++) {
        seed[i] = i;
	}

    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)seed,
                      48);

    const uint32_t mlen = 33;
    unsigned long long smlen = 0, mlen1;
    //unsigned char *m  = (unsigned char *)calloc(mlen, sizeof(unsigned char));
    unsigned char *m1 = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char *sm = (unsigned char *)calloc(mlen+CRYPTO_BYTES, sizeof(unsigned char));
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES] = {0}, sk[CRYPTO_SECRETKEYBYTES] = {0};

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

    return 0;
}


#define NUM_TEST_ITERATIONS 1
int main(int argc, char* argv[]){
    (void)argc;
    (void)argv;
    return LESS_sign_verify_test_multiple();
    //LESS_sign_verify_test_KAT();

    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"012345678912345",
                      16);
    fprintf(stderr,"LESS reference implementation functional testbench\n");
    info();

    int tests_ok = 0;
    for (int i = 0; i < NUM_TEST_ITERATIONS; i++) {
        fputc('.',stderr);
      // fprintf(stderr,"test %d: ",i);
       inverse_mod_tester();
       gen_by_monom_tester();
       monomial_tester();
       rref_gen_compress_tester();
       gausselim_tester();
       rref_gen_by_monom_tester();
        // rref_gen_byte_compress_tester();
        // mono_is_compress_tester();
      tests_ok += LESS_sign_verify_test();
    }
    fprintf(stderr,"%d tests functional out of %d\n",tests_ok,NUM_TEST_ITERATIONS);    
    return 0;
}
