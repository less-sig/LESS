#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fq_arith.h"
#include "monomial_mat.h"
#include "timing_and_stat.h"
#include "codes.h"
#include "LESS.h"
#include "rng.h"

#define GRN "\e[0;32m"
#define WHT "\e[0;37m"

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif

/* Exhaustive testing of inverses mod Q */
void inverse_mod_tester(){
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
void monomial_tester(){
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
void gausselim_tester(){
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
void gen_by_monom_tester(){
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


/* draw random full rank G and pack it */
/* compute mu = Q_a^-1 Q_b */
/* compute G2 = G Q_a */
/* compute G3 = G Q_b */
/* compute Gcheck = G2 mu */

void rref_gen_by_monom_tester(){
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

void rref_gen_compress_tester(){
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

void mono_is_compress_tester(){
    monomial_action_IS_t Q_a, Qcheck;
    uint8_t compressed [MONO_ACTION_PACKEDBYTES];

    monomial_t mono_rnd;
    monomial_mat_rnd(&mono_rnd);

    // Create random q
    for (int i = 0; i < K; i++) {
        Q_a.coefficients[i] = mono_rnd.coefficients[i];
        Q_a.permutation[i] = mono_rnd.permutation[i];
    }

     compress_monom_action(compressed,&Q_a);
     expand_to_monom_action(&Qcheck,compressed);

    if( memcmp( &Qcheck,&Q_a,sizeof(monomial_action_IS_t)) !=0 ){
        printf("Monomial Action compression: ko\n");

       fprintf(stderr,"perm = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Q_a.permutation[i]);
       }
       fprintf(stderr,"%03u ]\n",Q_a.permutation[K-1]);
       fprintf(stderr,"coeffs = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Q_a.coefficients[i]);
       }
       fprintf(stderr,"%03u ]\n",Q_a.coefficients[K-1]);

       fprintf(stderr,"\n\n\n");
       fprintf(stderr,"perm = [");
       for(int i = 0; i < K-1; i++
        ) {
          fprintf(stderr,"%03u, ",Qcheck.permutation[i]);
       }
       fprintf(stderr,"%03u ]\n",Qcheck.permutation[K-1]);
       fprintf(stderr,"coeffs = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Qcheck.coefficients[i]);
       }
       fprintf(stderr,"%03u ]\n",Qcheck.coefficients[K-1]);

    } else {
        printf("Monomial Action compression: ok\n");
    }

}


void rref_gen_byte_compress_tester(){
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
     expand_to_rref(&G,G_compressed);

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

void info(){
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB\n", sizeof(pubkey_t));
    fprintf(stderr,"Signature: %luB\n", sizeof(sig_t));

}

/* returns 1 if the test is successful, 0 otherwise */
int LESS_sign_verify_test(){
    pubkey_t pk;
    prikey_t sk;
    sig_t signature;
    char message[8] = "Signme!";
    LESS_keygen(&sk,&pk);
    LESS_sign(&sk,message,8,&signature);
    int is_signature_ok;
    is_signature_ok = LESS_verify(&pk,message,8,&signature);
    // fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
    return is_signature_ok;
}

#define NUM_TEST_ITERATIONS 20
int main(int argc, char* argv[]){
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"012345678912345",
                      16);
    fprintf(stderr,"LESS reference implementation functional testbench\n");
    info();
    int tests_ok = 0;
    for (int i = 0; i < NUM_TEST_ITERATIONS; i++) {
        fputc('.',stderr);
      // fprintf(stderr,"test %d: ",i);
//       inverse_mod_tester();
//       gen_by_monom_tester();
//       monomial_tester();
//       monomial_IS_tester();
//       rref_gen_compress_tester();
//       gausselim_tester();
//       rref_gen_by_monom_tester();
        // rref_gen_byte_compress_tester();
        // mono_is_compress_tester();
      tests_ok += LESS_sign_verify_test();
    }
    fprintf(stderr,"%d tests functional out of %d\n",tests_ok,NUM_TEST_ITERATIONS);    
    return 0;
}
