#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"

#define ITERS 1

// real test
int test_compute_canonical_form_type3_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r;

    for (uint32_t k = 0; k < ITERS; k++) {
        permutation_mat_rng_v2(&P_c, N-K);
        permutation_mat_rng_v2(&P_r, K);

        // generate data
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);
        permutation_apply_col(&G2, &P_c);
        permutation_apply_row(&P_r, &G2);
        normalized_copy(&G3, &G2);

        if (compute_canonical_form_type3(&G1) == 0) return 1;
        if (compute_canonical_form_type3(&G2) == 0) return 1;
        if (compute_canonical_form_type3_ct(&G3) == 0) return 1;


        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    // normalized_pretty_print(&G1);
                    // normalized_pretty_print(&G2);
                    printf("error type3\n");
                    return 1;
                }

                if (G1.values[i][j] != G3.values[i][j]) {
                    // normalized_pretty_print(&G1);
                    // normalized_pretty_print(&G3);
                    printf("error type3 ct\n");
                    return 1;
                }
            }
        }
    }

    return 0;
}

int test_compute_canonical_form_type4_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r;
    diagonal_t D_r;

    for (uint32_t k = 0; k < ITERS; k++) {
        diagonal_mat_rnd_v2(&D_r, K);
        permutation_mat_rng_v2(&P_c, N-K);
        permutation_mat_rng_v2(&P_r, K);

        // generate data
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);
        permutation_apply_col(&G2, &P_c);
        diagonal_apply_row(&D_r, &G2);
        permutation_apply_row(&P_r, &G2);

        normalized_copy(&G3, &G2);
        const int ret1 = compute_canonical_form_type4(&G1);
        const int ret2 = compute_canonical_form_type4(&G2);
        // TODO
        // const int ret3 = compute_canonical_form_type4_ct(&G3);

        if (ret1 != ret2) {
            printf("error ret cf4 2\n");
            return 1;
        }

        // if (ret1 != ret3) {
        //     printf("error ret cf4 3\n");
        //     return 1;
        // }
        if (ret1 == 1) {
            for (uint32_t i = 0; i < K; i++) {
                for (uint32_t j = 0; j < N-K; j++) {
                    if (G1.values[i][j] != G2.values[i][j]) {
                        printf("error cf4 2\n");
                        // normalized_pretty_print(&G1);
                        // normalized_pretty_print(&G2);
                        return 1;
                    }

                    // if (G3.values[i][j] != G2.values[i][j]) {
                    //     printf("error cf4 3\n");
                    //     return 1;
                    // }
                }
            }
        }
    }
    return 0;
}

int test_compute_canonical_form_type5_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r;
    diagonal_t D_r, D_c;

    for (uint32_t k = 0; k < ITERS; k++) {
        permutation_mat_rng_v2(&P_c, N-K);
        permutation_mat_rng_v2(&P_r, K);
        diagonal_mat_rnd_v2(&D_r, K);
        diagonal_mat_rnd_v2(&D_c, N-K);

        // generate data
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);
        diagonal_apply_col(&G2, &D_c);
        diagonal_apply_row(&D_r, &G2);
        permutation_apply_col(&G2, &P_c);
        permutation_apply_row(&P_r, &G2);
        normalized_copy(&G3, &G2);

        const int ret1 = compute_canonical_form_type5_ct(&G1);
        const int ret2 = compute_canonical_form_type5(&G2);
        const int ret3 = compute_canonical_form_type5_popcnt(&G3);
        if (ret1 != ret2) {
            printf("error ret cf5\n");
            return 1;
        }
        
        if (ret1 != ret3) {
            printf("error ret cf5, popcnt\n");
            return 1;
        }
        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    // normalized_pretty_print(&G1);
                    // normalized_pretty_print(&G2);
                    printf("error cf5\n");
                    return 1;
                }
                
                if (G2.values[i][j] != G3.values[i][j]) {
                    // normalized_pretty_print(&G1);
                    // normalized_pretty_print(&G2);
                    // normalized_pretty_print(&G3);
                    printf("error cf5 popcnt %d %d\n", i, j);
                    return 1;
                }
            }
        }
    }
    return 0;
}


int main(void) {
    // actual tests
    if (test_compute_canonical_form_type3_v2()) return 1;
    if (test_compute_canonical_form_type4_v2()) return 1;
    if (test_compute_canonical_form_type5_v2()) return 1;

    printf("Done, all worked\n");
    return 0;
}
