#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"
#include "test_helpers.c"
#include "cf_test.h"

#define ITERS 100


// just debugging test
int test_compute_canonical_form_type3(void) {
    normalized_IS_t G;
    normalized_sf(&G);
    uint8_t L[N] = {0};

    if (compute_canonical_form_type3(&G, L) == 0) return 1;
    return 0;
}

// real test
int test_compute_canonical_form_type3_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r;
    uint8_t L[N] = {0};

    for (uint32_t k = 0; k < ITERS; k++) {
        permutation_mat_rng_v2(&P_c, N-K);
        permutation_mat_rng_v2(&P_r, K);

        // generate data
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);
        permutation_apply_col(&G2, &P_c);
        permutation_apply_row(&P_r, &G2);
        normalized_copy(&G3, &G2);

        if (compute_canonical_form_type3(&G1, L) == 0) return 1;
        if (compute_canonical_form_type3(&G2, L) == 0) return 1;
        if (compute_canonical_form_type3_ct(&G3) == 0) return 1;


        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G2);
                    printf("error type3\n");
                    return 1;
                }

                if (G1.values[i][j] != G3.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G3);
                    printf("error type3 ct\n");
                    return 1;
                }
            }
        }
    }

    return 0;
}

int test_compute_canonical_form_type4(void) {
    normalized_IS_t G;
    normalized_sf(&G);
    uint8_t L[N] = {0};

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type4(&G, L) == 0) return 1;
    // normalized_pretty_print(&G);

    return 0;
}

int test_compute_canonical_form_type4_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r;
    diagonal_t D_r;
    uint8_t L[N] = {0};

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
        const int ret1 = compute_canonical_form_type4(&G1, L);
        const int ret2 = compute_canonical_form_type4(&G2, L);
        const int ret3 = compute_canonical_form_type4_ct(&G3);

        if (ret1 != ret2) {
            printf("error ret cf4 2\n");
            return 1;
        }

        if (ret1 != ret3) {
            printf("error ret cf4 3\n");
            return 1;
        }
        if (ret1 == 1) {
            for (uint32_t i = 0; i < K; i++) {
                for (uint32_t j = 0; j < N-K; j++) {
                    if (G1.values[i][j] != G2.values[i][j]) {
                        printf("error cf4 2\n");
                        normalized_pretty_print(&G1);
                        normalized_pretty_print(&G2);
                        return 1;
                    }

                    if (G3.values[i][j] != G2.values[i][j]) {
                        printf("error cf4 3\n");
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

int test_compute_canonical_form_type5(void) {
    normalized_IS_t G0, G1, G2, G3;
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < K; j++) {
            G0.values[i][j] = cf_out[i][j];
            G1.values[i][j] = cf_in[i][j];
            G2.values[i][j] = cf_in[i][j];
            G3.values[i][j] = cf_in[i][j];
        }
    }

    if (compute_canonical_form_type5(&G1) == 0)        { return 1; }
    if (compute_canonical_form_type5_popcnt(&G2) == 0) { return 1; }
    if (compute_canonical_form_type5_tony(&G3) == 0) { return 1; }

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G1.values[i][j] != cf_out[i][j]) {
                normalized_pretty_print(&G0);
                normalized_pretty_print(&G1);
                printf("5 1\n");
                return 1;
            }
            if (G2.values[i][j] != cf_out[i][j]) {
                normalized_pretty_print(&G0);
                normalized_pretty_print(&G2);
                printf("5 2: %d %d\n", i, j);
                return 1;
            }
            //if (G3.values[i][j] != cf_out[i][j]) {
            //    normalized_pretty_print(&G0);
            //    normalized_pretty_print(&G2);
            //    printf("5 3\n");
            //    return 1;
            //}
            // if (G4.values[i][j] != cf_out[i][j]) {
            //     printf("5 4\n");
            //     return 1;
            // }
        }
    }

    return 0;
}

int test_compute_canonical_form_type5_v2(void) {
    normalized_IS_t G3, G1, G2;
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

        const int ret1 = compute_canonical_form_type5(&G1);
        const int ret2 = compute_canonical_form_type5_popcnt(&G2);
        const int ret3 = compute_canonical_form_type5_tony(&G3);
        if (ret1 != ret2) {
            printf("error ret cf5 popcnt\n");
            return 1;
        }

        if (ret1 != ret3) {
            printf("error ret cf5, tony\n");
            return 1;
        }

        if (ret1 == 0) {
            printf("THIS SHOULD NOT HAPPEN\n");
            continue;
        }

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G2);
                    printf("error cf5 popcnt %d %d\n", i, j);
                    return 1;
                }
            }
        }

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G3.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G3);
                    printf("error cf5 %d tony %d %d\n", k, i, j);
                    return 1;
                }
            }
        }
    }
    return 0;
}

// tests: if CF(G) and CF(RREF(G)) are the same
uint32_t  test_compute_canonical_form_type5_gaus(void) {
    generator_mat_t G1, G2;
    normalized_IS_t V1, V2, V3;
    for (uint32_t k = 0; k < ITERS; ++k) {
        generator_sf(&G1);
        memcpy(&G2, &G1, sizeof(generator_mat_t));

        uint8_t is_pivot_column[N] = {0};
        generator_RREF(&G2, is_pivot_column);
        for (uint32_t i = 0; i < K; ++i) {
            for (uint32_t j = 0; j < N; ++j) {
                if (G2.values[i][j] >= Q) {
                    return 2;
                }
            }
        }


        generator_to_normalized(&V1, &G1);
        generator_to_normalized(&V2, &G2);
        generator_to_normalized(&V3, &G2);

        const int ret1 = compute_canonical_form_type5(&V1);
        const int ret2 = compute_canonical_form_type5_popcnt(&V2);
        const int ret3 = compute_canonical_form_type5_tony(&V3);

        if ((ret1 != ret2) && (ret1 != ret3)) {
            printf("error: cf5 gaus ret\n");
            return 1;
        }

        for (uint32_t i = 0; i < K; ++i) {
            for (uint32_t j = 0; j < N - K; ++j) {
                if (V1.values[i][j] != V2.values[i][j]) {
                    normalized_pretty_print(&V1);
                    normalized_pretty_print(&V2);
                    printf("error: cf5 popcnt gaus %d %d\n", i, j);
                    return 1;
                }

                if (V1.values[i][j] != V3.values[i][j]) {
                    normalized_pretty_print(&V1);
                    normalized_pretty_print(&V3);
                    printf("error: cf5 tony gaus %d %d \n", i, j);
                    return 1;
                }
            }
        }
    }

    return 0;
}

int main(void) {
    // generic matrices test, not really testing anything, just printing the result
    // if (test_compute_canonical_form_type3()) return 1;
    // if (test_compute_canonical_form_type4()) return 1;

    // actual tests
    // if (test_compute_canonical_form_type3_v2()) return 1;
    // if (test_compute_canonical_form_type4_v2()) return 1;
    // if (test_compute_canonical_form_type5()) return 1;
    if (test_compute_canonical_form_type5_v2()) return 1;
    if (test_compute_canonical_form_type5_gaus()) return 1;

    printf("Done, all worked\n");
    return 0;
}
