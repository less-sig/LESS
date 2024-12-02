#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"

#define ITERS 100

// TODO write tests for cf4 and cf5, which lead to an wrong computation

// just debugging test
int test_compute_canonical_form_type3(void) {
    normalized_IS_t G;
    normalized_sf(&G);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type3(&G) == 0) return 1;
    // normalized_pretty_print(&G);

    return 0;
}

// real test
int test_compute_canonical_form_type3_v2(void) {
    normalized_IS_t G1, G2;
    permutation_t P_c, P_r;

    for (uint32_t k = 0; k < ITERS; k++) {
        permutation_mat_rng_v2(&P_c, N-K);
        permutation_mat_rng_v2(&P_r, K);

        // generate data
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);
        permutation_apply_col(&G2, &P_c);
        permutation_apply_row(&P_r, &G2);

        if (compute_canonical_form_type3(&G1) == 0) return 1;
        if (compute_canonical_form_type3(&G2) == 0) return 1;

        // normalized_pretty_print(&G1);
        // normalized_pretty_print(&G2);

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    printf("error type3\n");
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

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type4(&G) == 0) return 1;
    // normalized_pretty_print(&G);

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
        const int ret3 = compute_canonical_form_type4_non_ct(&G3);

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

#if defined(CATEGORY_0)
// tests dummy data
int test_compute_canonical_form_type4_v3(void) {
    if (K != 8) {
        return 0; 
    }
    normalized_IS_t G1;

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            G1.values[i][j] = 120 - i - j;
        }
    }

    if (compute_canonical_form_type4(&G1) == 0) return 1;

    // data comes from cf.py
    const uint8_t data[K][N-K] = {
        {  3,  4, 28, 29, 54, 79, 80,105},
        {113, 86, 73, 46,  6, 93, 66, 26},
        { 36, 54,105,123, 65,  7, 25, 94},
        {  9,117, 42, 23, 56, 89, 70,103},
        { 45, 33,126,114, 68, 22, 10, 91},
        { 11, 70, 89, 21, 99, 50,109, 60},
        { 80, 36,123, 79,122, 38,121, 37},
        { 50,106, 53,109,112,115, 44, 47},
    };

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (data[i][j] != G1.values[i][j]) {
                printf("error: test_compute_canonical_form_type4_v3\n");
                return 1;
            }
        }
    }

    return 0;
}
#endif

int test_compute_canonical_form_type5(void) {
    normalized_IS_t G;
    normalized_sf(&G);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type5(&G) == 0) return 1;
    // normalized_pretty_print(&G);
    return 0;
}

int test_compute_canonical_form_type5_v2(void) {
    normalized_IS_t G1, G2;
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

        const int ret1 = compute_canonical_form_type5(&G1);
        const int ret2 = compute_canonical_form_type5(&G2);
        if (ret1 != ret2) {
            printf("error ret cf5\n");
            return 1;
        }

        // normalized_pretty_print(&G1);
        // normalized_pretty_print(&G2);

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G2);
                    printf("error cf5\n");
                    return 1;
                }
            }
        }
    }
    return 0;
}

#if defined(CATEGORY_0)
int test_compute_canonical_form_type5_v3(void) {
    normalized_IS_t G1;

    /// just generate some deterministic numbers
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            G1.values[i][j] = 120 - i - j;
        }
    }

    // normalized_pretty_print(&G1);
    if (compute_canonical_form_type5(&G1) == 0) return 1;
    // normalized_pretty_print(&G1);
    // data comes from cf.py
    const uint8_t data[K][N-K] = {
        { 1,  1,  1,  1,  1,  1,  1,  1},
        { 1,  2, 17, 51, 53, 80, 85, 93},
        {84, 71,  3, 69, 43, 73,  8, 31},
        {46, 44, 14, 73, 69, 15,  5,116},
        {83, 87, 20, 29, 37, 18, 38, 70},
        {18,111,109, 96, 28,126, 83, 65},
        {28, 78, 66,115, 88, 41, 37, 56},
        {47, 28,124,113, 75, 70,102, 77},
    };

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (data[i][j] != G1.values[i][j]) {
                printf("error: test_compute_canonical_form_type5_v3\n");
                return 1;
            }
        }
    }

    return 0;
}
#endif

int main(void) {
    // generic matrices test, not really testing anything, just printing the result
    // if (test_compute_canonical_form_type3()) return 1;
    // if (test_compute_canonical_form_type4()) return 1;
    // if (test_compute_canonical_form_type5()) return 1;

    // actual tests
    if (test_compute_canonical_form_type3_v2()) return 1;
    if (test_compute_canonical_form_type4_v2()) return 1;
    if (test_compute_canonical_form_type5_v2()) return 1;

#if defined(CATEGORY_0)
    // value tests, (values taken from cf.py)
    if (test_compute_canonical_form_type4_v3()) return 1;
    if (test_compute_canonical_form_type5_v3()) return 1;
#endif

    return 0;
}
