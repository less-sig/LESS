#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"


int test_compute_canonical_form_type2(void) {
    normalized_IS_t G;
    diagonal_t D_c;
    permutation_t P_c;

    permutation_mat_id(&P_c);
    diagonal_mat_id(&D_c);
    normalized_sf(&G);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type2(&G, &D_c, &P_c) == 0) return 1;
    // normalized_pretty_print(&G);

    return 0;
}

// just debugging test
int test_compute_canonical_form_type3(void) {
    normalized_IS_t G;
    permutation_t P_c, P_r, P_c_ret, P_r_ret;

    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);

    normalized_sf(&G);
    permutation_mat_id(&P_c_ret);
    permutation_mat_id(&P_r_ret);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type3(&G, &P_r_ret, &P_c_ret) == 0) return 1;
    // normalized_pretty_print(&G);

    return 0;
}

// real test
int test_compute_canonical_form_type3_v2(void) {
    normalized_IS_t G1, G2;
    permutation_t P_c, P_r, P_c1, P_r1, P_c2, P_r2;

    permutation_mat_rng_v2(&P_c, N-K);
    permutation_mat_rng_v2(&P_r, K);

    // generate data
    normalized_sf(&G1);
    normalized_copy(&G2, &G1);
    permutation_apply_col(&G2, &P_c);
    permutation_apply_row(&P_r, &G2);

    permutation_mat_id(&P_c1);
    permutation_mat_id(&P_r1);
    permutation_mat_id(&P_c2);
    permutation_mat_id(&P_r2);

    if (compute_canonical_form_type3(&G1, &P_r1, &P_c1) == 0) return 1;
    if (compute_canonical_form_type3(&G2, &P_r2, &P_c2) == 0) return 1;

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

    return 0;
}

int test_compute_canonical_form_type4(void) {
    normalized_IS_t G;
    permutation_t P_c, P_r;
    diagonal_t D_c;

    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);
    diagonal_mat_id(&D_c);

    normalized_sf(&G);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type4(&G, &P_r, &D_c, &P_c) == 0) return 1;
    // normalized_pretty_print(&G);

    return 0;
}

int test_compute_canonical_form_type4_v2(void) {
    normalized_IS_t G1, G2, G3;
    permutation_t P_c, P_r, P_c1, P_r1, P_c2, P_r2, P_c3, P_r3;
    diagonal_t D_r, D_r1, D_r2, D_r3;

    diagonal_mat_rnd_v2(&D_r, K);
    permutation_mat_rng_v2(&P_c, N-K);
    permutation_mat_rng_v2(&P_r, K);

    // generate data
    normalized_sf(&G1);
    normalized_copy(&G2, &G1);
    permutation_apply_col(&G2, &P_c);
    diagonal_apply_row(&D_r, &G2);
    permutation_apply_row(&P_r, &G2);

    memcpy((void *)G3.values, (void*)G2.values, sizeof(normalized_IS_t));
    permutation_mat_id(&P_c1);
    permutation_mat_id(&P_r1);
    permutation_mat_id(&P_c2);
    permutation_mat_id(&P_r2);
    permutation_mat_id(&P_c3);
    permutation_mat_id(&P_r3);

    if (compute_canonical_form_type4(&G1, &P_r1, &D_r1, &P_c1) == 0) return 1;
    if (compute_canonical_form_type4(&G2, &P_r2, &D_r2, &P_c2) == 0) return 1;
    if (compute_canonical_form_type4_non_ct(&G3, &P_r3, &D_r3, &P_c3) == 0) return 1;

    // normalized_pretty_print(&G1);
    // normalized_pretty_print(&G2);

    // diagonal_pretty_print(&D_r2);
    // diagonal_pretty_print(&D_r3);

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G1.values[i][j] != G2.values[i][j]) {
                printf("error cf4 2\n");
                return 1;
            }

            if (G3.values[i][j] != G2.values[i][j]) {
                printf("error cf4 3\n");
                return 1;
            }
        }
    }

    for (uint32_t i = 0; i < N-K; i++) {
        if (P_c2.permutation[i] != P_c3.permutation[i]) {
            printf("error cf4 p_c\n");
            return 2;
        }
    }

    for (uint32_t i = 0; i < K; i++) {
        if (P_r2.permutation[i] != P_r3.permutation[i]) {
            printf("error cf4 p_3\n");
            return 3;
        }
    }

    for (uint32_t i = 0; i < K; i++) {
        if (D_r2.coefficients[i] != D_r3.coefficients[i]) {
            printf("error cf4 D\n");
            return 3;
        }
    }

    return 0;
}

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

    if (compute_canonical_form_type4(&G1, NULL, NULL, NULL) == 0) return 1;

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

int test_compute_canonical_form_type5(void) {
    normalized_IS_t G;
    permutation_t P_c, P_r;
    diagonal_t D_c, D_r;

    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);
    diagonal_mat_id(&D_c);
    diagonal_mat_id(&D_r);

    normalized_sf(&G);

    // normalized_pretty_print(&G);
    if (compute_canonical_form_type5(&G, &D_r, &P_r, &D_c, &P_c) == 0) return 1;
    // normalized_pretty_print(&G);
    return 0;
}

int test_compute_canonical_form_type5_v2(void) {
    normalized_IS_t G1, G2;
    permutation_t P_c, P_r, P_c1, P_r1, P_c2, P_r2;
    diagonal_t D_r, D_c, D_r1, D_c1, D_r2, D_c2;

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

    permutation_mat_id(&P_c1);
    permutation_mat_id(&P_r1);
    permutation_mat_id(&P_c2);
    permutation_mat_id(&P_r2);

    if (compute_canonical_form_type5(&G1, &D_r1, &P_r1, &D_c1, &P_c1) == 0) return 1;
    if (compute_canonical_form_type5(&G2, &D_r2, &P_r2, &D_c2, &P_c2) == 0) return 1;

    // normalized_pretty_print(&G1);
    // normalized_pretty_print(&G2);

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G1.values[i][j] != G2.values[i][j]) {
                printf("error cf5\n");
                return 1;
            }
        }
    }

    return 0;
}

int test_compute_canonical_form_type5_v3(void) {
    normalized_IS_t G1;

    /// just generate some deterministic numbers
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            G1.values[i][j] = 120 - i - j;
        }
    }

    // normalized_pretty_print(&G1);
    if (compute_canonical_form_type5(&G1, NULL, NULL, NULL, NULL) == 0) return 1;
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

int main(void) {
    // generic matrices test, not really testing anything, just printing the result
    // if (test_compute_canonical_form_type3()) return 1;
    // if (test_compute_canonical_form_type4()) return 1;
    // if (test_compute_canonical_form_type5()) return 1;

    // actual tests
    if (test_compute_canonical_form_type3_v2()) return 1;
    if (test_compute_canonical_form_type4_v2()) return 1;
    if (test_compute_canonical_form_type5_v2()) return 1;

    // value tests, (values taken from cf.py)
    if (test_compute_canonical_form_type4_v3()) return 1;
    if (test_compute_canonical_form_type5_v3()) return 1;

    return 0;
}
