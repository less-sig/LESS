#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"
#include "test_helpers.c"
#include "cf_test.h"
#include "../bench/cycles.h"

#define ITERS 100
#define BENCH_ITERS 1000


// just debugging test
// int test_compute_canonical_form_type3(void) {
//     normalized_IS_t G;
//     normalized_sf(&G);
//     uint8_t L[N] = {0};
// 
//     if (compute_canonical_form_type3(&G, L) == 0) return 1;
//     return 0;
// }
// 
// // real test
// int test_compute_canonical_form_type3_v2(void) {
//     normalized_IS_t G1, G2, G3;
//     permutation_t P_c, P_r;
//     uint8_t L[N] = {0};
// 
//     for (uint32_t k = 0; k < ITERS; k++) {
//         permutation_mat_rng_v2(&P_c, N-K);
//         permutation_mat_rng_v2(&P_r, K);
// 
//         // generate data
//         normalized_sf(&G1);
//         normalized_copy(&G2, &G1);
//         permutation_apply_col(&G2, &P_c);
//         permutation_apply_row(&P_r, &G2);
//         normalized_copy(&G3, &G2);
// 
//         if (compute_canonical_form_type3(&G1, L) == 0) return 1;
//         if (compute_canonical_form_type3(&G2, L) == 0) return 1;
//         if (compute_canonical_form_type3_ct(&G3) == 0) return 1;
// 
// 
//         for (uint32_t i = 0; i < K; i++) {
//             for (uint32_t j = 0; j < N-K; j++) {
//                 if (G1.values[i][j] != G2.values[i][j]) {
//                     normalized_pretty_print(&G1);
//                     normalized_pretty_print(&G2);
//                     printf("error type3\n");
//                     return 1;
//                 }
// 
//                 if (G1.values[i][j] != G3.values[i][j]) {
//                     normalized_pretty_print(&G1);
//                     normalized_pretty_print(&G3);
//                     printf("error type3 ct\n");
//                     return 1;
//                 }
//             }
//         }
//     }
// 
//     return 0;
// }
// 
// int test_compute_canonical_form_type4(void) {
//     normalized_IS_t G;
//     normalized_sf(&G);
//     uint8_t L[N] = {0};
// 
//     // normalized_pretty_print(&G);
//     if (compute_canonical_form_type4(&G, L) == 0) return 1;
//     // normalized_pretty_print(&G);
// 
//     return 0;
// }
// 
// int test_compute_canonical_form_type4_v2(void) {
//     normalized_IS_t G1, G2, G3;
//     permutation_t P_c, P_r;
//     diagonal_t D_r;
//     uint8_t L[N] = {0};
// 
//     for (uint32_t k = 0; k < ITERS; k++) {
//         diagonal_mat_rnd_v2(&D_r, K);
//         permutation_mat_rng_v2(&P_c, N-K);
//         permutation_mat_rng_v2(&P_r, K);
// 
//         // generate data
//         normalized_sf(&G1);
//         normalized_copy(&G2, &G1);
//         permutation_apply_col(&G2, &P_c);
//         diagonal_apply_row(&D_r, &G2);
//         permutation_apply_row(&P_r, &G2);
// 
//         normalized_copy(&G3, &G2);
//         const int ret1 = compute_canonical_form_type4(&G1, L);
//         const int ret2 = compute_canonical_form_type4(&G2, L);
//         const int ret3 = compute_canonical_form_type4_ct(&G3);
// 
//         if (ret1 != ret2) {
//             printf("error ret cf4 2\n");
//             return 1;
//         }
// 
//         if (ret1 != ret3) {
//             printf("error ret cf4 3\n");
//             return 1;
//         }
//         if (ret1 == 1) {
//             for (uint32_t i = 0; i < K; i++) {
//                 for (uint32_t j = 0; j < N-K; j++) {
//                     if (G1.values[i][j] != G2.values[i][j]) {
//                         printf("error cf4 2\n");
//                         normalized_pretty_print(&G1);
//                         normalized_pretty_print(&G2);
//                         return 1;
//                     }
// 
//                     if (G3.values[i][j] != G2.values[i][j]) {
//                         printf("error cf4 3\n");
//                         return 1;
//                     }
//                 }
//             }
//         }
//     }
//     return 0;
// }
// 
// int test_compute_canonical_form_type5(void) {
//     normalized_IS_t G0, G1, G2, G3;
//     for (uint32_t i = 0; i < K; i++) {
//         for (uint32_t j = 0; j < K; j++) {
//             G0.values[i][j] = cf_out[i][j];
//             G1.values[i][j] = cf_in[i][j];
//             G2.values[i][j] = cf_in[i][j];
//             G3.values[i][j] = cf_in[i][j];
//         }
//     }
// 
//     if (compute_canonical_form_type5(&G1) == 0)        { return 1; }
//     if (compute_canonical_form_type5_popcnt(&G2) == 0) { return 1; }
//     if (compute_canonical_form_type5_tony(&G3) == 0) { return 1; }
// 
//     for (uint32_t i = 0; i < K; i++) {
//         for (uint32_t j = 0; j < N-K; j++) {
//             if (G1.values[i][j] != cf_out[i][j]) {
//                 normalized_pretty_print(&G0);
//                 normalized_pretty_print(&G1);
//                 printf("5 1\n");
//                 return 1;
//             }
//             if (G2.values[i][j] != cf_out[i][j]) {
//                 normalized_pretty_print(&G0);
//                 normalized_pretty_print(&G2);
//                 printf("5 2: %d %d\n", i, j);
//                 return 1;
//             }
//             //if (G3.values[i][j] != cf_out[i][j]) {
//             //    normalized_pretty_print(&G0);
//             //    normalized_pretty_print(&G2);
//             //    printf("5 3\n");
//             //    return 1;
//             //}
//             // if (G4.values[i][j] != cf_out[i][j]) {
//             //     printf("5 4\n");
//             //     return 1;
//             // }
//         }
//     }
// 
//     return 0;
// }
// 
// int test_compute_canonical_form_type5_v2(void) {
//     normalized_IS_t G3, G1, G2;
//     permutation_t P_c, P_r;
//     diagonal_t D_r, D_c;
// 
//     for (uint32_t k = 0; k < ITERS; k++) {
//         permutation_mat_rng_v2(&P_c, N-K);
//         permutation_mat_rng_v2(&P_r, K);
//         diagonal_mat_rnd_v2(&D_r, K);
//         diagonal_mat_rnd_v2(&D_c, N-K);
// 
//         // generate data
//         normalized_sf(&G1);
//         normalized_copy(&G2, &G1);
//         diagonal_apply_col(&G2, &D_c);
//         diagonal_apply_row(&D_r, &G2);
//         permutation_apply_col(&G2, &P_c);
//         permutation_apply_row(&P_r, &G2);
// 
//         normalized_copy(&G3, &G2);
// 
//         const int ret1 = compute_canonical_form_type5(&G1);
//         const int ret2 = compute_canonical_form_type5_popcnt(&G2);
//         const int ret3 = compute_canonical_form_type5_tony(&G3);
//         if (ret1 != ret2) {
//             printf("error ret cf5 popcnt\n");
//             return 1;
//         }
// 
//         if (ret1 != ret3) {
//             printf("error ret cf5, tony\n");
//             return 1;
//         }
// 
//         if (ret1 == 0) {
//             printf("THIS SHOULD NOT HAPPEN\n");
//             continue;
//         }
// 
//         for (uint32_t i = 0; i < K; i++) {
//             for (uint32_t j = 0; j < N-K; j++) {
//                 if (G1.values[i][j] != G2.values[i][j]) {
//                     normalized_pretty_print(&G1);
//                     normalized_pretty_print(&G2);
//                     printf("error cf5 popcnt %d %d\n", i, j);
//                     return 1;
//                 }
//             }
//         }
// 
//         for (uint32_t i = 0; i < K; i++) {
//             for (uint32_t j = 0; j < N-K; j++) {
//                 if (G1.values[i][j] != G3.values[i][j]) {
//                     normalized_pretty_print(&G1);
//                     normalized_pretty_print(&G3);
//                     printf("error cf5 %d tony %d %d\n", k, i, j);
//                     return 1;
//                 }
//             }
//         }
//     }
//     return 0;
// }
// 
// // tests: if CF(G) and CF(RREF(G)) are the same
// uint32_t  test_compute_canonical_form_type5_gaus(void) {
//     generator_mat_t G1, G2;
//     normalized_IS_t V1, V2, V3;
//     for (uint32_t k = 0; k < ITERS; ++k) {
//         generator_sf(&G1);
//         memcpy(&G2, &G1, sizeof(generator_mat_t));
// 
//         uint8_t is_pivot_column[N] = {0};
//         generator_RREF(&G2, is_pivot_column);
//         for (uint32_t i = 0; i < K; ++i) {
//             for (uint32_t j = 0; j < N; ++j) {
//                 if (G2.values[i][j] >= Q) {
//                     return 2;
//                 }
//             }
//         }
// 
// 
//         generator_to_normalized(&V1, &G1);
//         generator_to_normalized(&V2, &G2);
//         generator_to_normalized(&V3, &G2);
// 
//         const int ret1 = compute_canonical_form_type5(&V1);
//         const int ret2 = compute_canonical_form_type5_popcnt(&V2);
//         const int ret3 = compute_canonical_form_type5_tony(&V3);
// 
//         if ((ret1 != ret2) && (ret1 != ret3)) {
//             printf("error: cf5 gaus ret\n");
//             return 1;
//         }
// 
//         for (uint32_t i = 0; i < K; ++i) {
//             for (uint32_t j = 0; j < N - K; ++j) {
//                 if (V1.values[i][j] != V2.values[i][j]) {
//                     normalized_pretty_print(&V1);
//                     normalized_pretty_print(&V2);
//                     printf("error: cf5 popcnt gaus %d %d\n", i, j);
//                     return 1;
//                 }
// 
//                 if (V1.values[i][j] != V3.values[i][j]) {
//                     normalized_pretty_print(&V1);
//                     normalized_pretty_print(&V3);
//                     printf("error: cf5 tony gaus %d %d \n", i, j);
//                     return 1;
//                 }
//             }
//         }
//     }
// 
//     return 0;
// }

/// simple dot product
static inline
FQ_ELEM row_dot(const FQ_ELEM *row1, const FQ_ELEM *row2) {
    uint32_t acc = 0;
    for (uint32_t col = 0; col < N-K; col++) {
        acc += (uint32_t)row1[col] * row2[col];
    }

    return acc % Q;
}

/// length of the tuple psi(i)
#define EDO_D 3

/// CF'' from cf5pp.pdf, Algorithm 1
/// \return 0 on failure (no admissible row, tie in psi, or identical rows)
/// 		1 on success
int EdoardoCF(normalized_IS_t *G) {
    uint32_t adm[K];
    uint32_t v = 0;

    // 1: admissible rows, i.e. rows without zeros
    for (uint32_t row = 0; row < K; row++) {
        adm[v] = row;
        v += !row_contains_zero(G->values[row]);
    }
    if (v == 0) { return 0; }

    static FQ_ELEM inv[K][N_K_pad] __attribute__((aligned(32)));
    static FQ_ELEM P[K][K];

    // 2: entrywise inverse of A_Adm
    for (uint32_t i = 0; i < v; i++) {
        row_inv2(inv[i], G->values[adm[i]]);
    }

    // 3: P = A_Adm * inv^T
    for (uint32_t i = 0; i < v; i++) {
        for (uint32_t j = 0; j < v; j++) {
            P[i][j] = row_dot(G->values[adm[i]], inv[j]);
        }
    }

    // 4-6: psi(i) = (sum_j Q_ij^h)_{h=1..d} with Q_ij = P_ij * P_ji, unique arg min
    FQ_ELEM best[EDO_D] = {0};
    uint32_t best_i = 0, unique = 0;
    for (uint32_t i = 0; i < v; i++) {
        FQ_ELEM psi[EDO_D] = {0};
        for (uint32_t j = 0; j < v; j++) {
            const FQ_ELEM q = fq_mul(P[i][j], P[j][i]);
            FQ_ELEM qh = q;
            for (uint32_t h = 0; h < EDO_D; h++) {
                psi[h] = fq_add(psi[h], qh);
                qh = fq_mul(qh, q);
            }
        }

        uint32_t h = 0;
        while ((h < EDO_D) && (psi[h] == best[h])) { h++; }
        if ((i == 0) || ((h < EDO_D) && (psi[h] < best[h]))) {
            memcpy(best, psi, sizeof(psi));
            best_i = i;
            unique = 1;
        } else if (h == EDO_D) {
            unique = 0;
        }
    }
    if (!unique) { return 0; }

    // 7: Phi(A, i*)
    for (uint32_t row = 0; row < K; row++) {
        row_mul3(G->values[row], G->values[row], inv[best_i]);
        if (row_all_same(G->values[row])) { continue; }

        FQ_ELEM s = fq_red(row_acc(G->values[row]));
        if (s != 0) {
            s = fq_inv_non_ct(s);
        } else {
            s = row_acc_inv(G->values[row]);
            if (s == 0) { return 0; }
        }

        row_mul(G->values[row], s);
    }

    if (SortRows(G, K, NULL) == 0) { return 0; }
    SortCols(G, K);
    return 1;
}


int test_edo_cf(void) {
    normalized_IS_t G1, G2, G3;
    uint32_t fails = 0;
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

        const int ret1 = EdoardoCF(&G1);
        const int ret2 = EdoardoCF(&G2);
        if (ret1 != ret2) {
            printf("error ret edo cf\n");
            return 1;
        }

        if (ret1 == 0) {
            fails++;
            continue;
        }

        // the output is in the orbit of the input, so CF''(CF''(G)) == CF''(G)
        normalized_copy(&G3, &G1);
        if (EdoardoCF(&G3) == 0) {
            printf("error ret edo cf idempotent\n");
            return 1;
        }

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G2);
                    printf("error edo cf %d %d\n", i, j);
                    return 1;
                }

                if (G1.values[i][j] != G3.values[i][j]) {
                    printf("error edo cf idempotent %d %d\n", i, j);
                    return 1;
                }
            }
        }
    }

    printf("edo cf: %d/%d failed\n", fails, ITERS);
    return 0;
}

int bench_edo_cf(void) {
    normalized_IS_t G1, G2;
    uint64_t c_edo = 0, c_cf = 0;
    uint32_t f_edo = 0, f_cf = 0;

    setup_cycle_counter();
    for (uint32_t k = 0; k < BENCH_ITERS; k++) {
        normalized_sf(&G1);
        normalized_copy(&G2, &G1);

        uint64_t t = read_cycle_counter();
        f_edo += EdoardoCF(&G1) == 0;
        c_edo += read_cycle_counter() - t;

        t = read_cycle_counter();
        f_cf += CF(&G2) == 0;
        c_cf += read_cycle_counter() - t;
    }

    printf("EdoardoCF:  %llu cycles, %d/%d failed\n", (unsigned long long)(c_edo / BENCH_ITERS), f_edo, BENCH_ITERS);
    printf("ImprovedCF: %llu cycles, %d/%d failed\n", (unsigned long long)(c_cf / BENCH_ITERS), f_cf, BENCH_ITERS);
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
    // if (test_compute_canonical_form_type5_v2()) return 1;
    // if (test_compute_canonical_form_type5_gaus()) return 1;

    if (test_edo_cf()) { return 1; }
    if (bench_edo_cf()) { return 1; }

    printf("Done, all worked\n");
    return 0;
}
