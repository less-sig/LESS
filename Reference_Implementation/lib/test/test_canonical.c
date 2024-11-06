#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


int test_sorting_network(void) {
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);
    for (size_t i = 0; i < s; i++) { d1[i] = s-1-i; }

    int8_sort(d1, s);

    
    for (size_t i = 1; i < s; i++) { 
        if (d1[i-1] > d1[i]) { return 1; }
    }

    free(d1);
    return 0;
}

int test_sorting_network_matrix(void) {
    normalized_IS_t G1, G2;
    permutation_t P_c1, P_c2;

    normalized_sf(&G1);
    permutation_mat_id(&P_c1);
    permutation_mat_id(&P_c2);

    memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));

    row_bubble_sort(&G1, &P_c1);
    row_bitonic_sort(&G2, &P_c2);

    normalized_pretty_print(&G1);
    normalized_pretty_print(&G2);

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G1.values[i][j] != G2.values[i][j]) {
                printf("error\n");
                return 1;
            }
        }
    }

    for (uint32_t i = 0; i < N; i++) {
        if (P_c1.permutation[i] != P_c2.permutation[i]) {
            return 2;
        }
    }

    return 0;
}


int test_compute_information_set_monomial(void) {
    monomial_t M;
    permutation_t P_is;

    monomial_mat_rnd(&M);
    // monomial_mat_pretty_print_name("M", &M);
    permutation_mat_id(&P_is);
    compute_information_set_monomial(&P_is, &M);
    // permutation_pretty_print(&P_is);
    // monomial_mat_pretty_print_name("M", &M);
    for (uint32_t i = 0; i < K; ++i) {
        if (M.permutation[i] != i) return 1;
    }
    return 0;
}

int test_compute_canonical_form_type2(void) {
    normalized_IS_t G;
    diagonal_t D_c;
    permutation_t P_c;

    permutation_mat_id(&P_c);
    diagonal_mat_id(&D_c);
    normalized_sf(&G);

    normalized_pretty_print(&G);
    if (compute_canonical_form_type2(&G, &D_c, &P_c) == 0) return 1;
    normalized_pretty_print(&G);

    return 0;
}

int test_compute_canonical_form_type3(void) {
    normalized_IS_t G;
    permutation_t P_c, P_r;

    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);
    normalized_sf(&G);

    normalized_pretty_print(&G);
    if (compute_canonical_form_type3(&G, &P_r, &P_c) == 0) return 1;
    normalized_pretty_print(&G);

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

    normalized_pretty_print(&G);
    if (compute_canonical_form_type4(&G, &P_r, &D_c, &P_c) == 0) return 1;

    normalized_pretty_print(&G);
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

    normalized_pretty_print(&G);
    if (compute_canonical_form_type5(&G, &D_r, &P_r, &D_c, &P_c) == 0) return 1;

    normalized_pretty_print(&G);
    return 0;
}

int main(void) {
    // generic matrices test
    // if (test_compute_canonical_form_type2()) return 1;
    // if (test_compute_canonical_form_type3()) return 1;
    // if (test_compute_canonical_form_type4()) return 1;
    // if (test_compute_canonical_form_type5()) return 1;
    
    // if (test_sorting_network()) return 1;
    if (test_sorting_network_matrix()) return 1;

    // monomial test
    //if (test_compute_information_set_monomial()) return 1;

    return 0;
}
