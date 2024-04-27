#include "monomial_mat.h"
#include "codes.h"
#include "canonical.h"


int test_compute_information_set_monomial() {
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

int test_compute_canonical_form_type2() {
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

int test_compute_canonical_form_type3() {
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

int test_compute_canonical_form_type4() {
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

int test_compute_canonical_form_type5() {
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

int main() {
    // generic matrices test
    // if (test_compute_canonical_form_type2()) return 1;
    // if (test_compute_canonical_form_type3()) return 1;
    // if (test_compute_canonical_form_type4()) return 1;
    if (test_compute_canonical_form_type5()) return 1;

    // monomial test
    //if (test_compute_information_set_monomial()) return 1;

    return 0;
}