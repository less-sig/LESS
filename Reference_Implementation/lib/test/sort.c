#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "sort.h"
#include "fq_arith.h"
#include "monomial_mat.h"
#include "codes.h"


int test_counting_sort(void) {
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);
    for (size_t i = 0; i < s; i++) { d1[i] = s-1-i; }

    counting_sort_u8(d1, s);
    for (size_t i = 1; i < s; i++) {
        if (d1[i-1] > d1[i]) {
            printf("error counting sort\n");
            return 1;
        }
    }

    free(d1);
    return 0;
}

int test_sorting_network(void) {
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);
    for (size_t i = 0; i < s; i++) { d1[i] = s-1-i; }

    bitonic_sort_i8(d1, s);
    for (size_t i = 1; i < s; i++) {
        if (d1[i-1] > d1[i]) {
            printf("error sorting_network\n");
            return 1;
        }
    }

    free(d1);
    return 0;
}

int test_sorting_network_matrix(void) {
    normalized_IS_t G1, G2, G3;

    normalized_sf(&G1);
    memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));
    memcpy((void *)G3.values, (void*)G1.values, sizeof(normalized_IS_t));

    row_bitonic_sort(&G2);
    row_quick_sort(&G3, K);

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G2.values[i][j] != G3.values[i][j]) {
                normalized_pretty_print(&G2);
                normalized_pretty_print(&G3);
                printf("error test_sorting_network_matrix\n");
                return 1;
            }
        }
    }

    return 0;
}


int test_col_sorting_network_matrix(void) {
    normalized_IS_t G1, G2, G3, G4;

    normalized_sf(&G1);
    memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));
    memcpy((void *)G3.values, (void*)G1.values, sizeof(normalized_IS_t));
    memcpy((void *)G4.values, (void*)G1.values, sizeof(normalized_IS_t));

    //normalized_pretty_print(&G1);
    lex_sort_cols(&G1);
    col_quicksort_transpose(&G3, K);//, 0, N - K - 1);
    col_bitonic_sort_transpose(&G4);

    // normalized_pretty_print(&G1);
    // normalized_pretty_print(&G2);
    // normalized_pretty_print(&G3);
    // normalized_pretty_print(&G4);

    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < N-K; j++) {
            if (G1.values[i][j] != G3.values[i][j]) {
                printf("error col sorting network3\n");
                return 1;
            }

            if (G1.values[i][j] != G4.values[i][j]) {
                printf("error col sorting network4\n");
                return 1;
            }
        }
    }
    return 0;
}

int main(void) {
    if (test_counting_sort()) return 1;
    if (test_sorting_network()) return 1;
    if (test_sorting_network_matrix()) return 1;
    if (test_col_sorting_network_matrix()) return 1;

    printf("Done, all worked\n");
    return 0;
}
