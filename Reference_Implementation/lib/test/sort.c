#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "sort.h"
#include "fq_arith.h"
#include "monomial_mat.h"
#include "codes.h"

#include "test_helpers.c"
#define TESTS 100

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

    uint8_t L[N] = {0};
    for (uint32_t k = 0; k < TESTS; k++) {
        normalized_sf(&G1);
        memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));
        memcpy((void *)G3.values, (void*)G1.values, sizeof(normalized_IS_t));

        row_bitonic_sort(&G2);
        SortRows(&G3, K, L);

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
    }

    return 0;
}

int test_col_sorting_network_matrix(void) {
    normalized_IS_t G1, G2;

    for (uint32_t k = 0; k < TESTS; k++) {
        normalized_sf(&G1);
        memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));

        lex_sort_cols(&G1);
        SortCols(&G2, K_pad);//, 0, N - K - 1);

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G2.values[i][j]) {
                    normalized_pretty_print(&G1);
                    normalized_pretty_print(&G2);
                    printf("error col sorting network3: %d %d\n", i, j);
                    return 1;
                }
            }
        }
    }
    return 0;
}

int main(void) {
    // if (test_counting_sort()) return 1;
    // if (test_sorting_network()) return 1;
    // if (test_sorting_network_matrix()) return 1;
    if (test_col_sorting_network_matrix()) return 1;

    printf("Done, all worked\n");
    return 0;
}
