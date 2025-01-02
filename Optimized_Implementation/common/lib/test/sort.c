#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "sort.h"
#include "fq_arith.h"
#include "monomial_mat.h"
#include "codes.h"

#define TESTS 100

// int test_sorting_network(void) {
//     const size_t s = K;
//     FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);
//     for (uint32_t k = 0; k < TESTS; k++) {
//         fq_star_rnd_elements(d1, K);
//
//         sortingnetwork(d1, s);
//         for (size_t i = 1; i < s; i++) {
//             if (d1[i-1] > d1[i]) {
//                 printf("error sorting_network: %d\n", i);
//                 return 1;
//             }
//         }
//     }
//
//     free(d1);
//     return 0;
// }

int test_sorting_network_matrix(void) {
    normalized_IS_t G1, G2, G3;

    for (uint32_t k = 0; k < TESTS; k++) {
        normalized_sf(&G1);
        memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));
        memcpy((void *)G3.values, (void*)G1.values, sizeof(normalized_IS_t));

        row_bitonic_sort(&G2);
        row_quick_sort(&G3, K);

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G2.values[i][j] != G3.values[i][j]) {
                    printf("error test_sorting_network_matrix\n");
                    return 1;
                }
            }
        }
    }

    return 0;
}


int test_col_sorting_network_matrix(void) {
    normalized_IS_t G1, G2, G3, G4;
    int ret = 0;

    for (uint32_t k = 0; k < TESTS; k++) {
        normalized_sf(&G1);
        memcpy((void *)G2.values, (void*)G1.values, sizeof(normalized_IS_t));
        memcpy((void *)G3.values, (void*)G1.values, sizeof(normalized_IS_t));
        memcpy((void *)G4.values, (void*)G1.values, sizeof(normalized_IS_t));

        lex_sort_cols(&G1);
        col_quicksort_transpose(&G3, K);
        col_bitonic_sort_transpose(&G4);

        for (uint32_t i = 0; i < K; i++) {
            for (uint32_t j = 0; j < N-K; j++) {
                if (G1.values[i][j] != G3.values[i][j]) {
                    printf("error col sorting quick\n");
                    ret = 1;
                }

                if (G1.values[i][j] != G4.values[i][j]) {
                    printf("error col sorting network bitonic\n");
                    ret = 1;
                }
            }
        }
    }

    return ret;
}

int main(void) {
    // if (test_sorting_network()) return 1;
    if (test_sorting_network_matrix()) return 1;
    if (test_col_sorting_network_matrix()) return 1;

    printf("ok\n");
    return 0;
}
