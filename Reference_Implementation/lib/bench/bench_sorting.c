#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"
#include "sort.h"
#include "../test/test_helpers.c"

#define ITERS (1u << 15u)


int bench_sorting(void) {
    printf("int : \n");
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *) malloc(sizeof(FQ_ELEM) * s);

    unsigned c = 0, c1;
    for (unsigned i = 0; i < ITERS; i++) {
        for (unsigned j = 0; j < s; j++) { d1[j] = s-1-i; }

        c -= read_cycle_counter();
        counting_sort_u8(d1, s);
        c += read_cycle_counter();
    }
    c1 = c / ITERS;
    printf("int8_sort: %u cyc\n", c1);

    free(d1);
    return 0;
}


int bench_row_sorting(void) {
    printf("row: \n");
    normalized_IS_t G1;
    uint8_t L[N] = {0};

    unsigned c = 0, c1;
    uint32_t ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
        normalized_rng(&G1);
        c -= read_cycle_counter();
        ctr += SortRows(&G1, K, L);
        c += read_cycle_counter();
    }
    c1 = c / ITERS;
    printf("quick: %u cyc, %d\n", c1, ctr);

    c=0; ctr=0;
    for (unsigned i = 0; i < ITERS; i++) {
        normalized_rng(&G1);
        c -= read_cycle_counter();
        ctr += row_quick_sort_recursive(&G1, K);
        c += read_cycle_counter();
    }
    c = c / ITERS;
    printf("quickr : %u cyc, %d\n", c, ctr);
    printf("factor %lf\n", ((double) c) / (double) c1);
    return 0;
}

int bench_col_sorting(void) {
    printf("col: \n");
    normalized_IS_t G1;

    unsigned c = 0, c1;
    for (unsigned i = 0; i < ITERS; i++) {
        normalized_rng(&G1);

        c -= read_cycle_counter();
        SortCols(&G1, K);
        c += read_cycle_counter();
    }
    c1 = c / ITERS;
    printf("quickT: %u cyc\n", c1);


    // c = 0;
    // for (unsigned i = 0; i < ITERS; i++) {
    //     normalized_rng(&G1);

    //     c -= read_cycle_counter();
    //     lex_sort_cols(&G1);
    //     c += read_cycle_counter();
    // }
    // c = c / ITERS;
    // printf("normal:  %u cyc\n", c);
    // printf("factor %lf\n", ((double) c) / (double) c1);


    return 0;
}

int main(void) {
    // if (bench_sorting()) return 1;
    // if (bench_row_sorting()) return 1;
    if (bench_col_sorting()) return 1;
    return 0;
}
