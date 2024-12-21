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
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);

    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-1-i; }

        c -= x86_64_rtdsc();
        counting_sort_u8(d1, s);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("int8_sort: %ld cyc\n", c1);

#ifdef USE_AVX2 // TODO remove
    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-j; }

        c -= x86_64_rtdsc();
        sortingnetwork(d1, s);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("network:     %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);
#endif
    free(d1);
    return 0;
}


int bench_row_sorting(void) {
    printf("row: \n");
    normalized_IS_t G1;

    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_rng(&G1);
        c -= x86_64_rtdsc();
        row_bitonic_sort(&G1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("bitonic: %ld cyc\n", c1);

    c=0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_rng(&G1);
        c -= x86_64_rtdsc();
        row_quick_sort(&G1, K);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("quick: %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);
    return 0;
}

int bench_col_sorting(void) {
    printf("col: \n");
    normalized_IS_t G1;

    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_rng(&G1);

        c -= x86_64_rtdsc();
        col_bitonic_sort_transpose(&G1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("bitonicT: %ld cyc\n", c1);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_rng(&G1);

        c -= x86_64_rtdsc();
        col_quicksort_transpose(&G1, K);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("quickT:  %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);


    return 0;
}

int main(void) {
    if (bench_sorting()) return 1;
    if (bench_row_sorting()) return 1;
    if (bench_col_sorting()) return 1;
    return 0;
}
