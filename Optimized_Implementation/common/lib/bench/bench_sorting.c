#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"

#define ITERS (1u << 13u)


int bench_sorting(void) {
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);

    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-1-i; }

        c -= x86_64_rtdsc();
        int8_sort(d1, s);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("int8_sort: %ld cyc\n", c1);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-1-j; }

        c -= x86_64_rtdsc();
        qsort((void *)d1, s, sizeof(FQ_ELEM), fqcmp);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("qsort:     %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-j; }

        c -= x86_64_rtdsc();
        counting_sort_u8(d1, s);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("count:     %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);
    free(d1);
    return 0;
}


int bench_row_sorting(void) {
    normalized_IS_t G1;
    permutation_t P_c1;

    normalized_sf(&G1);
    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_sf(&G1);
        c -= x86_64_rtdsc();
        row_bitonic_sort(&G1, &P_c1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("bitonic: %ld cyc\n", c1);

    normalized_sf(&G1);
    c=0;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_sf(&G1);

        c -= x86_64_rtdsc();
        row_quick_sort(&G1, &P_c1);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("quick: %ld cyc\n", c);

    //c = 0;
    //for (uint64_t i = 0; i < ITERS; i++) {
    //    // normalized_sf(&G1);

    //    c -= x86_64_rtdsc();
    //    row_bubble_sort(&G1, &P_c1);
    //    c += x86_64_rtdsc();
    //}
    //c = c/ITERS;
    //printf("bubble:  %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);
    return 0;
}

int bench_col_sorting(void) {
    normalized_IS_t G1;
    permutation_t P_c1;

    uint64_t c = 0, c1;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_sf(&G1);

        c -= x86_64_rtdsc();
        col_bitonic_sort(&G1, &P_c1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("bitonic: %ld cyc\n", c1);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_sf(&G1);

        c -= x86_64_rtdsc();
        canonical_col_lex_quicksort(&G1, 0, N-K-1, &P_c1);
        c += x86_64_rtdsc();
    }
    c = c/ITERS;
    printf("quick:  %ld cyc\n", c);
    printf("factor %lf\n", (double)c/(double)c1);
    return 0;
}

int main(void) {
    //if (bench_sorting()) return 1;
    if (bench_row_sorting()) return 1;
    // if (bench_col_sorting()) return 1;
    return 0;
}
