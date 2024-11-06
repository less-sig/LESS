#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"

#define ITERS (1u << 17u)


int bench_sorting(void) {
    const size_t s = K;
    FQ_ELEM *d1 = (FQ_ELEM *)malloc(sizeof(FQ_ELEM) * s);

    uint64_t c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-1-i; }

        c -= x86_64_rtdsc();
        int8_sort(d1, s);
        c += x86_64_rtdsc();
    }
    printf("int8_sort: %ld cyc\n", c);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        for (size_t j = 0; j < s; j++) { d1[j] = s-1-i; }

        c -= x86_64_rtdsc();
        qsort((void *)d1, s, sizeof(FQ_ELEM), fqcmp);
        c += x86_64_rtdsc();
    }
    printf("qsort:     %ld cyc\n", c);
    free(d1);
    return 0;
}


int bench_row_sorting(void) {
    normalized_IS_t G1;
    permutation_t P_c1;

    uint64_t c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_sf(&G1);
        c -= x86_64_rtdsc();
        row_bitonic_sort(&G1, &P_c1);
        c += x86_64_rtdsc();
    }
    printf("bitonic: %ld cyc\n", c);

    c = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_sf(&G1);

        c -= x86_64_rtdsc();
        row_bubble_sort(&G1, &P_c1);
        c += x86_64_rtdsc();
    }
    printf("bubble:  %ld cyc\n", c);
    return 0;
}

int main(void) {
    // if (bench_sorting()) return 1;
    if (bench_row_sorting()) return 1;
    return 0;
}
