#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"

#define ITERS (1u << 24u)


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

int main(void) {
    if (bench_sorting()) return 1;
    return 0;
}
