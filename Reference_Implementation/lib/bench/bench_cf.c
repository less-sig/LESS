#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"
#include "sort.h"

#define ITERS (1u << 4u)


/// NOTE: only for testing
void normalized_rng(normalized_IS_t *V) {
    randombytes((uint8_t *)V->values, K*(N-K));
    const uint8_t mask = 0x7F;
    for (uint32_t i = 0; i < K; ++i) {
        for (uint32_t j = 0; j < K; ++j) {
            V->values[i][j] &= mask;
        }
    }
}

int bench_cf3(void) {
    normalized_IS_t G1, G2;

	printf("cf3:\n");
    uint64_t c = 0, c1, ctr = 0;
	normalized_sf(&G2);
    for (uint64_t i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

    	c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type3(&G1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("non-ct: %ld cyc, ctr: %ld\n", c1, ctr);

    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        normalized_copy(&G1, &G2);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type3_ct(&G1);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("ct: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}

int bench_cf4(void) {
    normalized_IS_t G2, G1;

	printf("cf4:\n");
	normalized_sf(&G2);
    uint64_t c = 0, c1, ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type4(&G1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("non ct: %ld cyc, ctr: %ld\n", c1, ctr);

    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type4_ct(&G1);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("ct: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}

int bench_cf5(void) {
    normalized_IS_t G1, G2;

	normalized_sf(&G2);
	printf("cf5:\n");
    uint64_t c = 0, c1, ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
		// normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type5(&G1);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("non ct: %ld cyc, ctr: %ld\n\n", c1, ctr);


    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type5_popcnt(&G1);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("pop: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);


    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type5_fastest(&G1);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("fast: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}

int main(void) {
    // if (bench_cf3()) return 1;
    // if (bench_cf4()) return 1;
    if (bench_cf5()) return 1;
    return 0;
}
