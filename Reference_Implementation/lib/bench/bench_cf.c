#include <stdio.h>
#include <math.h>

#include "canonical.h"
#include "cycles.h"
#include "sort.h"

#define ITERS (1u << 14u)


/// NOTE: only for testing
void normalized_rng(normalized_IS_t *V) {
    randombytes((uint8_t *)V->values, K*(N-K));
    const uint8_t mask = 0x7F;
    for (unsigned i = 0; i < K; ++i) {
        for (unsigned j = 0; j < K; ++j) {
            V->values[i][j] &= mask;
        }
    }
}

int bench_cf3(void) {
    normalized_IS_t G1, G2;

	printf("cf3:\n");
    uint64_t c = 0;
    unsigned ctr = 0;
	normalized_sf(&G2);
    for (unsigned i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

    	c -= read_cycle_counter();
        ctr += compute_canonical_form_type3(&G1);
        c += read_cycle_counter();
    }
    unsigned c1 = (unsigned)(c / ITERS);
    printf("non-ct: %u cyc, ctr: %u\n", c1, ctr);

    c = 0; ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
        normalized_copy(&G1, &G2);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type3_ct(&G1);
        c += read_cycle_counter();
    }

    unsigned c2 = c / ITERS;
    printf("ct: %u cyc, ctr: %u\n", c2, ctr);
    printf("factor %lf\n\n", (double)c2/(double)c1);
    return 0;
}

int bench_cf4(void) {
    normalized_IS_t G2, G1;

	printf("cf4:\n");
	normalized_sf(&G2);
    uint64_t c = 0;
    unsigned ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type4(&G1);
        c += read_cycle_counter();
    }
    unsigned c1 = (unsigned) (c / ITERS);
    printf("non ct: %u cyc, ctr: %u\n", c1, ctr);

    c = 0; ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
		normalized_copy(&G1, &G2);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type4_ct(&G1);
        c += read_cycle_counter();
    }

    unsigned c2 = (unsigned) (c / ITERS);
    printf("ct: %u cyc, ctr: %u\n", c2, ctr);
    printf("factor %lf\n\n", (double)c2/(double)c1);
    return 0;
}

int bench_cf5(void) {
    normalized_IS_t G1, G2;

	normalized_sf(&G2);
	printf("cf5:\n");
    uint64_t c = 0;
    unsigned ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
		// normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type5(&G1);
        c += read_cycle_counter();
    }
    unsigned c1 = (unsigned) (c / ITERS);
    printf("non ct: %u cyc, ctr: %u\n\n", c1, ctr);


    c = 0; ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
        // normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type5_popcnt(&G1);
        c += read_cycle_counter();
    }

    unsigned c2 = (unsigned) (c / ITERS);
    printf("pop: %u cyc, ctr: %u\n", c2, ctr);
    printf("factor %lf\n\n", (double) c2 / (double) c1);


    c = 0; ctr = 0;
    for (unsigned i = 0; i < ITERS; i++) {
        // normalized_copy(&G1, &G2);
        normalized_rng(&G1);

        c -= read_cycle_counter();
        ctr += compute_canonical_form_type5_fastest(&G1);
        c += read_cycle_counter();
    }

    unsigned c3 = (unsigned) (c / ITERS);
    printf("fast: %u cyc, ctr: %u\n", c3, ctr);
    printf("factor %lf\n\n", (double) c3 / (double) c1);
    return 0;
}

int main(void) {
    // if (bench_cf3()) return 1;
    // if (bench_cf4()) return 1;
    if (bench_cf5()) return 1;
    return 0;
}
