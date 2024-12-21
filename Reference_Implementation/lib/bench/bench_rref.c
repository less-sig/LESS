#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"
#include "sort.h"
#include "test_helpers.h"

#define ITERS (1u << 10u)


void generator_copy(generator_mat_t *V1,
                    const generator_mat_t *V2) {
    memcpy(V1->values, V2->values, sizeof(generator_mat_t));
}

void generator_rnd_fullrank(generator_mat_t *G,
                            uint8_t *is_pivot_column) {
     do {
         generator_rnd(G);
         memset(is_pivot_column,0, N);
     } while ( generator_RREF(G,is_pivot_column) == 0);
}

int bench_rref(void) {
    generator_mat_t G2, G1;
    uint8_t is_pivot_column[N] = {0};
    uint8_t g_initial_pivot_flags [N] = {0};
    uint8_t g_permuted_pivot_flags [N];

    monomial_t q;
    generator_rnd_fullrank(&G2, is_pivot_column);

	printf("rref:\n");
    uint64_t c = 0, c1, ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        monomial_mat_rnd(&q);
        generator_monomial_mul(&G1, &G2, &q);

    	c -= read_cycle_counter();
        ctr += generator_RREF(&G1, is_pivot_column);
        c += read_cycle_counter();
    }
    c1 = c/ITERS;
    printf("normal: %ld cyc, ctr: %ld\n", c1, ctr);

    generator_rnd_fullrank(&G2, g_initial_pivot_flags);
    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        memset(is_pivot_column, 0, N);
        monomial_mat_rnd(&q);
        generator_monomial_mul(&G1, &G2, &q);
        for (uint32_t j = 0; j < N; j++) {
            g_permuted_pivot_flags[q.permutation[j]] = g_initial_pivot_flags[j];
        }

        c -= read_cycle_counter();
        ctr += generator_RREF_pivot_reuse(&G1, is_pivot_column, g_permuted_pivot_flags, VERIFY_PIVOT_REUSE_LIMIT);
        c += read_cycle_counter();
    }

    c = c/ITERS;
    printf("reuse: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}


int main(void) {
    if (bench_rref()) return 1;
    return 0;
}
