#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "canonical.h"
#include "cycles.h"
#include "sort.h"

#define ITERS (1u << 5u)



int cf3_bubble(normalized_IS_t *G,
               permutation_t *P_r,
               permutation_t *P_c) {
    if (row_bubble_sort(G, P_r) == 0) { return 0; }
    canonical_col_lex_quicksort(G, 0, N-K-1, P_c);
    return 1;
}



int cf4_bubble(normalized_IS_t *G,
               permutation_t *P_r, diagonal_t *D_c,
               permutation_t *P_c) {
	for (uint32_t row = 0; row < K; row++) {
        // if we cant find a power
		FQ_ELEM s = 0, sp = 0, tmp, q2=Q-2;
		for (uint32_t col = 0; col < (N-K); col++) {
			s = fq_add(s, G->values[row][col]);
			tmp = fq_pow(G->values[row][col], q2);
			sp = fq_add(sp, tmp);
		}

		if (s != 0) {
			s = fq_inv(s);
		} else {
			s = sp;
			if (s == 0) {
				return -1;
			}
		}

		for (uint32_t col = 0; col < (N-K); col++) {
			G->values[row][col] = fq_mul(s, G->values[row][col]);
		}

	}

	return cf3_bubble(G, P_r, P_c);
}

int cf5_bubble(normalized_IS_t *G,
              diagonal_t *D_r, permutation_t *P_r,
              diagonal_t *D_c, permutation_t *P_c) {
	normalized_IS_t Aj, smallest;
    int touched = 0;

	// init the output matrix to some `invalid` data
	memset(&smallest, -1, K*(N-K));


	for (uint32_t col = 0; col < N-K; col++) {
        if (normalized_is_zero_in_column(G, col)) { continue; }
		memcpy((void *)&Aj, G, K*(N-K));
        touched = 1;

		// first scale all rows
		for (uint32_t row = 0; row < K; row++) {
            FQ_ELEM tmp = fq_inv(Aj.values[row][col]);
			normalized_mat_scale_row(&Aj, row, tmp);
            // D_r->coefficients[row] = fq_mul(D_r->coefficients[row], tmp);
		}

		cf4_bubble(&Aj, P_r, D_c, P_c);

		if (compare_matrices(&Aj, &smallest) < 0) {
			memcpy(&smallest, &Aj, K*(N-K));
		}
	}

    if (!touched) { return 0; }

	memcpy(G, &smallest, K*(N-K));
	return 1;
}



int bench_cf3(void) {
    normalized_IS_t G1;
    permutation_t P_c, P_r;
    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);

	printf("cf3:\n");
    uint64_t c = 0, c1, ctr = 0;
	normalized_sf(&G1);
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_sf(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type3(&G1, &P_r, &P_c);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("ct: %ld cyc, ctr: %ld\n", c1, ctr);

	normalized_sf(&G1);
    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_sf(&G1);

        c -= x86_64_rtdsc();
        ctr += cf3_bubble(&G1, &P_r, &P_c);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("bubble: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}

int bench_cf4(void) {
    normalized_IS_t G1;
    permutation_t P_c, P_r;
    diagonal_t D_c;
    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);
    diagonal_mat_id(&D_c);

	printf("cf4:\n");
	normalized_sf(&G1);
    uint64_t c = 0, c1, ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        //normalized_sf(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type4(&G1, &P_r, &D_c, &P_c);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("ct: %ld cyc, ctr: %ld\n", c1, ctr);

	normalized_sf(&G1);
    c = 0; ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        //normalized_sf(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type4_non_ct(&G1, &P_r, &D_c, &P_c);
        c += x86_64_rtdsc();
    }

    c = c/ITERS;
    printf("non_ct: %ld cyc, ctr: %ld\n", c, ctr);
    printf("factor %lf\n\n", (double)c/(double)c1);

	// normalized_sf(&G1);
    // c = 0; ctr = 0;
    // for (uint64_t i = 0; i < ITERS; i++) {
    //     //normalized_sf(&G1);

    //     c -= x86_64_rtdsc();
    //     ctr += cf4_bubble(&G1, &P_r, &D_c, &P_c);
    //     c += x86_64_rtdsc();
    // }

    // c = c/ITERS;
    // printf("bubble: %ld cyc, ctr: %ld\n", c, ctr);
    // printf("factor %lf\n\n", (double)c/(double)c1);
    return 0;
}

int bench_cf5(void) {
    normalized_IS_t G1;
    permutation_t P_c, P_r;
    diagonal_t D_c, D_r;
    permutation_mat_id(&P_c);
    permutation_mat_id(&P_r);
    diagonal_mat_id(&D_c);
    diagonal_mat_id(&D_r);

	normalized_sf(&G1);
	printf("cf5:\n");
    uint64_t c = 0, c1, ctr = 0;
    for (uint64_t i = 0; i < ITERS; i++) {
        // normalized_sf(&G1);

        c -= x86_64_rtdsc();
        ctr += compute_canonical_form_type5(&G1, NULL, &P_r, NULL, &P_c);
        c += x86_64_rtdsc();
    }
    c1 = c/ITERS;
    printf("ct: %ld cyc, ctr: %ld\n\n", c1, ctr);

	// normalized_sf(&G1);
    // c = 0; ctr = 0;
    // for (uint64_t i = 0; i < ITERS; i++) {
    //     // normalized_sf(&G1);

    //     c -= x86_64_rtdsc();
    //     ctr += cf5_bubble(&G1, NULL, &P_r, NULL, &P_c);
    //     c += x86_64_rtdsc();
    // }

    // c = c/ITERS;
    // printf("bubble: %ld cyc, ctr: %ld\n", c, ctr);
    // printf("factor %lf\n", (double)c/(double)c1);
    return 0;
}

int main(void) {
    // if (bench_cf3()) return 1;
    // if (bench_cf4()) return 1;
    if (bench_cf5()) return 1;
    return 0;
}
