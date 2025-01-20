#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fq_arith.h"
#include "codes.h"
#include "test_helpers.h"
#include "test_helpers.c"

#define TESTS (1u << 14)

/// NOTE: these functions are outsourced to this file, to make the
/// optimizied implementation as easy as possible.
/// accumulates a row
/// \param d
/// \return sum(d) for _ in range(N-K)
static inline
FQ_ELEM org_row_acc(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < (N-K); col++) {
        s = fq_add(s, d[col]);
	 }

    return s;
}

/// accumulates the inverse of a row
/// \param d
/// \return sum(d) for _ in range(N-K)
static inline
FQ_ELEM org_row_acc_inv(const FQ_ELEM *d) {
    FQ_ELEM s = 0;
    for (uint32_t col = 0; col < (N-K); col++) {
        s = fq_add(s, fq_inv(d[col]));
	 }

    return s;
}

/// scalar multiplication of a row
/// /param row[in/out] *= s for _ in range(N-K)
/// /param s
static inline
void org_row_mul(FQ_ELEM *row, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (N-K); col++) {
        row[col] = fq_mul(s, row[col]);
    }
}

/// scalar multiplication of a row
/// \param out = s*in[i] for i in range(N-K)
/// \param in
/// \param s
static inline
void org_row_mul2(FQ_ELEM *out, const FQ_ELEM *in, const FQ_ELEM s) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_mul(s, in[col]);
    }
}

///
/// \param out = in1[i]*in2[i] for i in range(N-K)
/// \param in1[in]: vector of length N-K
/// \param in2[in]: vector of length N-K
static inline
void org_row_mul3(FQ_ELEM *out, const FQ_ELEM *in1, const FQ_ELEM *in2) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_mul(in1[col], in2[col]);
    }
}

/// invert a row
/// \param out = in[i]**-1 for i in range(N-K)
/// \param in[in]: vector of length N-K
static inline
void org_row_inv2(FQ_ELEM *out, const FQ_ELEM *in) {
    for (uint32_t col = 0; col < (N-K); col++) {
        out[col] = fq_inv(in[col]);
    }
}

/// \param in[in]: vector of length N-K
/// \return 1 if all elements are the same
///         0 else
static inline
uint32_t org_row_all_same(const FQ_ELEM *in) {
    for (uint32_t col = 1; col < N-K; col++) {
        if (in[col] != in[col - 1]) {
            return 0;
        }
    }
    return 1;
}

/// \param in[in]: vector of length N-K
/// \return 0 if no zero was found
///         1 if the row contains at least a single 0
static inline
uint32_t org_row_contains_zero(const FQ_ELEM *in) {
    for (uint32_t col = 0; col < N-K; col++) {
        if (in[col] == 0) {
            return 1;
        }
    }
    return 0;
}

int test_row_acc() {
    FQ_ELEM a[N_K_pad]={0}, b[N_K_pad]={0};
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, N_K_pad);
        memcpy(b, a, N_K_pad);
        const FQ_ELEM c1 = row_acc(a);
        const FQ_ELEM c2 = org_row_acc(b);
        if (c1 != c2) {
            printf("test_row_acc\n");
            return 1;
        }
    }

    return 0;
}

int test_row_acc_inv() {
    FQ_ELEM a[N_K_pad], b[N_K_pad];
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        memcpy(b, a, K);
        const FQ_ELEM c1 = row_acc_inv(a);
        const FQ_ELEM c2 = org_row_acc_inv(b);
        if (c1 != c2) {
            printf("test_row_acc_inv\n");
            return 1;
        }
    }

    return 0;
}

int test_row_mul() {
    FQ_ELEM a[N_K_pad], b[N_K_pad];
    FQ_ELEM s = 1;
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        memcpy(b, a, K);
        row_mul(a,s);
        org_row_mul(b,s);
        for (uint32_t j = 0; j < K; j++) {
            if (a[j] != b[j]) {
                printf("test_row_mul\n");
                return 1;
            }
        }

        s = (s + 1) % Q;
    }

    return 0;
}

int test_row_mul2() {
    FQ_ELEM a[N_K_pad], c1[N_K_pad], c2[N_K_pad];
    FQ_ELEM s = 1;
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        row_mul2(c1, a,s);
        org_row_mul2(c2,a,s);
        for (uint32_t j = 0; j < K; j++) {
            if (c1[j] != c2[j]) {
                printf("test_row_mul2\n");
                return 1;
            }
        }

        s = (s + 1) % Q;
    }

    return 0;
}

int test_row_mul3() {
    FQ_ELEM a[N_K_pad], b[N_K_pad], c1[N_K_pad], c2[N_K_pad];
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        rand_range_q_elements(b, K);
        row_mul3(c1,a,b);
        org_row_mul3(c2,a,b);
        for (uint32_t j = 0; j < K; j++) {
            if (c1[j] != c2[j]) {
                printf("test_row_mul3\n");
                return 1;
            }
        }
    }

    return 0;
}

int test_row_inv2() {
    FQ_ELEM a[N_K_pad], c1[N_K_pad], c2[N_K_pad];
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        row_inv2(c1, a);
        org_row_inv2(c2,a);
        for (uint32_t j = 0; j < K; j++) {
            if (c1[j] != c2[j]) {
                printf("test_row_inv\n");
                return 1;
            }
        }
    }

    return 0;
}

int test_row_all_same() {
    FQ_ELEM a[N_K_pad];
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        const uint32_t t1 = row_all_same(a);
        const uint32_t t2 = org_row_all_same(a);
        for (uint32_t j = 0; j < K; j++) {
            if (t1 != t2) {
                printf("test_row_all_same\n");
                return 1;
            }
        }
    }

    return 0;
}

int test_row_contains_zero() {
    FQ_ELEM a[N_K_pad];
    for (uint32_t i = 0; i < TESTS; i++) {
        rand_range_q_elements(a, K);
        const uint32_t t1 = row_contains_zero(a);
        const uint32_t t2 = org_row_contains_zero(a);
        if (t1 != t2) {
            printf("test_row_contains_zero\n");
            return 1;
        }
    }

    return 0;
}

void matrix_transpose_simple(normalized_IS_t *o,
                      const normalized_IS_t *in) {
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < K; j++) {
            o->values[j][i] = in->values[i][j];
        }
    }
}
int test_transpose() {
    normalized_IS_t A = {0}, B1, B2;
    // normalized_ind(&A);
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < K; j++) {
            A.values[i][j] = (i+j);
        }
    }
    matrix_transpose_opt((uint8_t *)B1.values, (uint8_t *)A.values, K, K_pad);
    matrix_transpose_simple(&B2, &A);
    for (uint32_t i = 0; i < K; i++) {
        for (uint32_t j = 0; j < K; j++) {
            if (B1.values[i][j] != B2.values[i][j]) {
                printf("test_transpose: %d %d\n", i, j);
                return 1;
            }
        }
    }
    return 0;
}

int main(void) {
    // if (test_row_acc()) return 1;
    // if (test_row_acc_inv()) return 1;
    // if (test_row_mul()) return 1;
    // if (test_row_mul2()) return 1;
    // if (test_row_mul3()) return 1;
    // if (test_row_inv2()) return 1;
    // if (test_row_all_same()) return 1;
    // if (test_row_contains_zero()) return 1;

    if (test_transpose()) return 1;

    printf("ok\n");
    return 0;
}
