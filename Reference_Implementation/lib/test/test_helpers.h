#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "monomial_mat.h"

void monomial_mat_rnd(monomial_t *res);
void monomial_mat_rnd_unique(monomial_t *res);

///
static inline
FQ_ELEM fq_add(const FQ_ELEM x, const FQ_ELEM y) {
      return (x + y) % Q;
}


#endif //TEST_HELPERS_H
