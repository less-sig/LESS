#ifndef LESS_CANONICAL_H
#define LESS_CANONICAL_H

#include "monomial_mat.h"
#include "codes.h"

////////////////////////////////////////////////////////////////////////
///                        Canonical Forms                           ///
////////////////////////////////////////////////////////////////////////
int compute_canonical_form_type2(normalized_IS_t *G,
                diagonal_t *D_c, permutation_t *P_c);

int compute_canonical_form_type3(normalized_IS_t *G, permutation_t *P_r, permutation_t *P_c);

int compute_canonical_form_type4(normalized_IS_t *G,
                                 permutation_t *P_r, diagonal_t *D_c,
                                 permutation_t *P_c);

int compute_canonical_form_type5(normalized_IS_t *G,
                                 diagonal_t *D_r, permutation_t *P_r,
                                 diagonal_t *d_c, permutation_t *P_c);

int compute_information_set_monomial(permutation_t *P_is, monomial_t *G);
#endif //LESS_CANONICAL_H