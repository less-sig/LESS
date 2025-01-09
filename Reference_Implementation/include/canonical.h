#ifndef LESS_CANONICAL_H
#define LESS_CANONICAL_H

#include "monomial_mat.h"
#include "codes.h"

////////////////////////////////////////////////////////////////////////
///                        Canonical Forms                           ///
////////////////////////////////////////////////////////////////////////
int compute_canonical_form_type5(normalized_IS_t *G);
int compute_canonical_form_type5_popcnt(normalized_IS_t *G);
int CF(normalized_IS_t *G);
void blind(normalized_IS_t *G, SHAKE_STATE_STRUCT *prng);
#endif //LESS_CANONICAL_H
