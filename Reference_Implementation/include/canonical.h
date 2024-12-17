#ifndef LESS_CANONICAL_H
#define LESS_CANONICAL_H

#include "monomial_mat.h"
#include "codes.h"

////////////////////////////////////////////////////////////////////////
///                        Canonical Forms                           ///
////////////////////////////////////////////////////////////////////////

int compute_canonical_form_type3(normalized_IS_t *G);
int compute_canonical_form_type3_ct(normalized_IS_t *G);

int compute_canonical_form_type4(normalized_IS_t *G);
int compute_canonical_form_type4_ct(normalized_IS_t *G);

int compute_canonical_form_type5(normalized_IS_t *G);
int compute_canonical_form_type5_ct(normalized_IS_t *G);
int compute_canonical_form_type5_popcnt(normalized_IS_t *G);

int cf5(normalized_IS_t *G);
int cf5_nonct(normalized_IS_t *G);

void blind(normalized_IS_t *G,
           SHAKE_STATE_STRUCT *prng);
#endif //LESS_CANONICAL_H
