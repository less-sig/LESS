#ifndef LESS_CANONICAL_H
#define LESS_CANONICAL_H

#include "monomial_mat.h"
#include "codes.h"




void canonical_col_lex_quicksort(normalized_IS_t *V,
                                 const int start,
                                 const int end);

////////////////////////////////////////////////////////////////////////
///                        Canonical Forms                           ///
////////////////////////////////////////////////////////////////////////

int compute_canonical_form_type3(normalized_IS_t *G);

int compute_canonical_form_type4(normalized_IS_t *G);

int compute_canonical_form_type4_non_ct(normalized_IS_t *G);

int compute_canonical_form_type5(normalized_IS_t *G);

int cf5(normalized_IS_t *G);
#endif //LESS_CANONICAL_H
