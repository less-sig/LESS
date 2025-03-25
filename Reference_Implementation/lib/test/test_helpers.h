#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "monomial_mat.h"


void bitonic_sort_i8(FQ_ELEM *x, const long long n);
void counting_sort_u8(FQ_ELEM *arr, const uint32_t size);
int SortRows_internal(FQ_ELEM *ptr[K],
                            uint32_t P[K],
                            const uint32_t n);
// wrapper struct around the set S_n
typedef struct {
  /* considering the product GQ, permutation[...] stores into the cell with
   * index 0, the position of the DESTINATION of column 0 in G after the
   * computation of GQ.
   */
  POSITION_T permutation[N];
} permutation_t;

// wrapper struct around D_n
typedef struct {
  /* coefficients listed in order of appearance column-wise */
  FQ_ELEM coefficients[N];
} diagonal_t;


void monomial_mat_rnd(monomial_t *res);
void monomial_mat_rnd_unique(monomial_t *res);


int compute_canonical_form_type4(normalized_IS_t *G, const uint8_t *L);
int compute_canonical_form_type3(normalized_IS_t *G, const uint8_t *L);


/* samples a random monomial matrix */
void generator_rnd(generator_mat_t *res);
void generator_sf(generator_mat_t *res);

void normalized_ind(normalized_IS_t *V);
void normalized_sf(normalized_IS_t *V);
void normalized_rng(normalized_IS_t *V);

void generator_pretty_print(const generator_mat_t *const G);
void generator_pretty_print_name(char *name, const generator_mat_t *const G);
void generator_rref_pretty_print_name(char *name,
                                      const rref_generator_mat_t *const G);

void normalized_pretty_print(const normalized_IS_t *const G);


int SortCols_internal_compare(uint8_t *ptr[K],
                                                                 const POSITION_T row_idx,
                                                                 const uint8_t pivot[K]);

int SortRows_internal_compare(uint8_t *ptr[K],
                                               const uint32_t row_idx,
                                               const uint8_t pivot[K]);

int compare_matrices(const normalized_IS_t *V1,
                     const normalized_IS_t *V2,
                     const uint32_t z);

int compute_canonical_form_type3_sub(normalized_IS_t *G,
                                     const uint32_t z);

int compute_canonical_form_type4_sub_v2(normalized_IS_t *G,
                                        const uint32_t z);

int compute_canonical_form_type5_single_row(normalized_IS_t *G,
                                           const uint32_t row);

void lex_minimize(normalized_IS_t *V,
                  POSITION_T dst_col_idx,
                  const generator_mat_t *const G,
                  const POSITION_T col_idx);

//
int lex_compare_column(const generator_mat_t *G1,
                       const generator_mat_t *G2,
                       const POSITION_T col1,
                       const POSITION_T col2);

int lex_compare_col(const normalized_IS_t *G1,
                    const POSITION_T col1,
                    const POSITION_T col2);
//
int lex_compare_with_pivot(normalized_IS_t *V,
                           const POSITION_T col_idx,
                           const FQ_ELEM pivot[K]);

// in place quick sort
void col_lex_quicksort(normalized_IS_t *V,
                       int start,
                       int end);

/* performs lexicographic sorting of the IS complement */
void lex_sort_cols(normalized_IS_t *V);

void column_swap(normalized_IS_t *V,
                 const POSITION_T col1,
                 const POSITION_T col2);


////////////////////////////////////////////////////////////////////////
///                        Permutation                               ///
////////////////////////////////////////////////////////////////////////

void permutation_swap(permutation_t *P, uint32_t i, uint32_t j);
void permutation_cswap(permutation_t *P, uint32_t i, uint32_t j, uintptr_t mask);
void permutation_mat_id(permutation_t *P);
void permutation_mat_rng(permutation_t *P);
void permutation_mat_id_v2(permutation_t *P, const uint32_t max);
void permutation_mat_rng_v2(permutation_t *P, const uint32_t max);
void permutation_pretty_print(const permutation_t *P);


////////////////////////////////////////////////////////////////////////
///                             Diagonal                             ///
////////////////////////////////////////////////////////////////////////
void diagonal_mat_zero(diagonal_t *D);
void diagonal_mat_id(diagonal_t *D);
void diagonal_mat_rnd(diagonal_t *D);
void diagonal_mat_id_v2(diagonal_t *D, uint32_t max);
void diagonal_mat_rnd_v2(diagonal_t *D, uint32_t max);
void diagonal_pretty_print(const diagonal_t *const D);

// defined in monomial.c
void permutation_apply_row(const permutation_t *P, normalized_IS_t *G);
void permutation_apply_col(normalized_IS_t *G, const permutation_t *P);

void diagonal_apply_row(const diagonal_t *P, normalized_IS_t *G);
void diagonal_apply_col(normalized_IS_t *G, const diagonal_t *P);



#endif //TEST_HELPERS_H
