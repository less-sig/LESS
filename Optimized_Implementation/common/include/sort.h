#ifndef SORT_H
#define SORT_H

#include "monomial_mat.h"
#include "codes.h"

////////////////////////////////////////////////////////////////////////
///                             Sorting                              ///
////////////////////////////////////////////////////////////////////////
int fqcmp(const void *a,const void *b);

void bitonic_sort_i8(FQ_ELEM *x, const long long n);
void counting_sort_u8(uint8_t *arr, const size_t size);

#ifdef USE_AVX2
void sortingnetwork(uint8_t *arr, const size_t size);
#endif

int compare_matrices(const normalized_IS_t *V1,
                     const normalized_IS_t *V2);

int row_bubble_sort(normalized_IS_t *G);
int row_bitonic_sort(normalized_IS_t *G);
int row_quick_sort(normalized_IS_t *G);

void col_bitonic_sort_transpose(normalized_IS_t *V);
void col_quicksort_transpose(normalized_IS_t *V);

#endif //SORT_H
