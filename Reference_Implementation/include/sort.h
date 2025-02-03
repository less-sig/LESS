/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Floyd Zweydinger <zweydfg8+github@rub.de>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/#ifndef SORT_H
#define SORT_H

#include "codes.h"

////////////////////////////////////////////////////////////////////////
///                             Sorting                              ///
////////////////////////////////////////////////////////////////////////
int fqcmp(const void *a,const void *b);


void bitonic_sort_i8(FQ_ELEM *x, const long long n);
void counting_sort_u8(FQ_ELEM *arr, const uint32_t size);

#ifdef USE_AVX2
void sortingnetwork(uint8_t *arr, const size_t size);
#endif

void sort(uint8_t *out, const uint8_t *in, const uint32_t len);

int compare_rows(const FQ_ELEM *row1, const FQ_ELEM *row2);

int SortRows(normalized_IS_t *G, const uint32_t n, const uint8_t *L);
int SortRows_internal(FQ_ELEM *ptr[K],
                            uint32_t P[K],
                            const uint32_t n);

void col_bitonic_sort_transpose(normalized_IS_t *V);
void SortCols(normalized_IS_t *V,
                             const uint32_t z);
void col_quicksort(normalized_IS_t *V,
                   const uint32_t z);
#endif //SORT_H
