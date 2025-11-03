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


/// Does not really sort a single row into `out`. It computes the
/// histogram of the input array `in`;
/// \param out[out]: buffer of 127 elements count the number occurrences
///     of each Fq element in `in`
/// \param in[in]: row of a matrix containing Fq elements
/// \param len[in]: length of the row
void sort(uint8_t *out,
          const uint8_t *in,
          uint32_t len);

/// \param row1[in]: pointer to the first row
/// \param row2[in]: pointer to the second row
/// \return: 0 if multiset(row1) == multiset(row2)
///          x if multiset(row1) > multiset(row2)
///         -x if multiset(row1) < multiset(row2)
int compare_rows(const FQ_ELEM *row1,
                 const FQ_ELEM *row2);

/// sort the rows of G
/// \param G[in]: pointer to a generator matrix (only non IS part) to sort.
/// \param n[in]: number of rows to sort
/// \param L[in]: pointer to the currently shortest row
/// \return 1 on success
///			0 if two rows generate the same multiset
int SortRows(normalized_IS_t *G,
             uint32_t n,
             const uint8_t *L);

/// sort the cols of G
/// \param G[in]: pointer to a generator matrix (only non IS part) to sort.
/// \param z[in]: number of columns to sort
void SortCols(normalized_IS_t *G,
              uint32_t z);

#endif //SORT_H
