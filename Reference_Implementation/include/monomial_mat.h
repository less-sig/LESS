/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
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
 **/

#pragma once

#include "parameters.h"
#include "fq_arith.h"

/* Structure representing a monomial matrix */
/* consider a n=3 example
 *     [ 2  0  0 ]   [             ]
 * M = [ 0  0  1 ] = [ m_0 m_1 m_2 ]
 *     [ 0  3  0 ]   [             ]
 *
 * Computing GM shuffles the columns of G and rescales them. Assume G is:
 *     [             ]            [                   ]
 * G = [ g_0 g_1 g_2 ]  then GM = [ 2*g_0 3*g_2 1*g_1 ]
 *     [             ]            [                   ]
 *
 * M is stored in compact form as two arrays, the first one stores the scaling
 * coefficients of the columns (i.e., [2 3 1] in the example), while the second
 * stores, for each element, the position of the column placed at its index
 * after the permutation.
 * In the example we have g_0 -> 0, g_1 -> 2, g_2 -> 1 hence we obtain [0 2 1]
 *
 */


typedef struct {
   /* coefficients listed in order of appearance columnwise */
   FQ_ELEM coefficients[N];
   /* considering the product GQ, permutation[...] stores into the cell with
    * index 0, the position of the DESTINATION of column 0 in G after the
    * computation of GQ.
    */
   POSITION_T permutation[N];
} monomial_t;

typedef struct {
   /* coefficients listed in order of appearance of the colums of the
    * target IS */
   FQ_ELEM coefficients[K];
   /* considering the product GQ, permutation[...] stores into the cell with
    * index 0, the position of the ORIGINAL column in G that goes into 0-th 
    * column of the IS after after the computation of GQ.
    */
   POSITION_T permutation[K];
} monomial_action_IS_t;

typedef struct {
   unsigned char value[SEED_LENGTH_BYTES];
} monomial_seed_t;

/* multiplies two monomial matrices*/
void monomial_mat_mul(monomial_t *res,
                      const monomial_t *const A,
                      const monomial_t *const B);

/* computes the inverse of the monomial matrix */
void monomial_mat_inv(monomial_t *res,
                      const monomial_t *const to_invert);

/* samples a random monomial matrix from the systemwide csprng*/
void monomial_mat_rnd(monomial_t *res);

/* expands a monomial matrix, given a PRNG seed */
void monomial_mat_seed_expand(monomial_t *res,
                              const unsigned char seed[SEED_LENGTH_BYTES]);

/* yields the identity matrix */
void monomial_mat_id(monomial_t *res);

/* composes a compressed action on an IS with the action of a monomial
 * matrix */
void monomial_compose_action(monomial_action_IS_t * out, 
                             const monomial_t * to_compose, 
                             const monomial_action_IS_t * in);

/* Compress MonomialAction object to byte array */
void compress_monom_action(uint8_t *compressed,
                            const monomial_action_IS_t * mono);

/* Decompress byte array to MonomialAction object */
void expand_to_monom_action(monomial_action_IS_t *mono,
                            const uint8_t *compressed);

/* pretty_print for monomial matrices */
void monomial_mat_pretty_print_name(char *name, const monomial_t *to_print);
// void comp_monomial_mat_pretty_print_name(char *name,
//       const compressed_monomial_t *to_print);

/* pretty_print for monomial matrices in their expanded form */
void monomial_mat_print_exp_name(char *name,const monomial_t *to_print);

