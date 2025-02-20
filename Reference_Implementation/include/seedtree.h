/**
 *
 * Reference ISO-C11 Implementation of LESS.
 *
 * @version 1.2 (February 2025)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Patrick Karl <patrick.karl@tum.de> 
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

#include <stdint.h>
#include "parameters.h"

/******************************************************************************/

void BuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
              const unsigned char root_seed[SEED_LENGTH_BYTES],
              const unsigned char salt[HASH_DIGEST_LENGTH]) ;

/******************************************************************************/

/* returns the number of seeds which have been published */
uint32_t GGMPath(const unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                // binary array denoting if node has to be released (cell == 0) or not
                const unsigned char indices_to_publish[T],
                unsigned char *seed_storage);

/******************************************************************************/

/* returns the number of seeds which have been used to regenerate the tree */
uint32_t RebuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                    const unsigned char indices_to_publish[T],
                    const unsigned char *stored_seeds,
                    const unsigned char salt[HASH_DIGEST_LENGTH]);   // input

void seed_leaves(unsigned char rounds_seeds[T*SEED_LENGTH_BYTES],
                 unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES]);
