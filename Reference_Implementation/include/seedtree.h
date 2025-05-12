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

/// \param seed_tree[out]: full seed tree
/// \param root_seed[in]: main seed fed into the prng
/// \param salt[in]: additionally fed into the prng
void BuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
              const unsigned char root_seed[SEED_LENGTH_BYTES],
              const unsigned char salt[HASH_DIGEST_LENGTH]) ;

/// \param seed_tree[in]: full seed tree
/// \param indices_to_publish[in]: binary array denoting if node has to be 
///     released (cell == 0) or not.
/// \param seed_storage[out]: byte array: sequence of seeds to be published
/// \returns the number of seeds which have been published 
uint32_t GGMPath(const unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                 const unsigned char indices_to_publish[T],
                 unsigned char *seed_storage);

/// \param seed_tree[out]: full seed tree
/// \param root_seed[in]: main seed fed into the prng
/// \param indices_to_publish[in]: binary array denoting if node has to be 
///     released (cell == 0) or not.
/// \param salt[in]: additional salt fed into the prng
/// \returns 1 on success, 0 on failure
uint32_t RebuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                    const unsigned char indices_to_publish[T],
                    const unsigned char *stored_seeds,
                    const unsigned char salt[HASH_DIGEST_LENGTH]);

/// \param rounds_seeds[out]: output parameter
/// \param seed_tree[in]: full seed tree
void seed_leaves(unsigned char rounds_seeds[T*SEED_LENGTH_BYTES],
                 const unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES]);
