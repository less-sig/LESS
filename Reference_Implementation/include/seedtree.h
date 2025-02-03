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
