#pragma once

#include <stdint.h>
#include <string.h> // memcpy(...), memset(...)
#include "parameters.h"
#include "rng.h"

/******************************************************************************/

void generate_seed_tree_from_root(unsigned char
                                  seed_tree[NUM_NODES_OF_SEED_TREE *
                                                               SEED_LENGTH_BYTES],
                                  const unsigned char root_seed[SEED_LENGTH_BYTES],
                                  const unsigned char salt[HASH_DIGEST_LENGTH]) ;

/******************************************************************************/

/* returns the number of seeds which have been published */
uint32_t seed_tree_path(const unsigned char
                  seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES],
                  // binary array denoting if node has to be released (cell == 0) or not
                  const unsigned char indices_to_publish[T],
                  unsigned char *seed_storage);

/******************************************************************************/

/* returns the number of seeds which have been used to regenerate the tree */
uint32_t rebuild_seed_tree_leaves(unsigned char
                      seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES],
                      const unsigned char indices_to_publish[T],
                      const unsigned char *stored_seeds,
                      const unsigned char salt[HASH_DIGEST_LENGTH]);   // input
