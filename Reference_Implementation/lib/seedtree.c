#include "seedtree.h"
#include <stdio.h>

#define LEFT_CHILD(i)  (2u*(uint32_t)(i)+1u)
#define RIGHT_CHILD(i) (2u*(uint32_t)(i)+2u)
#define PARENT(i) ((i-1)/2)

/* Seed tree implementation. The binary seed tree is linearized into an array
 * from root to leaves, and from left to right */

/**
 * unsigned char *seed_tree:
 * it is intended as an output parameter;
 * it is an array of uchars that is going to store a sequence of SEED_LENGTH_BYTES,
 * with length: 2*leaves-1.
 *
 *
 * The root seed is taken as a parameter.
 * The seed of its TWO children are computed expanding (i.e., shake128...) the
 * entropy in "salt" + "seedBytes of the parent" associated to each node
 *             from roots to leaves layer-by-layer from left to right,
 *             counting from 0 (the integer bound with the root node)"
 *
 */
void generate_seed_tree_from_root(unsigned char
                                  seed_tree[NUM_NODES_OF_SEED_TREE * SEED_LENGTH_BYTES],
                                  const unsigned char root_seed[SEED_LENGTH_BYTES],
                                  const unsigned char salt[HASH_DIGEST_LENGTH]) {
   /* input buffer to the CSPRNG, contains a salt, the seed to be expanded
    * and the integer index of the node being expanded for domain separation */
   const uint32_t csprng_input_len = HASH_DIGEST_LENGTH +
                                     SEED_LENGTH_BYTES;
   unsigned char csprng_input[csprng_input_len];
   SHAKE_STATE_STRUCT tree_csprng_state;

   memcpy(csprng_input, salt, HASH_DIGEST_LENGTH);

   /* Set the root seed in the tree from the received parameter */
   memcpy(seed_tree, root_seed, SEED_LENGTH_BYTES);
   for (uint32_t i = 0; i < NUM_LEAVES_OF_SEED_TREE-1; i++) {
      /* prepare the CSPRNG input to expand the children of node i */
      memcpy(csprng_input + HASH_DIGEST_LENGTH,
             seed_tree + i*SEED_LENGTH_BYTES,
             SEED_LENGTH_BYTES);

      /* expand the children (stored contiguously) */
      initialize_csprng(&tree_csprng_state, csprng_input, csprng_input_len);
      csprng_randombytes(seed_tree + LEFT_CHILD(i)*SEED_LENGTH_BYTES,
                         2*SEED_LENGTH_BYTES, &tree_csprng_state);
   }
} /* end generate_seed_tree */

/*****************************************************************************/

/**
 * const unsigned char *indices: input parameter denoting an array
 * with a number of binary cells equal to "leaves" representing
 * the labels of the nodes identified as leaves of the tree[...]
 * passed as second parameter.
 * A label = 0 means that the byteseed of the node having the same index
 * has to be released; = 1, otherwise.
 *
 * unsigned char *tree: input/output parameter denoting an array
 * with a number of binary cells equal to "2*leaves-1";
 * the first "leaves" cells (i.e., the ones with positions from 0 to leaves-1)
 * are the ones that will be modified by the current subroutine,
 * the last "leaves" cells will be a copy of the input array passed as first
 * parameter.
 *
 * uint64_t leaves: input parameter;
 *
 */

#define TO_PUBLISH 0
#define NOT_TO_PUBLISH 1

static void compute_seeds_to_publish(
   /* OUTPUT:
	* linearized binary tree of boolean nodes containing
    * flags for each node.
	* NOTE:1-filled nodes are not to be released */
   unsigned char flags_tree_to_publish[NUM_NODES_OF_SEED_TREE],
   /* INPUT:
	* Boolean Array indicating which of the T seeds must be
    * released convention as per the above defines */
   const unsigned char indices_to_publish[T]) {
   /* the indices to publish may be less than the full leaves, copy them
    * into the linearized tree leaves */
   memcpy(flags_tree_to_publish + NUM_LEAVES_OF_SEED_TREE-1, 
		  indices_to_publish, T);

   /* NOTE: do not change the `int`.
    * compute the value for the internal nodes of the tree starting from the
    * fathers of the leaves, right to left */
   for (int i = NUM_LEAVES_OF_SEED_TREE-2; i >= 0; i--) {
      flags_tree_to_publish[i] = ( flags_tree_to_publish[LEFT_CHILD(i)]  == NOT_TO_PUBLISH) ||
                                 ( flags_tree_to_publish[RIGHT_CHILD(i)] == NOT_TO_PUBLISH);
   }
} /* end compute_seeds_to_publish */

/*****************************************************************************/

uint32_t seed_tree_path(const unsigned char
                  seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES],
                  // INPUT: binary array storing in each cell a binary value (i.e., 0 or 1),
                  //        which in turn denotes if the seed of the node with the same index
                  //        must be released (i.e., cell == 0) or not (i.e., cell == 1).
                  //        Indeed the seed will be stored in the sequence computed as a result into the out[...] array.
                  const unsigned char
                  indices_to_publish[T], // INPUT: binary array denoting which node has to be released (cell == 0) or not
                  unsigned char
                  *seed_storage)             // OUTPUT: sequence of seeds to be released
{
   /* complete linearized binary tree containing boolean values determining
    * if a node is to be released or not. Nodes set to 1 are not to be released
    * oldest ancestor of sets of nodes equal to 0 are to be released */
   unsigned char flags_tree_to_publish[NUM_NODES_OF_SEED_TREE] = {0};
   compute_seeds_to_publish(flags_tree_to_publish, indices_to_publish);

   uint32_t num_seeds_published = 0;
   for (uint32_t i = 0; i < 2*NUM_LEAVES_OF_SEED_TREE-1; i++) {
      if (flags_tree_to_publish[i]         == TO_PUBLISH &&
          flags_tree_to_publish[PARENT(i)] == NOT_TO_PUBLISH ) {
         memcpy(seed_storage + num_seeds_published*SEED_LENGTH_BYTES,
                seed_tree + i*SEED_LENGTH_BYTES,
                SEED_LENGTH_BYTES );
         num_seeds_published++;
      }
   }

   return num_seeds_published;
} /* end seed_tree_path */

/*****************************************************************************/

uint32_t rebuild_seed_tree_leaves(unsigned char
                      seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES],
                      const unsigned char indices_to_publish[T],
                      const unsigned char *stored_seeds,
                      const unsigned char salt[HASH_DIGEST_LENGTH]) {
   /* complete linearized binary tree containing boolean values determining
    * if a node is to be released or not. Nodes set to 1 are not to be released
    * oldest ancestor of sets of nodes equal to 0 are to be released */
   unsigned char flags_tree_to_publish[2*NUM_LEAVES_OF_SEED_TREE-1] = {0};
   compute_seeds_to_publish(flags_tree_to_publish, indices_to_publish);

   const uint32_t csprng_input_len = HASH_DIGEST_LENGTH +
                                     SEED_LENGTH_BYTES;
   unsigned char csprng_input[csprng_input_len];
   SHAKE_STATE_STRUCT tree_csprng_state;

   memcpy(csprng_input, salt, HASH_DIGEST_LENGTH);

   uint32_t nodes_used = 0;
   for (uint32_t i = 0; i < NUM_NODES_OF_SEED_TREE; i++) {
      /* if the current node is a seed which was published, memcpy it in place */
      if ( flags_tree_to_publish[i] == TO_PUBLISH ) {
         if ( flags_tree_to_publish[PARENT(i)] == NOT_TO_PUBLISH ) {
            memcpy(seed_tree + SEED_LENGTH_BYTES*i,
                   stored_seeds + SEED_LENGTH_BYTES*nodes_used,
                   SEED_LENGTH_BYTES );
            nodes_used++;
         }
         /* if the current node is not a leaf, CSPRNG-expand its children */
         if ( i < NUM_LEAVES_OF_SEED_TREE-1 ) {
            /* prepare the CSPRNG input to expand the children of node i */
            memcpy(csprng_input + HASH_DIGEST_LENGTH,
                   seed_tree + i*SEED_LENGTH_BYTES,
                   SEED_LENGTH_BYTES);
            initialize_csprng(&tree_csprng_state, csprng_input, csprng_input_len);
            csprng_randombytes(seed_tree + LEFT_CHILD(i)*SEED_LENGTH_BYTES,
                               2*SEED_LENGTH_BYTES, &tree_csprng_state);
         }
      }
   }

   return nodes_used;
} /* end rebuild_seed_tree_leaves */
