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

#include <stdint.h>
#include <string.h>

#include "parameters.h"

#include "seedtree.h"
#include "sha3.h"
#include "rng.h"

#define LEFT_CHILD(i) (2*(i)+1)
#define RIGHT_CHILD(i) (2*(i)+2)
#define PARENT(i) ( ((i)%2) ? (((i)-1)/2) : (((i)-2)/2) )
#define SIBLING(i) ( ((i)%2) ? (i)+1 : (i)-1 )


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
void BuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
              const unsigned char root_seed[SEED_LENGTH_BYTES],
              const unsigned char salt[HASH_DIGEST_LENGTH]) {
    /* input buffer to the CSPRNG, contains the seed to be expanded, a salt,
     * and the integer index of the node being expanded for domain separation */
    const uint32_t csprng_input_len = SALT_LENGTH_BYTES +
        SEED_LENGTH_BYTES;
    unsigned char csprng_input[csprng_input_len];
    SHAKE_STATE_STRUCT tree_csprng_state;
    memcpy(csprng_input+SEED_LENGTH_BYTES, salt, SALT_LENGTH_BYTES);

    /* Set the root seed in the tree from the received parameter */
    memcpy(seed_tree,root_seed,SEED_LENGTH_BYTES);

    /* off contains the offsets required to move between two layers in order
     * to compensate for the truncation.
     * npl contains the number of nodes per level.
     * lpl contains the number of leaves per level
     * */
    const uint16_t off[LOG2(T)+1] = TREE_OFFSETS;
    const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;
    const uint16_t lpl[LOG2(T)+1] = TREE_LEAVES_PER_LEVEL;

    /* Generate the log_2(t) layers from the root, each iteration generates a tree
     * level; iterate on nodes of the parent level; the leaf nodes on each level
     * don't need to be expanded, thus only iterate to npl[level]-lpl[level] */
    int start_node = 0;
    for (int level = 0; level < LOG2(T); level++){
        for (int node_in_level = 0; node_in_level < npl[level]-lpl[level]; node_in_level++ ) {
            uint16_t father_node = start_node + node_in_level;
            uint16_t left_child_node = LEFT_CHILD(father_node) - off[level];

            /* prepare the CSPRNG input to expand the father node */
            memcpy(csprng_input,
                    seed_tree + father_node*SEED_LENGTH_BYTES,
                    SEED_LENGTH_BYTES);

            /* Domain separation using father node index */
            uint16_t domain_sep = father_node;

            /* Generate the children (stored contiguously).
             * By construction, the tree has always two children */
            initialize_csprng_ds(&tree_csprng_state, csprng_input, csprng_input_len, domain_sep);
            csprng_randombytes(seed_tree + left_child_node*SEED_LENGTH_BYTES,
                    2*SEED_LENGTH_BYTES,
                    &tree_csprng_state);
        }
        start_node += npl[level];
    }
}

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

static
void label_leaves(unsigned char flag_tree[NUM_NODES_SEED_TREE],
                     const unsigned char indices_to_publish[T])
{
    const uint16_t cons_leaves[TREE_SUBROOTS] = TREE_CONSECUTIVE_LEAVES;
    const uint16_t leaves_start_indices[TREE_SUBROOTS] = TREE_LEAVES_START_INDICES;

    unsigned int cnt = 0;
    for (size_t i=0; i<TREE_SUBROOTS; i++) {
        for (size_t j=0; j<cons_leaves[i]; j++) {
            flag_tree[leaves_start_indices[i]+j] = indices_to_publish[cnt];
            cnt++;
        }
    }
}

static void compute_seeds_to_publish(
    /* linearized binary tree of boolean nodes containing
     * flags for each node 1-filled nodes are to be released */
    unsigned char flags_tree_to_publish[NUM_NODES_SEED_TREE],
    /* Boolean Array indicating which of the T seeds must be
     * released convention as per the above defines */
    const unsigned char indices_to_publish[T]) {
    /* the indices to publish may be less than the full leaves, copy them
     * into the linearized tree leaves */
    label_leaves(flags_tree_to_publish, indices_to_publish);

    const uint16_t off[LOG2(T)+1] = TREE_OFFSETS;
    const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;
    const uint16_t leaves_start_indices[TREE_SUBROOTS] = TREE_LEAVES_START_INDICES;

    /* compute the value for the internal nodes of the tree starting from
     * the leaves, right to left */
    unsigned int start_node = leaves_start_indices[0];
    for (int level=LOG2(T); level>0; level--) {
        for (int i=npl[level]-2; i>=0; i-=2) {
            uint16_t current_node = start_node + i;
            uint16_t parent_node = PARENT(current_node) + (off[level-1] >> 1);
            if ((flags_tree_to_publish[current_node] == TO_PUBLISH) &&
                (flags_tree_to_publish[SIBLING(current_node)] == TO_PUBLISH)){
                 flags_tree_to_publish[parent_node] = TO_PUBLISH;
            } else {
                 flags_tree_to_publish[parent_node] = NOT_TO_PUBLISH;
            }
        }
        start_node -= npl[level-1];
    }
}

/*****************************************************************************/

uint32_t GGMPath(const unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                 // INPUT: binary array storing in each cell a binary value (i.e., 0 or 1),
                 //        which in turn denotes if the seed of the node with the same index
                 //        must be released (i.e., cell == 0) or not (i.e., cell == 1).
                 //        Indeed, the seed will be stored in the sequence computed as a result into the out[...] array.
                 const unsigned char indices_to_publish[T], // INPUT: binary array denoting which node has to be released (cell == 0) or not
                 unsigned char *seed_storage)             // OUTPUT: sequence of seeds to be released
{
    /* complete linearized binary tree containing boolean values determining
     * if a node is to be released or not according to convention above.
     * */
    unsigned char flags_tree_to_publish[NUM_NODES_SEED_TREE] = {NOT_TO_PUBLISH};
    compute_seeds_to_publish(flags_tree_to_publish, indices_to_publish);

    const uint16_t off[LOG2(T)+1] = TREE_OFFSETS;
    const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;

    /* no sense in trying to publish the root node, start examining from level 1 */
    int start_node = 1;
    int num_seeds_published = 0;

    for (int level = 1; level <= LOG2(T); level++){
        for (int node_in_level = 0; node_in_level < npl[level]; node_in_level++ ) {
            uint16_t current_node = start_node + node_in_level;
            uint16_t father_node = PARENT(current_node) + (off[level-1] >> 1);

            /* if seed is to be published and its ancestor/parent node is not,
             * add it to the seed storage */
            if ( (flags_tree_to_publish[current_node] == TO_PUBLISH) &&
                 (flags_tree_to_publish[father_node] == NOT_TO_PUBLISH) ) {
                memcpy(seed_storage + num_seeds_published*SEED_LENGTH_BYTES,
                        seed_tree + current_node*SEED_LENGTH_BYTES,
                        SEED_LENGTH_BYTES);
                num_seeds_published++;
            }
        }
        start_node += npl[level];
    }
   return num_seeds_published;
}

/*****************************************************************************/

// \return 1 on success
//         0 on failure
uint32_t RebuildGGM(unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES],
                    const unsigned char indices_to_publish[T],
                    const unsigned char *stored_seeds,
                    const unsigned char salt[HASH_DIGEST_LENGTH]) {
   /* complete linearized binary tree containing boolean values determining
     * if a node is to be released or not according to aboves convention
     */
    unsigned char flags_tree_to_publish[NUM_NODES_SEED_TREE] = {0};
    compute_seeds_to_publish(flags_tree_to_publish, indices_to_publish);

    const uint32_t csprng_input_len = SALT_LENGTH_BYTES +
        SEED_LENGTH_BYTES;
    unsigned char csprng_input[csprng_input_len];
    SHAKE_STATE_STRUCT tree_csprng_state;

    const uint16_t off[LOG2(T)+1] = TREE_OFFSETS;
    const uint16_t npl[LOG2(T)+1] = TREE_NODES_PER_LEVEL;
    const uint16_t lpl[LOG2(T)+1] = TREE_LEAVES_PER_LEVEL;

    memcpy(csprng_input + SEED_LENGTH_BYTES, salt, SALT_LENGTH_BYTES);

    /* regenerating the seed tree never starts from the root, as it is never
     * disclosed */
    int nodes_used = 0;
    int start_node = 1;
    for (int level = 1; level <= LOG2(T); level++){
        for (int node_in_level = 0; node_in_level < npl[level]; node_in_level++ ) {
            uint16_t current_node = start_node + node_in_level;
            uint16_t father_node = PARENT(current_node) + (off[level-1] >> 1);
            uint16_t left_child = LEFT_CHILD(current_node) - off[level];

            /* if the current node is a seed which was published (thus its father
             * was not), memcpy it in place */
            if ( flags_tree_to_publish[current_node] == TO_PUBLISH ) {
                if ( flags_tree_to_publish[father_node] == NOT_TO_PUBLISH ) {
                    memcpy(seed_tree + current_node*SEED_LENGTH_BYTES,
                            stored_seeds + nodes_used*SEED_LENGTH_BYTES,
                            SEED_LENGTH_BYTES );
                    nodes_used++;
                }
            }

            /* If the current node is published and not a leaf, CSPRNG-expand its children.
             * Since there is no reason of expanding leaves, only iterate to nodes per level (npl)
             * minus leaves per level (lpl) in each level */
            if ( ( flags_tree_to_publish[current_node] == TO_PUBLISH ) && (node_in_level < npl[level]-lpl[level] ) ) {
                /* prepare the CSPRNG input to expand the children of node current_node */
                memcpy(csprng_input,
                        seed_tree + current_node*SEED_LENGTH_BYTES,
                        SEED_LENGTH_BYTES);

                /* Domain separation using father node index */
                uint16_t domain_sep = current_node;

                /* expand the children (stored contiguously), by construction always two children */
                initialize_csprng_ds(&tree_csprng_state, csprng_input, csprng_input_len, domain_sep);
                csprng_randombytes(seed_tree + left_child*SEED_LENGTH_BYTES,
                        2*SEED_LENGTH_BYTES,
                        &tree_csprng_state);
            }
        }
        start_node += npl[level];
    }

    return 1;
}


/*****************************************************************************/
void seed_leaves(unsigned char rounds_seeds[T*SEED_LENGTH_BYTES],
                 unsigned char seed_tree[NUM_NODES_SEED_TREE*SEED_LENGTH_BYTES])
{
    const uint16_t cons_leaves[TREE_SUBROOTS] = TREE_CONSECUTIVE_LEAVES;
    const uint16_t leaves_start_indices[TREE_SUBROOTS] = TREE_LEAVES_START_INDICES;

    unsigned int cnt = 0;
    for (size_t i=0; i<TREE_SUBROOTS; i++) {
        for (size_t j=0; j<cons_leaves[i]; j++) {
            memcpy(rounds_seeds + cnt*SEED_LENGTH_BYTES,
                   seed_tree + (leaves_start_indices[i]+j)*SEED_LENGTH_BYTES,
                   SEED_LENGTH_BYTES);
            cnt++;
        }
    }
}
