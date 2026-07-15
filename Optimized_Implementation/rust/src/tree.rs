use sha3::Shake128;
use crate::config::{
    T, T_LOG2,
    SEED_LENGTH_BYTES, SALT_LENGTH_BYTES, NUM_NODES_SEED_TREE,
    TREE_NODES_PER_LEVEL, TREE_LEAVES_PER_LEVEL, TREE_OFFSETS, TREE_SUBROOTS, TREE_LEAVES_START_INDICES, TREE_CONSECUTIVE_LEAVES
};
use crate::prng::{
    csprng_randombytes,
    csprng_initialize_ds
};

const TO_PUBLISH: u8 = 0;
const NOT_TO_PUBLISH: u8 = 1;

#[inline(always)]
const fn left_child(i: usize) -> usize {
    2 * i + 1
}

#[inline(always)]
const fn right_child(i: usize) -> usize {
    2 * i + 2
}

#[inline(always)]
const fn parent(i: usize) -> usize {
    if i & 1 == 1 {
        (i - 1) / 2
    } else {
        (i - 2) / 2
    }
}

#[inline(always)]
const fn sibling(i: usize) -> usize {
    if i & 1 == 1 {
        i + 1
    } else {
        i - 1
    }
}
#[inline(always)]
fn seed_range(node: usize) -> std::ops::Range<usize> {
    let start = node * SEED_LENGTH_BYTES;
    start..start + SEED_LENGTH_BYTES
}

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
pub fn build_ggm(
    seed_tree: &mut [u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
    root_seed: &[u8; SEED_LENGTH_BYTES],
    salt: &[u8; SALT_LENGTH_BYTES],
) {
    let mut csprng_input = [0u8; SEED_LENGTH_BYTES + SALT_LENGTH_BYTES];

    csprng_input[SEED_LENGTH_BYTES..].copy_from_slice(salt);

    seed_tree[..SEED_LENGTH_BYTES].copy_from_slice(root_seed);

    let mut start_node = 0usize;

    for level in 0..T_LOG2 {
        for node_in_level in 0..(TREE_NODES_PER_LEVEL[level] - TREE_LEAVES_PER_LEVEL[level]) {

            let father = start_node + node_in_level;
            let left = left_child(father) - TREE_OFFSETS[level];

            csprng_input[..SEED_LENGTH_BYTES]
                .copy_from_slice(&seed_tree[seed_range(father)]);

            let mut state = csprng_initialize_ds::<Shake128>(&csprng_input, father as u16);

            csprng_randombytes(
                &mut seed_tree[left * SEED_LENGTH_BYTES
                    ..left * SEED_LENGTH_BYTES + 2 * SEED_LENGTH_BYTES],
                &mut state,
            );
        }

        start_node += TREE_NODES_PER_LEVEL[level];
    }
}

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
fn label_leaves(
    flag_tree: &mut [u8; NUM_NODES_SEED_TREE],
    indices: &[u8; T],
) {
    let mut cnt = 0;

    for i in 0..TREE_SUBROOTS {
        for j in 0..TREE_CONSECUTIVE_LEAVES[i] {
            flag_tree[TREE_LEAVES_START_INDICES[i] + j] = indices[cnt];
            cnt += 1;
        }
    }
}
fn compute_seeds_to_publish(
    flags: &mut [u8; NUM_NODES_SEED_TREE],
    indices: &[u8; T],
) {
    label_leaves(flags, indices);

    let mut start = TREE_LEAVES_START_INDICES[0];

    for level in (1..=T_LOG2).rev() {

        let n = TREE_NODES_PER_LEVEL[level];

        let mut i = n - 2;

        loop {
            let current = start + i;

            let parent =
                parent(current) + (TREE_OFFSETS[level - 1] >> 1);

            flags[parent] =
                if flags[current] == TO_PUBLISH
                    && flags[sibling(current)] == TO_PUBLISH
                {
                    TO_PUBLISH
                } else {
                    NOT_TO_PUBLISH
                };

            if i < 2 {
                break;
            }

            i -= 2;
        }

        start -= TREE_NODES_PER_LEVEL[level - 1];
    }
}

pub fn ggm_path(
    seed_tree: &[u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
    indices: &[u8; T],
    out: &mut [u8],
) -> usize {
    let mut flags = [NOT_TO_PUBLISH; NUM_NODES_SEED_TREE];

    compute_seeds_to_publish(&mut flags, indices);

    // TODO out.clear();

    let mut start = 1usize;

    for level in 1..=T_LOG2 {
        for node in 0..TREE_NODES_PER_LEVEL[level] {
            let current = start + node;
            let father = parent(current) + (TREE_OFFSETS[level - 1] >> 1);

            if flags[current] == TO_PUBLISH  && flags[father] == NOT_TO_PUBLISH  {
                out.copy_from_slice(&seed_tree[seed_range(current)]
                );
            }
        }

        start += TREE_NODES_PER_LEVEL[level];
    }

    out.len() / SEED_LENGTH_BYTES
}

pub fn seed_leaves(
    out: &mut [u8; T * SEED_LENGTH_BYTES],
    tree: &[u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
) {
    let mut cnt = 0;

    for i in 0..TREE_SUBROOTS {
        for j in 0..TREE_CONSECUTIVE_LEAVES[i] {

            let node = TREE_LEAVES_START_INDICES[i] + j;

            out[cnt * SEED_LENGTH_BYTES
                ..(cnt + 1) * SEED_LENGTH_BYTES]
                .copy_from_slice(&tree[seed_range(node)]);

            cnt += 1;
        }
    }
}