use crate::config::{
    T, T_LOG2, Xof,
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
#[allow(dead_code)]
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

            let mut state = csprng_initialize_ds::<Xof>(&csprng_input, father as u16);

            seed_tree[left * SEED_LENGTH_BYTES..left * SEED_LENGTH_BYTES + 2 * SEED_LENGTH_BYTES]
                .copy_from_slice(&csprng_randombytes(&mut state, 2 * SEED_LENGTH_BYTES));
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

    out.fill(0);

    /* no sense in trying to publish the root node, start examining from level 1 */
    let mut start = 1usize;
    let mut num_seeds_published = 0usize;

    for level in 1..=T_LOG2 {
        for node in 0..TREE_NODES_PER_LEVEL[level] {
            let current = start + node;
            let father = parent(current) + (TREE_OFFSETS[level - 1] >> 1);

            /* if seed is to be published and its ancestor/parent node is not,
             * add it to the seed storage */
            if flags[current] == TO_PUBLISH  && flags[father] == NOT_TO_PUBLISH  {
                let dst = num_seeds_published * SEED_LENGTH_BYTES;
                out[dst..dst + SEED_LENGTH_BYTES]
                    .copy_from_slice(&seed_tree[seed_range(current)]);
                num_seeds_published += 1;
            }
        }

        start += TREE_NODES_PER_LEVEL[level];
    }

    num_seeds_published
}

/// regenerates the seed tree from the seeds which were published by
/// `ggm_path`. The root is never disclosed, so the rebuild starts at level 1.
/// \param seed_tree[out]: the rebuilt tree
/// \param indices[in]: binary array, 0 for the leaves whose seed was released
/// \param stored_seeds[in]: the published seeds
/// \param salt[in]: salt which was used to build the tree
/// \return 1 on success
///         0 on failure
pub fn rebuild_ggm(
    seed_tree: &mut [u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES],
    indices: &[u8; T],
    stored_seeds: &[u8],
    salt: &[u8; SALT_LENGTH_BYTES],
) -> u32 {
    let mut flags = [0u8; NUM_NODES_SEED_TREE];
    compute_seeds_to_publish(&mut flags, indices);

    let mut csprng_input = [0u8; SEED_LENGTH_BYTES + SALT_LENGTH_BYTES];
    csprng_input[SEED_LENGTH_BYTES..].copy_from_slice(salt);

    let mut nodes_used = 0usize;
    let mut start = 1usize;

    for level in 1..=T_LOG2 {
        for node_in_level in 0..TREE_NODES_PER_LEVEL[level] {
            let current = start + node_in_level;
            let father = parent(current) + (TREE_OFFSETS[level - 1] >> 1);
            let left = left_child(current) - TREE_OFFSETS[level];

            /* if the current node is a seed which was published (thus its father
             * was not), memcpy it in place */
            if flags[current] == TO_PUBLISH && flags[father] == NOT_TO_PUBLISH {
                let src = nodes_used * SEED_LENGTH_BYTES;
                if src + SEED_LENGTH_BYTES > stored_seeds.len() {
                    return 0;
                }
                let s: [u8; SEED_LENGTH_BYTES] =
                    stored_seeds[src..src + SEED_LENGTH_BYTES].try_into().unwrap();
                seed_tree[seed_range(current)].copy_from_slice(&s);
                nodes_used += 1;
            }

            /* If the current node is published and not a leaf, CSPRNG-expand its
             * children. */
            if flags[current] == TO_PUBLISH
                && node_in_level < TREE_NODES_PER_LEVEL[level] - TREE_LEAVES_PER_LEVEL[level]
            {
                csprng_input[..SEED_LENGTH_BYTES]
                    .copy_from_slice(&seed_tree[seed_range(current)]);

                let mut state = csprng_initialize_ds::<Xof>(&csprng_input, current as u16);

                seed_tree[left * SEED_LENGTH_BYTES..left * SEED_LENGTH_BYTES + 2 * SEED_LENGTH_BYTES]
                    .copy_from_slice(&csprng_randombytes(&mut state, 2 * SEED_LENGTH_BYTES));
            }
        }

        start += TREE_NODES_PER_LEVEL[level];
    }

    1
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{MAX_PUBLISHED_SEEDS, SEED_TREE_MAX_PUBLISHED_BYTES, W};

    fn build() -> [u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES] {
        let mut tree = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
        let root = [0x42u8; SEED_LENGTH_BYTES];
        let salt = [0x17u8; SALT_LENGTH_BYTES];
        build_ggm(&mut tree, &root, &salt);
        tree
    }

    #[test]
    fn node_indices() {
        for i in 1..64usize {
            assert_eq!(parent(left_child(i)), i);
            assert_eq!(parent(right_child(i)), i);
            assert_eq!(sibling(left_child(i)), right_child(i));
            assert_eq!(sibling(right_child(i)), left_child(i));
        }
        assert_eq!(left_child(0), 1);
        assert_eq!(right_child(0), 2);
        assert_eq!(parent(1), 0);
        assert_eq!(parent(2), 0);
    }

    #[test]
    fn seed_range_is_one_seed() {
        let r = seed_range(0);
        assert_eq!(r.start, 0);
        assert_eq!(r.end, SEED_LENGTH_BYTES);
        let r = seed_range(5);
        assert_eq!(r.start, 5 * SEED_LENGTH_BYTES);
        assert_eq!(r.len(), SEED_LENGTH_BYTES);
    }

    #[test]
    fn build_ggm_root_and_expansion() {
        let tree = build();
        assert_eq!(&tree[..SEED_LENGTH_BYTES], &[0x42u8; SEED_LENGTH_BYTES]);

        // every node that is reachable from the root must have been written
        let mut zero_nodes = 0;
        for node in 0..NUM_NODES_SEED_TREE {
            if tree[seed_range(node)].iter().all(|&b| b == 0) {
                zero_nodes += 1;
            }
        }
        assert_eq!(zero_nodes, 0, "{} nodes were never expanded", zero_nodes);
    }

    #[test]
    fn build_ggm_is_deterministic_and_salted() {
        let a = build();

        let mut b = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
        build_ggm(&mut b, &[0x42u8; SEED_LENGTH_BYTES], &[0x17u8; SALT_LENGTH_BYTES]);
        assert_eq!(a, b);

        let mut c = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
        build_ggm(&mut c, &[0x42u8; SEED_LENGTH_BYTES], &[0x18u8; SALT_LENGTH_BYTES]);
        assert_ne!(a[SEED_LENGTH_BYTES..], c[SEED_LENGTH_BYTES..]);

        let mut d = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
        build_ggm(&mut d, &[0x43u8; SEED_LENGTH_BYTES], &[0x17u8; SALT_LENGTH_BYTES]);
        assert_ne!(a[SEED_LENGTH_BYTES..], d[SEED_LENGTH_BYTES..]);
    }

    #[test]
    fn seed_leaves_collects_t_distinct_seeds() {
        let tree = build();
        let mut leaves = [0u8; T * SEED_LENGTH_BYTES];
        seed_leaves(&mut leaves, &tree);

        let mut seen = std::collections::HashSet::new();
        for i in 0..T {
            let s = &leaves[i * SEED_LENGTH_BYTES..(i + 1) * SEED_LENGTH_BYTES];
            assert!(s.iter().any(|&b| b != 0));
            assert!(seen.insert(s.to_vec()), "leaf {} is a duplicate", i);
        }
    }

    #[test]
    fn label_leaves_places_every_index() {
        let mut flags = [0xffu8; NUM_NODES_SEED_TREE];
        let indices: [u8; T] = core::array::from_fn(|i| (i % 2) as u8);
        label_leaves(&mut flags, &indices);

        let mut cnt = 0;
        for i in 0..TREE_SUBROOTS {
            for j in 0..TREE_CONSECUTIVE_LEAVES[i] {
                assert_eq!(flags[TREE_LEAVES_START_INDICES[i] + j], indices[cnt]);
                cnt += 1;
            }
        }
        assert_eq!(cnt, T);
    }

    #[test]
    fn compute_seeds_to_publish_all_or_nothing() {
        // every leaf published => the root is published
        let mut flags = [NOT_TO_PUBLISH; NUM_NODES_SEED_TREE];
        compute_seeds_to_publish(&mut flags, &[TO_PUBLISH; T]);
        assert_eq!(flags[0], TO_PUBLISH);

        // no leaf published => the root is not published
        let mut flags = [NOT_TO_PUBLISH; NUM_NODES_SEED_TREE];
        compute_seeds_to_publish(&mut flags, &[NOT_TO_PUBLISH; T]);
        assert_eq!(flags[0], NOT_TO_PUBLISH);
    }

    #[test]
    fn ggm_path_publishes_nothing_when_all_leaves_are_kept() {
        let tree = build();
        let mut out = [0u8; SEED_TREE_MAX_PUBLISHED_BYTES];
        let n = ggm_path(&tree, &[NOT_TO_PUBLISH; T], &mut out);
        assert_eq!(n, 0);
        assert!(out.iter().all(|&b| b == 0));
    }

    #[test]
    fn ggm_path_publishes_nothing_when_all_leaves_are_published() {
        let tree = build();
        let mut out = [0u8; SEED_TREE_MAX_PUBLISHED_BYTES];
        let n = ggm_path(&tree, &[TO_PUBLISH; T], &mut out);

        // the whole tree collapses into the root, and the root is never
        // emitted, so there is nothing to publish
        assert_eq!(n, 0);
        assert!(out.iter().all(|&b| b == 0));
    }

    /// keeping a single leaf secret costs one seed per level on the path from
    /// that leaf to the root
    #[test]
    fn ggm_path_publishes_the_siblings_of_one_kept_leaf() {
        let tree = build();

        let mut indices = [TO_PUBLISH; T];
        indices[0] = NOT_TO_PUBLISH;

        let mut out = [0u8; SEED_TREE_MAX_PUBLISHED_BYTES];
        let n = ggm_path(&tree, &indices, &mut out);

        assert!(n >= 1);
        assert!(n <= T_LOG2, "{} > T_LOG2", n);
        assert!(n <= MAX_PUBLISHED_SEEDS);

        // every emitted seed is a real node of the tree
        for i in 0..n {
            let s = &out[i * SEED_LENGTH_BYTES..(i + 1) * SEED_LENGTH_BYTES];
            assert!(s.iter().any(|&b| b != 0));
        }
        // and nothing past the reported count was written
        assert!(out[n * SEED_LENGTH_BYTES..].iter().all(|&b| b == 0));
    }

    #[test]
    fn ggm_path_respects_max_published_seeds() {
        let tree = build();

        // a challenge of weight W, i.e. W leaves kept and T-W published
        let mut indices = [TO_PUBLISH; T];
        for i in 0..W {
            indices[(i * T) / W] = NOT_TO_PUBLISH;
        }

        let mut out = [0u8; SEED_TREE_MAX_PUBLISHED_BYTES];
        let n = ggm_path(&tree, &indices, &mut out);
        assert!(n > 0);
        assert!(n <= MAX_PUBLISHED_SEEDS, "{} > MAX_PUBLISHED_SEEDS", n);

        // every published seed must be a real node of the tree
        for i in 0..n {
            let s = &out[i * SEED_LENGTH_BYTES..(i + 1) * SEED_LENGTH_BYTES];
            assert!(s.iter().any(|&b| b != 0));
        }
        // and nothing past the reported count was touched
        assert!(out[n * SEED_LENGTH_BYTES..].iter().all(|&b| b == 0));
    }
}