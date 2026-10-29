#!/usr/bin/env python3
from config import (
    T, T_LOG2, SEED_LENGTH_BYTES, NUM_NODES_SEED_TREE, TREE_NODES_PER_LEVEL,
    TREE_LEAVES_PER_LEVEL, TREE_OFFSETS, TREE_SUBROOTS, TREE_LEAVES_START_INDICES,
    TREE_CONSECUTIVE_LEAVES,
)
from prng import csprng_initialize_ds, csprng_randombytes

TO_PUBLISH = 0
NOT_TO_PUBLISH = 1


def left_child(i: int) -> int:
    return 2 * i + 1


def right_child(i: int) -> int:
    return 2 * i + 2


def parent(i: int) -> int:
    return (i - 1) // 2 if i & 1 else (i - 2) // 2


def sibling(i: int) -> int:
    return i + 1 if i & 1 else i - 1


def seed_range(node: int) -> slice:
    return slice(node * SEED_LENGTH_BYTES, (node + 1) * SEED_LENGTH_BYTES)


def expand_children(seed_tree: bytearray, node: int, left: int, salt: bytes):
    state = csprng_initialize_ds(bytes(seed_tree[seed_range(node)]) + bytes(salt), node)
    seed_tree[left * SEED_LENGTH_BYTES:(left + 2) * SEED_LENGTH_BYTES] = \
        csprng_randombytes(state, 2 * SEED_LENGTH_BYTES)


def build_ggm(root_seed: bytes, salt: bytes) -> bytearray:
    """
    Seed tree implementation. The binary seed tree is linearized into an array
    from root to leaves, and from left to right

    The root seed is taken as a parameter.
    The seed of its TWO children are computed expanding (i.e., shake128...) the
    entropy in "salt" + "seedBytes of the parent" associated to each node
                from roots to leaves layer-by-layer from left to right,
                counting from 0 (the integer bound with the root node)"
    """
    seed_tree = bytearray(NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES)
    seed_tree[seed_range(0)] = root_seed
    start_node = 0
    for level in range(T_LOG2):
        for node_in_level in range(TREE_NODES_PER_LEVEL[level] - TREE_LEAVES_PER_LEVEL[level]):
            father = start_node + node_in_level
            expand_children(seed_tree, father, left_child(father) - TREE_OFFSETS[level], salt)
        start_node += TREE_NODES_PER_LEVEL[level]
    return seed_tree


def label_leaves(flag_tree: list[int], indices: list[int]):
    cnt = 0
    for i in range(TREE_SUBROOTS):
        for j in range(TREE_CONSECUTIVE_LEAVES[i]):
            flag_tree[TREE_LEAVES_START_INDICES[i] + j] = indices[cnt]
            cnt += 1


def compute_seeds_to_publish(indices: list[int], init: int) -> list[int]:
    flags = [init] * NUM_NODES_SEED_TREE
    label_leaves(flags, indices)
    start = TREE_LEAVES_START_INDICES[0]
    for level in range(T_LOG2, 0, -1):
        for i in range(TREE_NODES_PER_LEVEL[level] - 2, -1, -2):
            current = start + i
            father = parent(current) + (TREE_OFFSETS[level - 1] >> 1)
            both = flags[current] == TO_PUBLISH and flags[sibling(current)] == TO_PUBLISH
            flags[father] = TO_PUBLISH if both else NOT_TO_PUBLISH
        start -= TREE_NODES_PER_LEVEL[level - 1]
    return flags


def ggm_path(seed_tree: bytearray, indices: list[int]) -> bytes:
    flags = compute_seeds_to_publish(indices, NOT_TO_PUBLISH)
    out = bytearray()

    # no sense in trying to publish the root node, start examining from level 1
    start = 1
    for level in range(1, T_LOG2 + 1):
        for node in range(TREE_NODES_PER_LEVEL[level]):
            current = start + node
            father = parent(current) + (TREE_OFFSETS[level - 1] >> 1)

            # if seed is to be published and its ancestor/parent node is not,
            # add it to the seed storage
            if flags[current] == TO_PUBLISH and flags[father] == NOT_TO_PUBLISH:
                out += seed_tree[seed_range(current)]
        start += TREE_NODES_PER_LEVEL[level]
    return bytes(out)


def rebuild_ggm(indices: list[int], stored_seeds: bytes, salt: bytes):
    """
    regenerates the seed tree from the seeds which were published by
    `ggm_path`. The root is never disclosed, so the rebuild starts at level 1.
    \\param indices[in]: binary array, 0 for the leaves whose seed was released
    \\param stored_seeds[in]: the published seeds
    \\param salt[in]: salt which was used to build the tree
    \\return the rebuilt tree on success
             None on failure
    """
    flags = compute_seeds_to_publish(indices, TO_PUBLISH)
    seed_tree = bytearray(NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES)
    nodes_used = 0
    start = 1
    for level in range(1, T_LOG2 + 1):
        for node_in_level in range(TREE_NODES_PER_LEVEL[level]):
            current = start + node_in_level
            father = parent(current) + (TREE_OFFSETS[level - 1] >> 1)

            # if the current node is a seed which was published (thus its father
            # was not), memcpy it in place
            if flags[current] == TO_PUBLISH and flags[father] == NOT_TO_PUBLISH:
                src = seed_range(nodes_used)
                if src.stop > len(stored_seeds):
                    return None
                seed_tree[seed_range(current)] = stored_seeds[src]
                nodes_used += 1

            # If the current node is published and not a leaf, CSPRNG-expand its
            # children.
            if flags[current] == TO_PUBLISH and \
                    node_in_level < TREE_NODES_PER_LEVEL[level] - TREE_LEAVES_PER_LEVEL[level]:
                expand_children(seed_tree, current, left_child(current) - TREE_OFFSETS[level], salt)
        start += TREE_NODES_PER_LEVEL[level]
    return seed_tree


def seed_leaves(tree: bytearray) -> list[bytes]:
    return [bytes(tree[seed_range(TREE_LEAVES_START_INDICES[i] + j)])
            for i in range(TREE_SUBROOTS) for j in range(TREE_CONSECUTIVE_LEAVES[i])]


def test():
    for i in range(1, 64):
        assert parent(left_child(i)) == i and parent(right_child(i)) == i
        assert sibling(left_child(i)) == right_child(i)

    tree = build_ggm(bytes([0x42] * SEED_LENGTH_BYTES), bytes([0x17] * 2 * SEED_LENGTH_BYTES))
    # every node that is reachable from the root must have been written
    assert all(any(tree[seed_range(node)]) for node in range(NUM_NODES_SEED_TREE))
    assert len(set(seed_leaves(tree))) == T

    assert ggm_path(tree, [NOT_TO_PUBLISH] * T) == b""
    # the whole tree collapses into the root, and the root is never
    # emitted, so there is nothing to publish
    assert ggm_path(tree, [TO_PUBLISH] * T) == b""

    # keeping a single leaf secret costs one seed per level on the path from
    # that leaf to the root
    indices = [TO_PUBLISH] * T
    indices[0] = NOT_TO_PUBLISH
    out = ggm_path(tree, indices)
    assert 1 <= len(out) // SEED_LENGTH_BYTES <= T_LOG2

    rebuilt = rebuild_ggm(indices, out, bytes([0x17] * 2 * SEED_LENGTH_BYTES))
    assert seed_leaves(rebuilt)[1:] == seed_leaves(tree)[1:]


if __name__ == "__main__":
    test()
