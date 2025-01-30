from math import log, log2, ceil

def to_c_array(values, name="foo", formatter=str):
    values = [formatter(v) for v in values]
    body = ', '.join(values)
    return '#define {} {{{}}}'.format(name, body)

def clog2(a):
    return max(int(ceil(log(a,2))), 1)

def l_child(a):
    return 2*a + 1

def r_child(a):
    return 2*a + 2

def tree_offsets_and_nodes(T):

    # Full trees on the left half, so we can already count (i.e. subtract) these values as well as the root node
    missing_nodes_per_level = [2**(i-1) for i in range(1, clog2(T)+1)]
    missing_nodes_per_level.insert(0,0)

    remaining_leaves = T - 2**(clog2(T)-1)
    level = 1

    # Starting from the first level, we construct the tree in a way that the left
    # subtree is always a full binary tree.
    while(remaining_leaves > 0):
        depth = 0
        stree_found = False
        while not stree_found:
            if (remaining_leaves <= 2**depth):
                for i in range(depth, 0, -1):
                    missing_nodes_per_level[level+i] -= 2**(i-1)
                remaining_leaves -= (2**clog2(remaining_leaves)) // 2

                # Subtract root and increase level for next iteration
                missing_nodes_per_level[level] -= 1
                level += 1
                stree_found = True
            else:
                depth += 1

    # The offsets are the missing nodes per level subtracted by the missing nodes of all previous levels, as this
    # is already included
    offsets = [missing_nodes_per_level[i] for i in range(len(missing_nodes_per_level))]
    for i in range(clog2(T), -1, -1):
        for j in range(i):
            offsets[i] -= offsets[j]

    nodes_per_level = [2**i - missing_nodes_per_level[i] for i in range(clog2(T)+1)]
    return offsets, nodes_per_level


def tree_leaves(T, offsets):
    leaves = [0]*T
    leaves_per_level = [0]*(clog2(T)+1)
    start_index_per_level = [0]*(clog2(T)+1)
    ctr = 0

    remaining_leaves = T
    depth = 0
    level = 0
    root_node = 0
    left_child = l_child(root_node) - offsets[level+depth]

    while (remaining_leaves > 0):
        depth = 1
        subtree_found = False
        while not subtree_found:
            if (remaining_leaves <= 2**depth):
                for i in range(2**clog2(remaining_leaves)//2):
                    leaves[ctr] = root_node if remaining_leaves==1 else left_child+i
                    if (remaining_leaves==1):
                        leaves_per_level[level] += 1
                        start_index_per_level[level] = root_node if start_index_per_level[level] == 0 else start_index_per_level[level]
                    else:
                        leaves_per_level[level+depth] += 1
                        start_index_per_level[level+depth] = left_child if start_index_per_level[level+depth] == 0 else start_index_per_level[level+depth]
                    ctr += 1
                root_node = r_child(root_node) - offsets[level]
                left_child = l_child(root_node) - offsets[level]
                level += 1
                remaining_leaves -= 2**clog2(remaining_leaves)//2
                subtree_found = True
            else:
                left_child = l_child(left_child) - offsets[level+depth]
                depth += 1

    # Now create array with start idx and number of leaves by removing zeros
    cons_leaves = [i for i in leaves_per_level if i != 0]
    start_index_per_level = [i for i in start_index_per_level if i != 0]

    return leaves_per_level, len(cons_leaves), start_index_per_level[::-1], cons_leaves[::-1]

def worst_case_tree_cost(t,w):
    num_not_to_reveal = w
    Hamming_weight_t = len(bin(t).strip('0b').replace('0',''))
    return ceil(num_not_to_reveal* log2(t/(num_not_to_reveal))) + Hamming_weight_t - 1


if __name__ == '__main__':

    NIST_cats = [1, 3, 5]
    OPT_CORNER_BY_CAT = {  1: ('BALANCED','INTERMEDIATE','SHORT_SIG'),
                           3: ('BALANCED','SHORT_SIG'),
                           5: ('BALANCED','SHORT_SIG')  }
    T_W_BY_CAT = {1: { 'BALANCED': (192,36),
                       'INTERMEDIATE': (68,42),
                       'SHORT_SIG': (45,34)},
                  3: { 'BALANCED': (220,68),
                       'SHORT_SIG': (102,61)},
                  5: { 'BALANCED': (345,75),
                       'SHORT_SIG': (137,79)}
                  }

    for cat in NIST_cats:
        print(f'Category {cat}')
        for corner in OPT_CORNER_BY_CAT[cat]:
            print(f'Corner {corner}')
            t,w = T_W_BY_CAT[cat][corner]
            print('--------------------------------------')
            offsets, nodes_per_level = tree_offsets_and_nodes(t)
            leaves_per_level, subroots, start_idx, cons_leaves = tree_leaves(t, offsets)
            nodes_to_store = worst_case_tree_cost(t,w)
            print(to_c_array(offsets, 'TREE_OFFSETS'))
            print(to_c_array(nodes_per_level, 'TREE_NODES_PER_LEVEL'))
            print(to_c_array(leaves_per_level, 'TREE_LEAVES_PER_LEVEL'))
            print(f'#define TREE_SUBROOTS {subroots}')
            print(to_c_array(start_idx, 'TREE_LEAVES_START_INDICES'))
            print(to_c_array(cons_leaves, 'TREE_CONSECUTIVE_LEAVES'))
            print(f'#define MAX_PUBLISHED_SEEDS {nodes_to_store}')
            print('--------------------------------------\n')
