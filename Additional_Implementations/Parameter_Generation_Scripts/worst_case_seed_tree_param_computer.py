#!/usr/bin/python3
from math import log2, ceil
from sys import argv,exit

# This function computes the worst-case seed tree cost by operatively
# determining the number of seeds to be sent assuming a t-long string of seeds
# and omega elements not to be sent.

# The worst case takes place whenever each additional seed not to be sent
# is a leaf of the highest remaining subtree. Once this element is added, the
# subtree is split into subtree-height smaller subtrees with all possible heights

def worst_case_seed_tree_cost(t,omega,sec_param):
    # If omega > t/2, every other leaf, and then more, must not be sent.
    # In the worst-case scenario, the seed-tree has the same cost as not using it
    if (omega > t//2):
        print("Warning, the constant weight challenge string is too dense, no gain.")
        return omega*sec_param

    tree_height = ceil(log2(t))
    available_subtrees_per_height = {tree_height: 1}
    num_seeds_to_send = 1
    #print(f"boostrap: {available_subtrees_per_height}, t:{t}, omega:{omega}")
    for i in range(omega):

        # removing a leaf to be sent adds the height of
        # the highest subtree to the seeds to be sent, unless only leaves are left
        #print(available_subtrees_per_height)
        height_highest = max(available_subtrees_per_height.keys())
        # remove one such subtree
        available_subtrees_per_height[height_highest] -= 1
        if (available_subtrees_per_height[height_highest] == 0):
            available_subtrees_per_height.pop(height_highest)
        # increment the number of seeds
        num_seeds_to_send += height_highest
        # add the subtrees created by the split: they have height in the 1
        # to height_highest-1 range
        for h in range(1,height_highest):
            if h in available_subtrees_per_height.keys():
                available_subtrees_per_height[h] += 1
            else:
                available_subtrees_per_height[h] = 1
    # each node of a tree is sec_param bits long, plus an index for the position
    # in the tree
    return num_seeds_to_send*(sec_param)


if __name__ == "__main__":
    if len(argv) < 2:
        print("Script to compute the maximum size (in bytes) of the reveals via seed-tree ")
        print(f"Usage: {argv[0]} num_rounds_t num_set_rounds_omega seed_size_bytes")
        exit(1)
    else:
        worst_case_size = worst_case_seed_tree_cost(int(argv[1]),
                                                    int(argv[2]),
                                                    int(argv[3])*8)
        print(f"worst case size {ceil(worst_case_size/8)}")
