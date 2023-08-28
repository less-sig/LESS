#!/usr/bin/python3
from math import log2, ceil
from sys import argv,exit

# This function computes the worst-case seed tree cost by operatively
# determining the number of seeds to be sent assuming a t-long string of seeds
# and omega elements not to be sent.

def worst_case_seed_tree_cost(t,omega,sec_param):
    # If omega > t/2, every other leaf, and then more, must not be sent.
    # In the worst-case scenario, the seed-tree has the same cost as not using it
    if (omega > t//2):
        print("Warning, the constant weight challenge string is too dense, no gain.")
        return omega*sec_param
    
    return (sec_param/8)*ceil(omega*log2(t/omega))


if __name__ == "__main__":
    if len(argv) < 2:
        print("Script to compute the maximum size (in bytes) of the reveals via seed-tree ")
        print(f"Usage: {argv[0]} num_rounds_t num_set_rounds_omega sec_param_lambda")
        exit(1)
    else:
        worst_case_size = worst_case_seed_tree_cost(int(argv[1]),
                                                    int(argv[2]),
                                                    int(argv[3]))
        print(f"worst case size {ceil(worst_case_size)}")
