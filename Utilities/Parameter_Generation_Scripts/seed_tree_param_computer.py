#!/usr/bin/env python3
from math import log2, ceil
from sys import argv,exit

def worst_case_seed_tree_cost(t,omega,sec_param):
    # If omega > t/2, every other leaf, and then more, must not be sent.
    # In the worst-case scenario, the seed-tree has the same cost as not using it
    if (omega > t//2):
        print("Warning, the constant weight challenge string is too dense, no gain.")
        return omega*sec_param
    
    return (sec_param/8)*(2**ceil(log2(omega))+omega*(ceil(log2(t))-ceil(log2(omega))-1))


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
