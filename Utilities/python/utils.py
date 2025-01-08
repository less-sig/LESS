#!/usr/bin/env python3
"""
"""

import math


def yt_shuffle_state(rng, 
                     permutation: list[int],
                     N: int):
    """
    """
    POS_BITS = math.ceil(math.log2(N-1))
    POS_MASK = (1 << POS_BITS) - 1
    rand_u64 =int(rng.read(8).hex(), 16)
    c = 0
    x = 0
    for i in range(N):
        while 1:
            if c == ((64/POS_BITS)-1):
                c = 0
                rand_u64 =int(rng.read(8).hex(), 16)

            x = rand_u64 & POS_MASK
            rand_u64 = rand_u64 >> POS_BITS
            c = c + 1
            if x < N:
                break

        tmp = permutation[i]
        permutation[i] = permutation[x]
        permutation[x] = tmp

