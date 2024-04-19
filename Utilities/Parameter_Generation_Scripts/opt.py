#!/usr/bin/env python3
import math
from CryptographicEstimators.cryptographic_estimators.LEEstimator import *


def clog2(x: float):
    """ computes ceil(log(x)) """
    return math.ceil(math.log2(x))


# some hardcoded constants
q = 127
clq = clog2(q - 1)
s = 2

l_tree_seed = { 1: 16, 3: 24, 5: 32 }
l_sec_seed = { 1: 32, 3: 48, 5: 64 }
l_pub_seed = { 1: 16, 3: 24, 5: 32 }
l_salt = { 1: 32, 3: 48, 5: 64 }
l_digest = { 1: 32, 3: 48, 5: 64 }
min_bit_complexity = { 1: 143, 3: 207, 5: 272}

def isobits(case: int, n: int, k: int) -> int:
    """ computes the needed bits for each canonical form """
    assert 1 < case < 6
    if case == 2:
        return k * (clog2(n) + clq)
    if case == 3:
        return n
    if case == 4:
        return n + (k * clq)
    if case == 5:
        return n 
    
    # should never be called
    assert False

def Lambda(t: int, w: int) -> float :
    """ return the number of seeds in the tree """
    return math.pow(2, clog2(w)) + w * (clog2(t) - clog2(w) - 1.)


def signature_size(cat: int, case: int, n: int, k: int, t: int, w: int):
    """
    :param int n: code length
    :param int k: code dimension
    :param int s: number of executions in parallel: number of G_i in the pk
    :param int s: number of protocol repetitions
    :param int w: number of non null challenges
    """
    return w * math.ceil(isobits(case, n, k) / 8) + Lambda(t, w)*l_tree_seed[cat] + \
            l_salt[cat] + l_digest[cat]

min_size = math.inf
bit_complexity = math.inf
min_params = {"t": 0, "w": 0}
cat = 1
for n in range(252, 253):
    k = n // 2
    for case in range(3, 4):
        for t in range(200, 241):
            for w in range(10, 19):
                time = LEEstimator(n, k, q).fastest_algorithm().time_complexity()
                time += math.log2(n)
                if time < min_bit_complexity[1]:
                    continue

                sig_size = signature_size(cat, case, n, k, t, w)
                if sig_size < min_size:
                    min_size = sig_size
                    bit_complexity = time
                    min_params["t"] = t
                    min_params["w"] = w


print(min_size, min_params, bit_complexity)
# print(signature_size(1, 3, 252, 126, 244, 20))
# print(signature_size(1, 2, 252, 126, 244, 20))
