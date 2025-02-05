
import numpy as np
import galois
import time

from CF_basic_impl import CF5
from CF_faster_impl import CF5_faster
from CF_super_faster_impl import CF5_super_faster

def sample_perm(Fq, k):
    perm_list = np.random.permutation(k)
    P = Fq(np.zeros((k,k),dtype = int))
    
    for i in range(k):
        P[i, perm_list[i]] = Fq(1)
    return P

def sample_diag(Fq, q, k):

    D = Fq(np.zeros((k,k), dtype = int))
    for i in range(k):
        a = Fq.Random(1)
        while a == 0:
            a = Fq.Random(1)
        D[i,i] = Fq(a)
    return D

#test if CF is indeed an invariant

q = 127
k = 126
num_test = 1000


Fq = galois.GF(q)

pr_failure = 0
pr_canonical_ok = 0

time_faster = 0
time_super = 0
num_skipped = 0
size_J = 0

num_full_CF = 0

for id in range(num_test):

    #sample random matrix
    A = Fq.Random((k,k))

    #compute invariant matrix
    P_row = sample_perm(Fq, k)

    #print(P_row)
    P_col = sample_perm(Fq, k)
    #print(P_col)
    new_A = P_row @ A @ P_col
    #print(new_A)

    D_row = sample_diag(Fq, q, k)
    new_A = D_row @ new_A

    D_col = sample_diag(Fq, q, k)
    new_A = new_A @ D_col

    #measure time for standard
    start = time.time()
    B, num_skipped, size_J = CF5_faster(new_A,k,q,Fq,num_skipped, size_J)
    end = time.time()

    time_faster += (end-start)

    #measure time for faster CF5
    start = time.time()
    new_B, num_full_CF = CF5_super_faster(new_A,k,q,Fq, num_full_CF)
    end = time.time()

    time_super += (end-start)

    if str(B)==str(-1):
        pr_failure += 1
    if str(B)==str(new_B):
        pr_canonical_ok += 1
 #   print(B)
    print("Num tests = ",id+1,", Pr failure = ",pr_failure/(id+1),", Pr CF ok = ",pr_canonical_ok/(id+1))
#    print(B)
#    print(" ")
#    print(new_B)
#    print(str(B) == str(-1))
#    print(str(B) == str(new_B))
#    print("---")
    print("---> T(Faster) = ",time_faster/(id+1),", T(Super) = ",(time_super/(id+1)))
    print("---> Average Num full CF5 computations  = ",(num_full_CF/(id+1)))