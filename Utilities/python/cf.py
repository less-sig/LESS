#!/usr/bin/env python3
""" super simple canonical form implementation. The only goal is to have zero dependencies """

import random
import copy
import bitarray
from ctypes import c_ushort 
from Crypto.Cipher import AES
from Crypto.Hash import SHAKE256

from fq import Fq
from matrix import Matrix
from permutation import Permutation

def random_diagonal_matrix(k):
    """
    :param k: 
    :return: a matrix with random entries [1, q-1] on the main diagonal
    """
    Dr = Matrix(k, k, q)
    for i in range(k):
        Dr[i,i] = random.randint(1, q-1)
    return Dr


def lex_min_multisets(a_multiset, b_multiset):
    """
    Works only if q is a prime
    Returns the index corresponding to the minimum multiset
    Return -1 if the multisets are the same

    :param a_multiset
    :param b_multiset
    """
    for i in range(len(a_multiset)):
        if a_multiset[i] < b_multiset[i]:
            return 0
        else:
            if a_multiset[i] > b_multiset[i]:
                return 1
            
    #At this point, this means the two multisets are the same: return -1
    return -1


def lex_min_vectors(a_multiset, b_multiset):
    """
    Works only if q is a prime
    Returns the index corresponding to the minimum multiset
    If the vectors are the same, return the first one by default
    :param a_multiset
    :param b_multiset
    """ 
    ret = lex_min_multisets(a_multiset, b_multiset)        
    if ret == -1:
        return 0
    return ret


def lex_min_matrices(A, B):
    """
    -1 equal 
     0 smaller:A < B 
     1 bigger :A > B
    """
    assert A.nrows == B.nrows

    i = 0
    ret = -1
    while ret == -1 and i < A.nrows:
        ret = lex_min_multisets(A[i], B[i])
        i += 1

    return ret 


def sort_multisets(row_multisets):
    """
    Sort multisets using lexicograph ordering
    Raises an error if, at some point, the multisets are equal
    The function returns the permutation sorting rows

    :param row_multisets:
    :return :
    """
    # Representation of permutation as list of indices
    indices = [i for i in range(len(row_multisets))]
    
    swap = True
    # Loop until swap is true
    while swap: 
        
        swap = False
        for i in range(len(row_multisets)-1):
            
            #lex ordering of multisets
            min_index = lex_min_multisets(row_multisets[i], row_multisets[i+1]) 
            
            #Report failure if min_index == -1
            if min_index == -1:
                return -1
            
            #Swap elements, if needed
            if min_index == 1:
#               specify that a swap was done
                swap = True 
                    
                #Swap multisets
                tmp = row_multisets[i+1]
                row_multisets[i+1] = row_multisets[i]
                row_multisets[i] = tmp
                
                #Swap indices
                tmp = indices[i+1]
                indices[i+1] = indices[i]
                indices[i] = tmp
                
    return indices


def sort_vectors(vectors):
    """
    Sort vectors using lexicograph ordering
    The function returns the permutation sorting columns

    :param vectors:
    :return :
    """
    #Representation of permutation as list of indices
    indices = [i for i in range(len(vectors))]
    
    swap = True
    while swap: #Loop until swap is true
        
        swap = False
        for i in range(len(vectors)-1):
            min_index = lex_min_vectors(vectors[i], vectors[i+1]) #lex ordering of vectors
                        
            #Swap elements, if needed
            if min_index == 1:

                swap = True #specify that a swap was done
                    
                #Swap multisets
                tmp = vectors[i+1]
                vectors[i+1] = vectors[i]
                vectors[i] = tmp
                
                #Swap indices
                tmp = indices[i+1]
                indices[i+1] = indices[i]
                indices[i] = tmp
                
    return indices


def case_3_CF(B):
    """
    :param B: input matrix
    """
    n = B.nrows 
    m = B.ncols 
    row_multisets = []
    for i in range(B.nrows):
        b_i = copy.copy(B[i])
        b_i.sort()
        row_multisets.append(b_i)

    #Compute row permutation sorting rows 
    row_indices = sort_multisets(row_multisets)
    
    #Report failure if row permutation is not defined
    if row_indices == -1:
        return -1, -1, -1
    
    #Apply row permutation
    row_sorted_B = Matrix(n, m, q).zero()
    for i in range(n):
        for j in range(m):
            row_sorted_B[i,j] = B[row_indices[i], j]
    
    #Now, sort columns
    vectors = [row_sorted_B.col(i) for i in range(m)]
    
    #Compute permutation sorting columns
    col_indices = sort_vectors(vectors)
    
    #Apply column permutation
    CF_B = Matrix(n, m, q)
    for i in range(m):
        for j in range(n):
            CF_B[j, i] = row_sorted_B[j][col_indices[i]]

    return row_indices, col_indices, CF_B


def case_4_CF(B):
    """
    """
    n = B.nrows 
    m = B.ncols 
    Ap = Matrix(n, m, q).zero()

    for i in range(n):
        if all(B[i, 0] == B[i, j] for j in range(m)):
            for j in range(m): Ap[i, j] = 1
            continue

        v = B[i]
        s = sum(v) % q
        sp = sum([pow(t, (q-2), q) for t in v]) % q
        if s != 0:
            s = pow(s, -1, q)
        else:
            s = sp 
            if s == 0: return -1, -1, -1

        for j in range(m):
            Ap[i, j] = (s*v[j]) % q 

    return case_3_CF(Ap)



def case_5_CF(B):
    """
    """
    n = B.nrows 
    m = B.ncols 
    A_j = Matrix(n, m, q)
    for j in range(n):
        for i in range(m):
            A_j[i, j] = q -1

    for i in range(n):
        cont = False 
        for j in range(m):
            if B[i, j] == 0:
                cont = True
                break
        if cont:
            continue

        A = Matrix(n, m, q).zero()
        for j in range(n):
            sc = pow(B[i, j], -1, q)
            for k in range(m):
                A[k, j] = (sc * B[k, j]) % q
        
        t, _, T = case_4_CF(A)
        if t != -1:
            if lex_min_matrices(A_j, T):
                A_j = T
    
    return 0, 0, A_j


    

q = 127
k = 8
n = 2*k

A = Matrix(k, n-k, q).random()
for i in range(k):
    for j in range(k):
        A[i, j] = 120 - i - j

#_, _, B = case_3_CF(A)
#print("B=CF3(A)")
#print(B)

#_, _, B = case_4_CF(A)
#print("B=CF4(A)")
#print(B)

_, _, B = case_5_CF(A)
print("B=CF5(A)")
print(B)
exit(1)

Pr = Permutation(k).random().to_matrix(q)
Pc = Permutation(n-k).random().to_matrix(q)
Dr = random_diagonal_matrix(k)
Dc = random_diagonal_matrix(n-k)
A_prime = Pr*Dr*A*Dc*Pc

row, col, B = case_5_CF(A)
row_p, col_p, B_p = case_5_CF(A_prime)
print(B)
print(B_p)

# print(A)
# print(A_prime)
