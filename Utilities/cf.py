#!/usr/bin/env python3
""" super simple canonical form implementation. The only goal is to have zero dependencies """

import copy
from Crypto.Cipher import AES
from Crypto.Hash import SHAKE256

from fq import Fq
from matrix import Matrix
from permutation import Permutation


q = 127
k = 3 
n = 7

A = Matrix(k, n-k, q).random()
Pr = Permutation(k).random().to_matrix(q)
Pc = Permutation(n-k).random().to_matrix(q)

#print(A)
#print(Pr)
#print(Pc)
A_prime = Pr*A*Pc
#print(A_prime)

#Works only if q is a prime
#Returns the index corresponding to the minimum multiset
#Return -1 if the multisets are the same
def lex_min_multisets(a_multiset, b_multiset):
    for i in range(len(a_multiset)):
        if a_multiset[i] < b_multiset[i]:
            return 0
        else:
            if a_multiset[i] > b_multiset[i]:
                return 1
            
    #At this point, this means the two multisets are the same: return -1
    return -1


#Works only if q is a prime
#Returns the index corresponding to the minimum multiset
#If the vectors are the same, return the first one by default
def lex_min_vectors(a_multiset, b_multiset):
    
    for i in range(len(a_multiset)):
        if a_multiset[i] < b_multiset[i]:
            return 0
        else:
            if a_multiset[i] > b_multiset[i]:
                return 1
            
    #At this point, this means the two multisets are the same: return 0
    return 0


#Sort multisets using lexicograph ordering
#Raises an error if, at some point, the multisets are equal
#The function returns the permutation sorting rows
def sort_multisets(row_multisets):
    
    #Representation of permutation as list of indices
    indices = [i for i in range(len(row_multisets))]
    
    swap = True
    while swap: #Loop until swap is true
        
        swap = False
        for i in range(len(row_multisets)-1):

            min_index = lex_min_multisets(row_multisets[i], row_multisets[i+1]) #lex ordering of multisets
            
            #Report failure if min_index == -1
            if min_index == -1:
                return -1
            
            #Swap elements, if needed
            if min_index == 1:

                swap = True #specify that a swap was done
                    
                #Swap multisets
                tmp = row_multisets[i+1]
                row_multisets[i+1] = row_multisets[i]
                row_multisets[i] = tmp
                
                #Swap indices
                tmp = indices[i+1]
                indices[i+1] = indices[i]
                indices[i] = tmp
                
    return indices


#Sort vectors using lexicograph ordering
#The function returns the permutation sorting columns
def sort_vectors(vectors):
    
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
    #print("B")
    #print(B)
    n = B.nrows 
    m = B.ncols 
    row_multisets = []
    for i in range(B.nrows):
        b_i = copy.copy(B[i])
        b_i.sort()
        row_multisets.append(b_i)

    #Compute row permutation sorting rows 
    #print(row_multisets)
    row_indices = sort_multisets(row_multisets)
    #print(row_indices)
    
    #Report failure if row permutation is not defined
    if row_indices == -1:
        return -1, -1, -1
    
    #Apply row permutation
    row_sorted_B = Matrix(n, m, q).zero()
    for i in range(n):
        for j in range(m):
            row_sorted_B[i,j] = B[row_indices[i], j]
    
    #print("row_sorted_B")
    #print(row_sorted_B)

    #Now, sort columns
    vectors = [row_sorted_B.col(i) for i in range(m)]
    
    #print("vectors")
    #print(vectors)
    #Compute permutation sorting columns
    col_indices = sort_vectors(vectors)
    #print("col_indices")
    #print(col_indices)
    
    #Apply column permutation
    CF_B = Matrix(n, m, q)
    for i in range(m):
        for j in range(n):
            CF_B[j, i] = row_sorted_B[j][col_indices[i]]

    return row_indices, col_indices, CF_B


row, col, B = case_3_CF(A)
row_p, col_p, B_p = case_3_CF(A_prime)
print(B)
print(B_p)
