import numpy as np
from code_utils import sample_monomial, KeyGen, build_full_generator_matrix, apply_monomial, inverse_perm, combine_perms, verify_rsp, compress
import random
import numpy as np
import galois

def compute_multisets(A,k,q):
    '''
    Computes multiset for each row of input A
    '''
    multisets = []
    for i in range(k):
        this_multiset = np.zeros(q,dtype=int)
        for j in range(k):
            a_ij = A[i,j]
            this_multiset[a_ij] += 1
        multisets.append(this_multiset)
    return multisets

def lex_min_multisets(multiset_1, multiset_2, q):
    '''
    Returns 0 if multiset_1 < multiset_2, 1 if multiset_1 > multiset_2, failure (-1) if they're equal 
    '''
    i = 0
    while (i<q):
        if multiset_1[i]>multiset_2[i]:
            return 0
        else:
            if multiset_1[i]<multiset_2[i]:
                return 1
        i+= 1
    return -1

def lex_min_vectors(vector_1, vector_2, k):
    '''
    Returns 0 if vector_1 < vector_2, 1 if vector_1 > vector_2, 2 if they're equal
    '''
    i = 0
    while (i<k):
        if vector_1[i]<vector_2[i]:
            return 0
        else:
            if vector_1[i]>vector_2[i]:
                return 1
        i+= 1
    return 2

def lex_min_matrices(A_1, A_2, k):
    '''
    Returns 0 if A_1 <= A_2, 1 if A_1 > A_2
    '''
    i = 0
    while (i<k):
        lex_result = lex_min_vectors(A_1[i,:], A_2[i,:], k)
        #if lex_result == 0, A_1 is minimum
        #if lex_result == 1, A_2 is minimum
        #if lex_result == 2, columns are equal, continue with next column
        if lex_result == 0:
            return A_1
        else:
            if lex_result == 1:
                return A_2
        i+=1

    return A_1


def bubble_sort_multisets(multisets, k, q):
    '''
    Sorts multisets, returns sorting permutation  (seen as a list of length k)
    '''
    sorting_indices = [i for i in range(k)]
    swap = True #swap = False when no swapping is done
    while swap:
        swap = False
        for i in range(k-1):
            compare_lex_i_i_plus_1 = lex_min_multisets(multisets[i], multisets[i+1], q)
            if compare_lex_i_i_plus_1 == -1: #multisets cannot be sorted, return failure
                return -1
            else:
                if compare_lex_i_i_plus_1 == 1: #do a swap
                    swap = True
                    tmp = multisets[i]
                    multisets[i] = multisets[i+1]
                    multisets[i+1] = tmp

                    #sort indices accordingly
                    tmp = sorting_indices[i]
                    sorting_indices[i] = sorting_indices[i+1]
                    sorting_indices[i+1] = tmp
    return sorting_indices

def bubble_sort_columns(X, k):
    '''
    Sorts columns so that they're in ascendind lex order
    '''
    sorting_indices = [i for i in range(k)]
    swap = True #swap = False when no swapping is done
    while swap:
        swap = False
        for i in range(k-1):
            compare_lex_i_i_plus_1 = lex_min_vectors(X[:,sorting_indices[i]], X[:, sorting_indices[i+1]], k)
            if compare_lex_i_i_plus_1 == 1: #do a swap
                swap = True
                tmp = sorting_indices[i]
                sorting_indices[i] = sorting_indices[i+1]
                sorting_indices[i+1] = tmp

    #apply sorting indices
    sorted_X = np.zeros((k,k), dtype = int)
    for i in range(k):
        sorted_X[:,i] = X[:, sorting_indices[i]]
    return sorted_X

def CF3(A,k,q):
    '''
    Computes CF (case 3) for input matrix A
    '''
    #compute multisets
    multisets = compute_multisets(A,k,q) 

    #sort multisets
    sorting_indices = bubble_sort_multisets(multisets, k, q) 
    if sorting_indices == -1: #report failure
        return -1
    
    #sort rows
    sorted_A = np.zeros((k,k), dtype=int)
    for i in range(k):
        sorted_A[i,:] = A[sorting_indices[i],:]

    #sort columns
    B = bubble_sort_columns(sorted_A, k)

    return B

def CF4(A, k, q, Fq):
    '''
    Computes CF (case 4) for input matrix A
    '''

    scaled_A = np.zeros((k,k), dtype = int)
    for i in range(k):
     #   print(i)
        #do stuff only if row i is non null
        #set row to (1,...,1) if it's all equal; return -1 if all values are 0
        all_same_vals = 1
        j = 0
        while (j<(k-1))&(all_same_vals):
            if A[i,j] != A[i,j+1]:
                all_same_vals = 0
            j+= 1
        if all_same_vals: #replace row with all ones
            if A[i,j] != Fq(0): #all values are equal and null
                for j in range(k):
                    scaled_A[i,j] = Fq(1)
        else: #values are different
            s = Fq(0)
            for j in range(k):
                s+= Fq(A[i,j])
        #    print("--> s = ",s)

            if s!=Fq(0):
                s_inv = s**-1
                for j in range(k):
                    scaled_A[i,j] = s_inv*A[i,j]
            else: #deal with s = 0
                s_prime = Fq(0)
                for j in range(k):
                    s_prime += Fq(A[i,j])**(q-2)
           #     print("--> s' = ",s_prime)

                if s_prime == Fq(0):
                    return -1
                else:
                    for j in range(k):
                        scaled_A[i,j] = s_prime*A[i,j]

    return CF3(scaled_A, k, q)

def CF5(A,k,q,Fq):

    #fill the matrix with maximum value in Fq
    B_max = np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            B_max[i,j] = Fq(q-1)

    CF_fail = 1 #this remains 1 if all matrices lead to a failure
    #consider all keys

    scaled_A = np.zeros((k,k),dtype=int)
    
    for i in range(k):
        
        #check if row i has 0; continue only if all values are non null
        j = 0
        zeros_exist = 0
        while (j<k)&(zeros_exist == 0):
            if A[i,j]==Fq(0):
                zeros_exist = 1
            j += 1

        if zeros_exist == 0: #all values are non null, scale and call CF4

            #scale columns
            for j in range(k):
                scalar_coeff_j = A[i,j]**-1
                for ell in range(k):
                    scaled_A[ell,j] = scalar_coeff_j*A[ell,j]

            #call CF4
            B_i = CF4(scaled_A, k, q, Fq)   

            if (str(B_i) == str(-1)) == False:
                CF_fail = 0
                B_max = lex_min_matrices(B_max, B_i, k)
    if CF_fail:
        return -1
    else:
        return B_max

if __name__ == "__main__":
    k = 8
    n = 2*8
    q = 127
    Fq = galois.GF(q); 

    A = Fq([
        [113, 99, 60, 37, 44, 36,  7,105], 
        [116, 90, 66, 37,  7, 43,111,111], 
        [  2, 12, 92, 96, 38, 41, 79, 49], 
        [109,119, 24, 36, 69,  0, 84, 99], 
        [ 36, 68, 64, 46, 27,124,107, 36], 
        [ 20,107, 63, 96,119, 83, 81, 27], 
        [ 59, 63, 91, 75, 37,  2, 33, 21], 
        [  0, 78, 89, 71,  4, 67, 12,  9]
    ])
    T = CF5(A, k, q, Fq)
    print(T)


