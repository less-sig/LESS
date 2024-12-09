import numpy as np

def compute_single_multiset(x,k,q):
    '''
    Computes multiset of single vector
    '''
    multiset = np.zeros(q,dtype=int)
    for j in range(k):
        a_ij = x[j]
        multiset[a_ij] += 1

    return multiset

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

def sub_CF4(sub_A, z, k, q, Fq, min_multiset):
    '''
    z: number of rows in sub_A
    min_multiset: minimum multiset for current CF
    '''
    CF_exists = 1
    scaled_sub_A = np.zeros((z,k), dtype = int)
    for i in range(z):

        s = Fq(0)
        for j in range(k):
            s+= Fq(sub_A[i,j])

        if s!=Fq(0):
            s_inv = s**-1
            for j in range(k):
                scaled_sub_A[i,j] = s_inv*sub_A[i,j]
        else: #deal with s = 0
            s_prime = Fq(0)
            for j in range(k):
                s_prime += Fq(sub_A[i,j])**(q-2)

            if s_prime == Fq(0):
                CF_exists = 0
            else:
                for j in range(k):
                    scaled_sub_A[i,j] = s_prime*sub_A[i,j]

    #comptue multisets and see if one is less than minimum
#    print("scaled sub A:")
#    print(scaled_sub_A)
#    print(" ")

#    print("min_multiset = ",min_multiset)
#    print(" ")
    minimum_found = 0
    for i in range(z):
        multiset_i = compute_single_multiset(scaled_sub_A[i],k,q)
        index_min_multiset = lex_min_multisets(multiset_i, min_multiset, q)
        if index_min_multiset!=1:
            minimum_found = 1

    return CF_exists, minimum_found


def CF5_faster(A, k, q, Fq):

    #fill the matrix with maximum value in Fq
    #this can be skipped and replaced with the first CF we compute below
    B_max = np.zeros((k,k), dtype = int)
    for i in range(k):
        for j in range(k):
            B_max[i,j] = Fq(q-1)

    min_multiset = compute_single_multiset(B_max[0,:],k,q)

    #compute hamming weights
    zero_rows = []
    J = []
    max_num_zeros = 0
    for i in range(k):
        #count zeros of row i
        num_zeros = 0
        for j in range(k):
            if A[i,j] == Fq(0):
                num_zeros += 1
        if num_zeros > 0:
            zero_rows.append(i)

        if num_zeros > max_num_zeros:
            J = [i]
            max_num_zeros = num_zeros
        else:
            if num_zeros == max_num_zeros:
                J.append(i)

    sub_A = A[J,:]

    scaled_sub_A = np.zeros((len(J),k),dtype=int)
    scaled_A = np.zeros((k,k), dtype = int)


    CF_fail = 1 #this remains 1 if all matrices lead to a failure
    #consider all keys
    for i in range(k):
        if i not in zero_rows: #skip if i is in zero_rows (i-th row contains zeros for sure)
            min_multiset = compute_single_multiset(B_max[0,:],k,q)
            #we first scale the rows indexed by J
            coeffs = [A[i,j]**-1 for j in range(k)]

            #scale columns

            for j in range(k):
                for ell in range(len(J)):
                    scaled_sub_A[ell,j] = coeffs[j]*sub_A[ell,j]

            #now, apply row scaling due to CF4
            
            CF_exists, min_found = sub_CF4(scaled_sub_A, len(J), k, q, Fq, min_multiset)

            #continue only if min_found = 1
            if min_found:
                #scale full matrix A
                for ell in range(k):
                    for j in range(k):
                        scaled_A[ell,j] = coeffs[j]*A[ell,j]

                #call CF4
                B_i = CF4(scaled_A, k, q, Fq)   

                if (str(B_i) == str(-1)) == False:
                    CF_fail = 0
                    B_max = lex_min_matrices(B_max, B_i, k)

    if CF_fail:
        return -1
    else:
        return B_max