#CAREFUL: THE CODE IS SPECIFIC FOR RATE 1/2, i.e., k = n/2

def compute_single_multiset(x,k,q):
    '''
    Computes multiset of input vector x of length k
    (actually, this is done by populating a vector of counters so that element i is number of times element i appears in x
    '''
    multiset = np.zeros(q,dtype=int) #vector filled with zeros
    for j in range(k):
        x_j = x[j]
        multiset[x_j] += 1

    return multiset

def compute_multisets(A,q):
    '''
    Computes multiset for each row of input A
    (multisets are done as counters, to speed-up computations)
    '''
    multisets = [] #list containing multisets
    z = A.shape[0] #number of rows of input A
    k = A.shape[1] #number of columns of input A
    
    for i in range(z):
        this_multiset = np.zeros(q,dtype=int)
        for j in range(k):
            a_ij = A[i,j]
            this_multiset[a_ij] += 1 #update counter

        multisets.append(this_multiset) #update new multiset

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

def lex_min_vectors(vector_1, vector_2, m):
    '''
    Returns 0 if vector_1 < vector_2, 1 if vector_1 > vector_2, 2 if they're equal
    '''
    i = 0
    while (i<m):
        if vector_1[i]<vector_2[i]:
            return 0
        else:
            if vector_1[i]>vector_2[i]:
                return 1
        i+= 1
    return 2

def lex_min_matrices(A_1, A_2):
    '''
    Returns:
    - 0 if A_1 < A_2 
    - 1 if A_1 > A_2 
    - 2 if A_1 = A_2
    Matrices are compared considering rows (e..g, A_1 is less than A_2 if first row of A_1 is less than first row of A_2)
    '''
    z = A_1.shape[0] #number of rows of A_1 and A_2
    k = A_1.shape[1] #number of columns of A_1 and A_2
    i = 0
    while (i<z):
        lex_result = lex_min_vectors(A_1[i,:], A_2[i,:], k) #compare rows 
        #if lex_result == 0, A_1 is minimum
        #if lex_result == 1, A_2 is minimum
        #if lex_result == 2, columns are equal, continue with next column
        if lex_result == 0:
            return 0
        else:
            if lex_result == 1:
                return 1
        i+=1

    return 2


def bubble_sort_multisets(multisets, z, q):
    '''
    Sorts multisets, returns sorting permutation  (seen as a list of length k)
    - z: number of multisets
    - q: finite field size
    '''
    sorting_indices = [i for i in range(z)]
    swap = True #swap = False when no swapping is done
    while swap:
        swap = False
        for i in range(z-1):
            compare_lex_i_i_plus_1 = lex_min_multisets(multisets[i], multisets[i+1], q)
            if compare_lex_i_i_plus_1 == -1: #multisets cannot be sorted, return failure
                return -1
            else:
                if compare_lex_i_i_plus_1 == 1: #do a swap

                    #Swap multisets
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
    Sorts columns so that they're in ascending lex order
    Input X must have k columns
    '''
    sorting_indices = [i for i in range(k)]
    swap = True #swap = False when no swapping is done
    z = X.shape[0]
    
    while swap:
        swap = False
        for i in range(k-1):
            compare_lex_i_i_plus_1 = lex_min_vectors(X[:,sorting_indices[i]], X[:, sorting_indices[i+1]], z)
            if compare_lex_i_i_plus_1 == 1: #do a swap
                swap = True
                tmp = sorting_indices[i]
                sorting_indices[i] = sorting_indices[i+1]
                sorting_indices[i+1] = tmp

    #apply sorting indices
    sorted_X = np.zeros((z, k), dtype = int)
    for i in range(k):
        sorted_X[:,i] = X[:, sorting_indices[i]]

    return sorted_X


def CF3(A, q):
    '''
    Computes CF (case 3) for input matrix A
    '''

    #get dimensions of input matrix A
    z = A.shape[0] #number of rows of A matrix
    k = A.shape[1] #number of rows of A matrix

    #compute multisets
    multisets = compute_multisets(A,q) 

    #sort multisets
    sorting_indices = bubble_sort_multisets(multisets, z, q) 
    if sorting_indices == -1: #report failure
        return -1
    
    #sort rows
    sorted_A = np.zeros((z,k), dtype=int)
    for i in range(z):
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

    return CF3(scaled_A, q)

def sub_CF4(sub_A, z, k, q, Fq):
    '''
    Scales rows of sub_A; returns failure if a row cannot be scaled
    Input parameters:
        - sub_A: matrix to be scaled
        - z: number of rows in sub_A
        - k: number of columns of sub_A
        - q: finite field size
        - Fq: finite field
    '''
    scaled_sub_A = np.zeros((z,k), dtype = int)
    for i in range(z):

        #s = sum of elements in i-th row of sub_A
        s = Fq(0)
        for j in range(k):
            s+= Fq(sub_A[i,j])

        #if s != 0, scale row i of sub_A
        if s!=Fq(0):
            s_inv = s**-1
            for j in range(k):
                scaled_sub_A[i,j] = s_inv*sub_A[i,j]

        else: #deal with s = 0: compute s'
            s_prime = Fq(0)
            for j in range(k):
                s_prime += Fq(sub_A[i,j])**(q-2)

            if s_prime == Fq(0):
                return -1 #returns -1 as first argument cause CF cannot be computed
            else:
                for j in range(k):
                    scaled_sub_A[i,j] = s_prime*sub_A[i,j]

    return scaled_sub_A


def sort_sub_CFs(all_CF_matrices, indices, q):
    '''
    Input:
    - all_matrices: list with matrices [A_i1, A_i2, ...]
    - indices: corresponding indices [i1, i2, ...]
    The function:
    - for each A_i1, computes B_i1 = CF3(A_i1)
    - sorts the B_i1 so that they are in ascending order
    Returns the sorting indices, together with the sorted B matrices
    '''

    #sort all the CFs
    swap = True
    while swap == True:
        swap = False
        for i in range(len(indices)-1):
            compare_lex_i_i_plus_1 = lex_min_matrices(all_CF_matrices[i], all_CF_matrices[i+1])
            if compare_lex_i_i_plus_1 == 1: #do a swap
                
                swap = True
                
                #swap matrices
                tmp = all_CF_matrices[i]
                all_CF_matrices[i] = all_CF_matrices[i+1]
                all_CF_matrices[i+1] = tmp
                
                #swap indices
                tmp = indices[i]
                indices[i] = indices[i+1]
                indices[i+1] = tmp
                
    return indices, all_CF_matrices

def CF5_single_row(A, i, q, Fq):
    '''
    Computes CF5 using the elements in row i
    '''
    
    k = A.shape[0]
        
    #we know for sure the row does not have zeros, we proceed to scale columns
    scaled_A = np.zeros((k,k),dtype=int)
    
    #scale columns
    for j in range(k):
        scalar_coeff_j = A[i,j]**-1
        for ell in range(k):
            scaled_A[ell,j] = scalar_coeff_j*A[ell,j]

    #call CF4
    print("----> calling CF4")
    B_i = CF4(scaled_A, k, q, Fq)   

    return B_i


def CF5_super_faster(A, k, q, Fq, num_full_CF):
    '''
    Paolo's + Tony's ideas
    '''

    #Find rows having largest number of zeros, save row indices in set J
    #Save also rows not having zeros, cause they will be skipped later on when 
    #scaling columns
    non_zero_rows = [] #set of rows having zeros
    J = [] #set with rows having largest number of zeros
    max_num_zeros = 0 #largest number of zeros in a row
    for i in range(k):

        #count zeros of row i
        num_zeros = 0
        for j in range(k):
            if A[i,j] == Fq(0):
                num_zeros += 1

        #if row doesn't have zeros, append index to non_zero_rows
        if num_zeros == 0:
            non_zero_rows.append(i)

        #Update max number of zeros, or set J
        if num_zeros > max_num_zeros:
            J = [i]
            max_num_zeros = num_zeros
        else:
            if num_zeros == max_num_zeros:
                J.append(i)

    #Extract rows indexed by J
    sub_A = A[J,:]

    #Prepare matrix for scaling of sub_A
    scaled_sub_A = np.zeros((len(J),k),dtype=int)

    CF_fail = 1 #this remains 1 if all matrices lead to a failure

    #consider all rows that do not have zeros
    #apply column and row scaling only to sub_A
    all_sub_CFs  = [] #collect all scaled versions of sub_A
    good_indices = [] #indices for which CF4 has not failed

    for i in non_zero_rows:

        #we first scale the rows indexed by J
        coeffs = [A[i,j]**-1 for j in range(k)]

        #scale columns of sub_A
        for j in range(k):
            for ell in range(len(J)):
                scaled_sub_A[ell,j] = coeffs[j]*sub_A[ell,j]
                
#        print(">>>>>>>>>>>>>>>", scaled_sub_A)

        #now, apply row scaling due to CF4
        this_sub_CF_4 = sub_CF4(scaled_sub_A, len(J), k, q, Fq)
#        print(">>>>>>>>>>>>>>>>>>>>>> sub_CF4 = ", this_sub_CF_4)
                    
        if str(this_sub_CF_4)!=str(-1):
            this_sub_CF_3 = CF3(this_sub_CF_4, q)
#            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CF3 = ", this_sub_CF_3)
            if str(this_sub_CF_3) !=str(-1):
                all_sub_CFs.append(this_sub_CF_3)
                good_indices.append(i)
    
#    print("All sub-CFs computed")
    #Now, sort the matrices in all_scaled_sub_A
    #FOR FLOYD: this is done with bubble sort, very shitty :)
    
    sorting_indices, sorted_sub_CFS = sort_sub_CFs(all_sub_CFs, good_indices, q)
    
    #compute CFs until a well defined CF is found
    #FOR FLOYD: full CFs are computed only here
    flag_CF_ok = False
    i = 0
#    print("---> ",sorted_sub_CFS[0])
#    print("Final computation of CF")
    while (flag_CF_ok == False)&(i<len(sorting_indices)):
 #       print("computing CF")
        B_i = CF5_single_row(A, sorting_indices[i], q, Fq)
        num_full_CF += 1
 #       print(" CF computed")        
        if str(B_i)!=str(-1):
            flag_CF_ok = True        
        else:
            i+= 1
    
    if i == len(sorting_indices):
        return -1, num_full_CF
    
    #deal with the fact that more CFs may need to be computed
    while (str(sorted_sub_CFS[i]) == str(sorted_sub_CFS[i+1]))&(i<(len(sorting_indices)-1)):
        new_B_i = CF5_single_row(A, sorting_indices[i+1], q, Fq)
        num_full_CF += 1

        #compare this with old CF, update it if necessary
        min_index = lex_min_matrices(B_i, new_B_i)
        if min_index == 1:
            B_i = new_B_i

    return B_i, num_full_CF


if __name__ == "__main__":
    from code_utils import sample_monomial, KeyGen, build_full_generator_matrix, apply_monomial, inverse_perm, combine_perms, verify_rsp, compress
    import random
    import galois
    # test_full_implementation()
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
    T = CF5_super_faster(A, k, q, Fq)
    print(T)
