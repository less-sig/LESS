import galois
import numpy as np
import random

from CF_faster_impl import CF5_faster

def sample_monomial(Fq, n):
    '''
    Returns the list of indices defining the permutation matrix of size n x n,
    and the list of scaling coefficients
    '''

    perm = np.array([i for i in range(n)])
    random.shuffle(perm)

    #sample scaling coefficients
    coeffs = []
    for i in range(n):
        a = Fq.Random(1)
        while a == 0:
            a = Fq.Random(1)
        coeffs.append(a)

    return [perm[i] for i in range(n)], coeffs

#########################################

def apply_monomial(Fq, n, k, G, perm, coeffs):
    '''
    Apply monomial matrix represented by perm and coeffs to input G,
    compute systematic form, return non systematic part
    '''

    #Apply monomial to G
    new_G = Fq.Zeros((k,n))
    for i in range(n):
        for j in range(k):
            new_G[j,i] = G[j, perm[i]]*coeffs[i]

#    print(new_G)
    #Compute systematic form
    new_A = SF(n, k, new_G)

    return new_A

#########################################

def SF(n, k, G):
    '''
    Brings matrix in SF form; for simplicity, we construct pivots only in the 
    first k columns. If they do not form a non singular matrix, the code doesn't work.
    The function returns only the non systematic part of the matrix
    '''
    sub_G = G[:,0:k]
#    print(sub_G)
    S = np.linalg.inv(sub_G)

    return S @ G[:,k:n]

##########################################

def KeyGen(Fq, n, k):
    '''
    Returns:
    - A: non systematic part of code C
    - new_A: non systematic part of code C'
    - perm, coeffs: secret key (monomial transformation)
    '''
    #Generate first code
    A = Fq.Random((k, n-k))

    G = build_full_generator_matrix(Fq, n, k, A)

    #Sample monomial transformation
    perm, coeffs = sample_monomial(Fq, n)

    #Get new code
    new_A = apply_monomial(Fq, n, k, G, perm, coeffs)

    return A, new_A, perm, coeffs

def build_full_generator_matrix(Fq, n, k, A):
    '''
    receives as input the non systematic part of a matrix, builds the full matrix
    '''

    G = np.zeros((k,n),dtype = int)
    for i in range(k):
        G[i,i] = Fq(1)
    for i in range(n-k):
        G[:, k+i] = A[:, i]

    return G

def inverse_perm(n, perm):
    '''
    takes as inpnut a list of indices representing a permutation,
    returns the list representing the inverse permutation
    '''

    inv_perm = []
    for i in range(n):
        j = perm.index(i)
        inv_perm.append(j)

    return inv_perm

####################################

def combine_perms(n, perm1, perm2):
    '''
    Combines perm1 with perm 2 (first applies perm1, then perm2)
    Returns a list representing the resulting permutation
    '''

    perm_result = []
    for i in range(n):
        j = perm1[perm2[i]]
        perm_result.append(j)
    return perm_result

################################

def compress(n, k, perm):
    '''
    Returns a vector with length n and ones only in the positions which are 
    employed for the information set (first k entries)
    '''
    rsp = np.zeros(n, dtype = int)
    for i in range(k):
        rsp[perm[i]] = 1
    
    return rsp

###############################

def verify_rsp(Fq, q, n, k, A_prime, rsp):
    '''
    Verifies response:
    - moves columns corresponding to ones in rsp to first k positions
    - applies SF
    - computes CF
    '''

    #Builds full matrix G_prime
    full_G_prime = build_full_generator_matrix(Fq, n, k, A_prime)

    #Applies permutation represented by rsp (moves only coordinates to left)
    perm_G_prime_left = Fq.Zeros((k,k))
    perm_G_prime_right = Fq.Zeros((k,n-k))
    
    num_left = 0; num_right = 0
    for i in range(n):
        if rsp[i] == 1:
            perm_G_prime_left[:,num_left] = full_G_prime[:,i]
            num_left += 1
        else:
            perm_G_prime_right[:,num_right] = full_G_prime[:,i]
            num_right += 1

    #computes non systematic part
    S = np.linalg.inv(perm_G_prime_left)
    B = S @ perm_G_prime_right

    #Computes CF
    verify_cmt = CF5_faster(B, k, q, Fq)

    return verify_cmt