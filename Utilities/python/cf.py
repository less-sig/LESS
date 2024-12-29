#!/usr/bin/env python3
""" super simple canonical form implementation. The only goal is to have zero dependencies """

import random
import copy
import bitarray
from ctypes import c_ushort

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
            
    # At this point, this means the two multisets are the same: return -1
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
            
            # lex ordering of multisets
            min_index = lex_min_multisets(row_multisets[i], row_multisets[i+1]) 
            
            # Report failure if min_index == -1
            if min_index == -1:
                return -1
            
            # Swap elements, if needed
            if min_index == 1:
                swap = True
                    
                # Swap multisets
                tmp = row_multisets[i+1]
                row_multisets[i+1] = row_multisets[i]
                row_multisets[i] = tmp
                
                # Swap indices
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
    # Representation of permutation as list of indices
    indices = [i for i in range(len(vectors))]
    
    swap = True
    while swap:
        swap = False
        for i in range(len(vectors)-1):
            # lex ordering of vectors
            min_index = lex_min_vectors(vectors[i], vectors[i+1])
                        
            # Swap elements, if needed
            if min_index == 1:
                # specify that a swap was done
                swap = True
                    
                # Swap multisets
                tmp = vectors[i+1]
                vectors[i+1] = vectors[i]
                vectors[i] = tmp
                
                # Swap indices
                tmp = indices[i+1]
                indices[i+1] = indices[i]
                indices[i] = tmp
                
    return indices


def compute_multisets(A, q):
    """
    NOTE: compute histograms
    Computes multiset for each row of input A
    (multisets are done as counters, to speed-up computations)
    """
    n = A.nrows
    m = A.ncols
    multisets = []
    for i in range(n):
        this_multiset = [0 for _ in range(m)]
        for j in range(m):
            a_ij = A[i, j]
            this_multiset[a_ij] += 1  # update counter

        multisets.append(this_multiset)  # update new multiset

    return multisets


def case_3_CF(B: Matrix):
    """
    :param B: input matrix
    """
    n = B.nrows 
    m = B.ncols 
    row_multisets = []
    for i in range(B.nrows):
        b_i = [B[i, j].get() for j in range(m)]
        b_i.sort()
        row_multisets.append(b_i)

    # Compute row permutation sorting rows
    row_indices = sort_multisets(row_multisets)
    
    # Report failure if row permutation is not defined
    if row_indices == -1:
        return -1, -1, -1
    
    # Apply row permutation
    row_sorted_B = Matrix(n, m, q).zero()
    for i in range(n):
        for j in range(m):
            row_sorted_B[i,j] = B[row_indices[i], j]
    
    # Now, sort columns
    vectors = [row_sorted_B.col(i) for i in range(m)]
    
    # Compute permutation sorting columns
    col_indices = sort_vectors(vectors)
    
    # Apply column permutation
    CF_B = Matrix(n, m, q)
    for i in range(m):
        for j in range(n):
            CF_B[j, i] = row_sorted_B[j][col_indices[i]]

    return row_indices, col_indices, CF_B


def CF3_tony(A: Matrix) -> (bool, Matrix):
    # get dimensions of input matrix A
    z = A.nrows
    k = A.ncols

    # compute multisets
    multisets = compute_multisets(A, q)

    # sort multisets
    sorting_indices = sort_multisets(multisets)
    # sorting_indices = bubble_sort_multisets(multisets, z, q)
    if sorting_indices == -1:  # report failure
        return False, A

    # Apply row permutation
    row_sorted_B = Matrix(z, k, q).zero()
    for i in range(z):
        for j in range(k):
            row_sorted_B[i,j] = A[sorting_indices[i], j]

    # Now, sort columns
    vectors = [row_sorted_B.col(i) for i in range(k)]

    # Compute permutation sorting columns
    col_indices = sort_vectors(vectors)

    # Apply column permutation
    CF_B = Matrix(z, k, q)
    for i in range(k):
        for j in range(z):
            CF_B[j, i] = row_sorted_B[j][col_indices[i]]

    return True, CF_B


def case_4_CF(B: Matrix):
    """
    standard version, taken from the paper
    """
    n = B.nrows 
    m = B.ncols 
    Ap = Matrix(n, m, q).zero()

    for i in range(n):
        if all(B[i, 0] == B[i, j] for j in range(m)):
            for j in range(m): Ap[i, j] = 1
            continue

        v = [t.get() for t in B[i]]
        assert type(v) == list
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


def case_4_CF_tony(B: Matrix,
                   nz:list ):
    """
    implements tonys improvement
    :param B: matrix
    :param nz: number of zeros in the matrix, in the current row
    """
    n = B.nrows 
    m = B.ncols 
    Ap = Matrix(n, m, q).zero()

    print(nz)

    for i in range(n):
        s = [0 for _ in range(n)]
        v = [t.get() for t in B[i]]
        touched = False
        for j in range(m):
            if v[j] != 0:
                touched = True
                s[nz[j]] += v[j]

        if not touched:
            return False
        
        print(s)
        t = 0
        for j in range(len(s)):
            if s[j] != 0:
                t = pow(s[j], (q-2), q)
                break

        for j in range(len(s)):
            Ap[i, j] = (t*v[j]) % q 

    return case_3_CF(Ap)


def sub_CF4(sub_A: Matrix, min_multiset):
    """
    faster CF4 for popcount cf5
    :param sub_A: sub matrix
    :param min_multiset: minimum multiset for current CF5
    """
    z = sub_A.nrows
    nc = sub_A.ncols
    exists = True
    min_found = False
    for i in range(z):
        # sum the current row
        s = sum([v.get() for v in sub_A[i]]) % q
        if s == 0:
            s = sum([pow(v.get(), (q - 2), q) for v in sub_A[i]]) % q
        else:
            s = pow(s, -1, q)

        if s == 0:
            exists = False
            continue

        # TODO:
        w = [(s * sub_A[i, j].get()) % q for j in range(nc)]
        w.sort()

        # compute the multiset and see if its less
        if lex_min_multisets(w, min_multiset) != 1:
            min_found = True

    return exists, min_found


def sub_CF4_tony(sub_A: Matrix) -> (bool, Matrix):
    """
    faster CF4 for popcount cf5
    :param sub_A: sub matrix
    :return
    """
    z = sub_A.nrows
    nc = sub_A.ncols

    scaled_sub_A = Matrix(z, nc, q)
    for i in range(z):
        # sum the current row
        s = sum([v.get() for v in sub_A[i]]) % q
        if s == 0:
            s = sum([pow(v.get(), (q-2), q) for v in sub_A[i]]) % q
        else:
            s = pow(s, -1, q)

        if s == 0:
            return False, scaled_sub_A

        s = Fq(s, q)
        for j in range(k):
           scaled_sub_A[i, j] = s*sub_A[i,j]

    return True, scaled_sub_A


def sub_CF4_tony_v2(sub_A: Matrix,
                 min_multiset: list,
                 nz: list):
    """
    NOTE: incompatible with the other CF4
    faster CF4 for popcount cf5 including tonys improvement
    :param sub_A: sub matrix 
    :param min_multiset: minimum multiset for current CF5
    """
    z = sub_A.nrows
    nc = sub_A.ncols
    exists = True
    min_found = False
    for i in range(z):
        s = [Fq(0, q) for _ in range(nc)]
        touched = False
        for j in range(nc):
            # can also remove this line
            if sub_A[i, j] != Fq(0, q): 
                touched = True
                s[nz[j]] += sub_A[i, j]
           
        if not touched:
            return False

        t = 0
        for j in range(len(s)):
            if s[j] != 0:
                # ret[r,:] = (v/s[j])[:]
                t = pow(s[j].get(), (q-2), q)

        w = [(t*sub_A[i, j].get())%q for j in range(nc)]
        w.sort()

        # compute the multiset and see if its less
        if lex_min_multisets(w, min_multiset) != 1:
            min_found = True
            return exists, min_found

    return exists, min_found


def case_5_CF(B: Matrix):
    """
    original version from the paper
    :param B:
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
            if B[i, j] == Fq(0, q):
                cont = True
                break
        if cont:
            continue

        A = Matrix(n, m, q).zero()
        for j in range(n):
            sc = pow(B[i, j].get(), -1, q)
            for k in range(m):
                A[k, j] = B[k, j] * sc
        
        t, _, T = case_4_CF(A)
        if t != -1 and lex_min_matrices(A_j, T):
            A_j = T
    
    return 0, 0, A_j


def case_5_CF_popcnt(B: Matrix):
    """
    version from paolo based on counting zeros
    """
    nr = B.nrows  
    nc = B.ncols
    # this remains 1 if all matrices lead to a failure
    CF_fail = 1 
    min_multiset = [q-1 for i in range(nc)]
    A_j = Matrix(nr, nc, q).set(q-1)

    # count zeros in each row
    max_zeros = 0
    row_has_zero = [False for _ in range(nr)]
    J = []
    for i in range(k):
        # count zeros of row i
        num_zeros = 0
        for j in range(k):
            if B[i,j] == Fq(0):
                num_zeros += 1
        if num_zeros > 0:
            row_has_zero[i] = True

        if num_zeros > max_zeros:
            J = [i]
            max_zeros = num_zeros
        else:
            if num_zeros == max_zeros:
                J.append(i)
  
    assert len(J) > 0
    sub_A = Matrix(len(J), nc, q)
    for i in range(len(J)):
        for j in range(nc):
            sub_A[i, j] = B[J[i], j]

    scaled_sub_A = Matrix(len(J), nc, q)
    scaled_A = Matrix(nr, nc, q)

    for i in range(k):
        # skip if i is in zero_rows (i-th row contains zeros for sure)
        if not row_has_zero[i]:
            min_multiset = [A_j[0, j].get() for j in range(nc)]
            min_multiset.sort()

            # we first scale the rows indexed by J
            coeffs = [pow(B[i, j].get(), (q-2), q) for j in range(nc)]
            
            # scale columns
            for j in range(k):
                for ell in range(len(J)):
                    scaled_sub_A[ell,j] = sub_A[ell,j]*coeffs[j]

            cf_exists, min_found = sub_CF4(scaled_sub_A, min_multiset)
            
            # continue only if min_found = 1
            if min_found:
                # scale full matrix A
                for ell in range(k):
                    for j in range(k):
                        scaled_A[ell,j] = B[ell,j]*coeffs[j]

                # call CF4
                t, _, B_i = case_4_CF(scaled_A)   
                if t != -1 and lex_min_matrices(A_j, B_i):
                    CF_fail = False
                    A_j = B_i

    if CF_fail:
        return -1, -1, A_j
    else:
        return 0, 0, A_j



def case_5_CF_tony(B: Matrix):
    """
    version from tony based on counting zeros
    """
    nr = B.nrows
    nc = B.ncols
    # this remains 1 if all matrices lead to a failure
    CF_fail = 1
    min_multiset = [q-1 for i in range(nc)]
    A_j = Matrix(nr, nc, q).set(q-1)

    # count zeros in each row
    max_zeros = 0
    non_zero_rows = [0 for _ in range(nr)]
    non_zero_rows_ctr = 0
    J = []
    for i in range(k):
        # count zeros of row i
        num_zeros = 0
        for j in range(k):
            if B[i,j] == Fq(0):
                num_zeros += 1
        if num_zeros == 0:
            non_zero_rows[non_zero_rows_ctr] = i
            non_zero_rows_ctr += 1

        if num_zeros > max_zeros:
            J = [i]
            max_zeros = num_zeros
        else:
            if num_zeros == max_zeros:
                J.append(i)

    assert len(J) > 0
    sub_A = Matrix(len(J), nc, q)
    for i in range(len(J)):
        for j in range(nc):
            sub_A[i, j] = B[J[i], j]

    scaled_sub_A = Matrix(len(J), nc, q)

    all_sub_CFs  = [] #collect all scaled versions of sub_A
    good_indices = [] #indices for which CF4 has not failed
    for i in non_zero_rows:
        # we first scale the rows indexed by J
        coeffs = [pow(B[i, j].get(), (q-2), q) for j in range(nc)]

        # scale columns
        for j in range(k):
            for ell in range(len(J)):
                scaled_sub_A[ell,j] = sub_A[ell,j]*coeffs[j]

        ok, this_scaled_sub_A = sub_CF4(scaled_sub_A, min_multiset)

        # continue only if min_found = 1
        if ok:
            ok, this_scaled_sub_A3 = CF3_tony(this_scaled_sub_A)
            if ok:
                all_sub_CFs.append(this_scaled_sub_A3)
                good_indices.append(i)


    # TODO unfinished
    return 0


q = 127
k = 8
n = 2*k

data = [
    [113, 99, 60, 37, 44, 36,  7,105], 
    [116, 90, 66, 37,  7, 43,111,111], 
    [  2, 12, 92, 96, 38, 41, 79, 49], 
    [109,119, 24, 36, 69,  0, 84, 99], 
    [ 36, 68, 64, 46, 27,124,107, 36], 
    [ 20,107, 63, 96,119, 83, 81, 27], 
    [ 59, 63, 91, 75, 37,  2, 33, 21], 
    [  0, 78, 89, 71,  4, 67, 12,  9]
]

correct = Matrix(k, n-k, q).set([
    [  0,   1,   7, 12, 34, 93,115,120],
    [ 39,  18,  40,  8,100,  0, 18, 32],
    [  1,   1,   1,  1,  1,  1,  1,  1],
    [ 98, 123, 107, 17, 52,  2,  3,107],
    [ 15,  87, 116, 49, 34,106, 82, 20],
    [ 75,  39,  22, 15, 23, 57, 89, 62],
    [114,  81, 106,103, 72,114, 28, 18],
    [ 75,  50, 112,117, 45,119, 97, 21],
])


def test_cf5():
    """ tests the first and simplest CF5 implementation """
    A = Matrix(k, n-k, q).random()
    A.set(data)
    _, _, B = case_5_CF(A)
    print("B=CF5(A)")
    print(B)
    print(B.eq(correct))


def test_cf5_paolo():
    A = Matrix(k, n-k, q).random()
    A.set(data)
    _, _, B = case_5_CF_popcnt(A)
    print("B=CF5_paolo(A)")
    print(B)
    print(B.eq(correct))


def test_cf5_tony():
    A = Matrix(k, n-k, q).random()
    A.set(data)
    _, _, B = case_5_CF_tony(A)
    print("B=CF5_tony(A)")
    print(B)
    print(B.eq(correct))


if __name__ == "__main__":
    # test_cf5()
    # test_cf5_paolo()
    test_cf5_tony()
    exit(1)

    A = Matrix(k, n-k, q).random()
    #_, _, B = case_3_CF(A)
    #print("B=CF3(A)")
    #print(B)
    
    #_, _, B = case_4_CF(A)
    #print("B=CF4(A)")
    #print(B)
    
    _, _, B = case_5_CF(A)
    print("B=CF5(A)")
    print(B)
    
    #_, _, B = case_5_CF_popcnt(A)
    #print("B=CF5_pop(A)")
    #print(B)
    #
    #_, _, B = case_5_CF_tony(A)
    #print("B=CF5_tony(A)")
    #print(B)
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
