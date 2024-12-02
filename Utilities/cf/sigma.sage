import sys

q = 11
F = GF(q)
k = 10
n = k*2 

assert (n-k) % F.characteristic() != 0

op_row_scaling = 1 # choose one of the two strategies for row scaling
assert op_row_scaling in [0,1]

######################################

# turning input vector into multiset
def into_multiset(v):

  ret = list(v)
  ret.sort()

  return vector(ret)

# scaling rows of input matrix
def rows_scaling_old(M):

  ret = M[:,:]

  for r in range(M.nrows()):

    v = vector(M[r,:])

    if not v.is_zero():

      s0 = sum([ a for a in v ])
      s1 = sum([ (a^(-1) if (a != 0) else 0) for a in v])
 
      if s0 != 0:
        ret[r,:] = (v/s0)[:]
 
      elif s1 != 0:
        ret[r,:] = (v*s1)[:]
 
      else:
        return False

  return ret

# scaling rows of input matrix
def rows_scaling_new(M, nz):

  ret = M[:,:]

  for r in range(M.nrows()):

    v = vector(M[r,:])

    if not v.is_zero():

       s = [ F(0) ] * (M.nrows()+1) # can also reduce the length
       for j in range(M.ncols()):
         if v[j] != 0: # can also remove this line
           s[ nz[j] ] += v[j]

       if s.count(0) == len(s):
         return False

       for j in range(len(s)):
         if s[j] != 0:
           ret[r,:] = (v/s[j])[:]
           break
       
  return ret

# scaling cols of input matrix so that row r becomes (1,1,...,1)
def cols_scaling(M, r):

  ret = M[:,:]
  for i in range(M.ncols()):
    assert M[r,i] != 0
    ret[:,i] = ret[:,i] / M[r,i]

  return ret

def key_vec(pair):

  return pair[1]

# sorting rows based on their multisets and then sorting columns
def rows_cols_sorting(M):

  # sorting rows

  L = [ (v, into_multiset(v)) for v in M.rows() ]
  L_sorted = sorted(L, key=key_vec)

  # checking if there are different rows with the same multisets
  for i in range(len(L_sorted)-1): 
    if (L_sorted[i][0] != L_sorted[i+1][0]):
      if (L_sorted[i][1] == L_sorted[i+1][1]):
        return False

  # sorting cols

  tmp = matrix([e[0] for e in L_sorted])

  L = [ v for v in tmp.columns() ]
  return matrix([ e for e in sorted(L) ]).transpose()


# canonical form computation
def canonical_form(M):

  # finding rows without 0

  L_r = [] 
  for i in range(M.nrows()):
    if F(0) not in M.row(i):
      L_r.append(i)

  cf = False

  if op_row_scaling == 1:
    nz = [ list(M[:,j]).count(0) for j in range(M.ncols()) ]

  for i in L_r:

    # scaling rows
 
    Ms = cols_scaling(M, i) # no failure

    # scaling cols

    if op_row_scaling == 0:
      result = rows_scaling_old(Ms)
    else:
      result = rows_scaling_new(Ms, nz)

    if result == False:
      continue

    Ms[:,:] =  result

    # sorting rows and cols

    result = rows_cols_sorting(Ms)
    if result == False:
      continue

    # compare with cf candidate

    cand = copy(result)

    if cf == False:
      cf = cand[:,:]

    else:

      if vector(cand) < vector(cf):
        cf[:,:] = cand[:,:]

  #

  return cf

######################################

def random_units(n):

  return [ F.unit_group().random_element() for i in range(n) ]

def random_permutation(n):

  return Permutations(list(range(n))).random_element()

def scale_cols(M, L):

  ret = copy(M)

  for i in range(len(L)):
    ret[:,i] = (ret.column(i) * L[i])[:]

  return ret

def permute_cols(M, p):

  ret = matrix(F, M.nrows(), M.ncols())
  for i in range(len(p)):
    ret[:,i] = M[:,p[i]]

  return ret

def random_mu(n, w):

  ret = vector(F, n)
  while len(ret.nonzero_positions()) < w:
    ret[ randint(0,n-1) ]= F(1)

  return vector(ret)

def apply_mu(M, mu):

  ret = copy(M)

  i0 = 0
  i1 = M.nrows()

  for i in range(len(mu)):

    if mu[i] != F(0): 
      ret[:,i0] = M[:,i]
      i0 += 1  

    else:
      ret[:,i1] = M[:,i]
      i1 += 1  

  return ret

def find_pivot_cols(M):

  ret = [ F(0) ] * M.ncols()

  i = 0
  for j in range(M.ncols()):
    if M[i,j] == F(1) and vector(M[:,j]).hamming_weight() == 1:
      ret[j] = F(1)
      i += 1
      if i >= M.nrows():
        break

  return ret

######################################

I_k = identity_matrix(F, k)
iters = 10^3
fail = 0

for i in range(iters):

  # key gen

  while True:
  
    G0 = random_matrix(F, k, n)
    G0 = G0.echelon_form()
  
    if G0.rank() == k:
      break
  
  G1 = copy(G0)
  L_s = random_units(n)
  G1 = scale_cols(G1, L_s)
  p = random_permutation(n)
  G1 = permute_cols(G1, p)
  G1 = G1.echelon_form()
  
  # one round

  for t in range(iters):
  
    print ("pair {0}, round {1}, fail {2}".format(i, t, fail))
    sys.stdout.flush()
   
    #################################################################
  
    p0 = random_permutation(n) # can be generated from a seed; for challenge=0
    A0 = permute_cols(G0, p0).echelon_form() # permute columns

    pivots = find_pivot_cols(A0)
    A0 = apply_mu(A0, pivots).echelon_form() # permute columns to turn the matrix into [I | *]
  
    mu1 = [ pivots[ p0.index( p[i] ) ] for i in range(n) ] # mu1 indicates which cols of G1 should be the first k ones; for challenge=1

    A1 = apply_mu(G1, mu1).echelon_form()

    assert A0.submatrix(0,0,k,k) == I_k
    assert A1.submatrix(0,0,k,k) == I_k
  
    #################################################################
  
    M0 = A0.submatrix(0,k,k,n-k)
    M1 = A1.submatrix(0,k,k,n-k)
  
    M0_c = canonical_form(M0)
    M1_c = canonical_form(M1)
  
    assert (M0 != M1) # sanity check
    assert (M0_c == M1_c) # results must be the same (can be both False)
 
    if M0_c == False:
      fail += 1

