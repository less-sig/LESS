#!/usr/bin/env python3
import random

from fq import Fq
from matrix import Matrix

class Monomial:
    def __init__(self, n: int, q: int = 127) -> None:
        """ container class for a monomial """
        self.n = n
        self.q = q

        self.values = [Fq(0, self.q) for _ in range(self.n)]
        self.perm = [i for i in range(self.n)]

    def random(self, seed=None) -> 'Monomial':
        """ samples a random monomial 
        """
        if seed is None:
            # choose a random permutation:
            for i in range(self.n):
                pos = random.randint(0, self.n - i - 1)
                tmp = self.perm[pos]
                self.perm[pos] = self.perm[i]
                self.perm[i] = tmp

            # choose random values
            for i in range(self.n):
                self.values[i] = Fq(0, self.q).random(1)
        else:
            assert(False) # not implementated
        return self

    def inv(self) -> 'Monomial':
        """ TODO not finished mod fq is missing """
        ret = Monomial(self.n, self.q)
        for i in range(self.n):
            k = 0
            for j in range(self.n):
                if self.perm[j] == i:
                    k = j
                    break
            ret.perm[i] = self.perm[k]
        return ret

    def comb(self, A: 'Monomial') -> 'Monomial':
        """ """
        assert self.q == A.q
        assert self.n == A.n
        ret = Monomial(self.n, self.q)
        for i in range(self.n):
            j = self.perm[A.perm[i]]
            ret.perm[i] = j
        return ret
    
    def apply(self, G: Matrix) -> Matrix:
        """Apply monomial matrix represented by perm and coeffs to input G
        NOTE: not inplace
        """
        n = G.nrows
        m = G.ncols
        assert m == self.n
        A = Matrix(n, m, self.q)
        for i in range(m):
            for j in range(n):
                A[j, i] = G[j, self.perm[i]] * self.values[i]
        return A

    def __str__(self) -> str:
        M = Matrix(self.n, self.n, self.q).zero()
        for i in range(self.n):
            M[self.perm[i], i] = self.values[i]

        return str(M)


# just some simple tests
if __name__ == "__main__":
    n, r, q = 10, 5, 2
    m1 = Monomial(n, q).random()
    print(m1)

    G1 = Matrix(r, n, q).id()
    print(G1)
    G2 = m1.apply(G1)
    print(G2)


