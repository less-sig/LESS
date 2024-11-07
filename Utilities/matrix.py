#!/usr/bin/env python3
""" super simple matrix implementation. The only goal is to have zero dependencies """

from typing import Union
import random

from fq import Fq


class Matrix:
    """ simple matrix class """
    def __init__(self, nrows: int, ncols: int, q: int = 2) -> None:
        """ zero initialized """
        self.nrows = nrows
        self.ncols = ncols
        self.q = q
        self.data = [[0 for _ in range(ncols)] for _ in range(nrows)] 

    def __getitem__(self, tup):
        """ nice access function
        NOTE: access it via:
            A[i, j] or  (returns field element)
            A[i]        (returns row)
        """
        if isinstance(tup, tuple):
            x, y = tup
            assert x < self.nrows and y < self.ncols
            return self.data[x][y]
        
        assert isinstance(tup, int)
        assert tup < self.nrows
        return self.data[tup]

    def __setitem__(self, tup, data):
        if isinstance(tup, tuple):
            x, y = tup
            assert x < self.nrows and y < self.ncols
            assert isinstance(data, int)
            self.data[x][y] = data % self.q
            return
        
        assert isinstance(tup, int)
        assert isinstance(data, list)
        assert tup < self.nrows
        self.data[tup] = data
        return

    def print(self, transpose: bool = False):
        """ printing """
        if transpose:
            for i in range(self.ncols):
                for j in range(self.nrows):
                    print(self.data[j][i], end='')
                print("")
            return
    
        for i in range(self.nrows):
            for j in range(self.ncols):
                print(self.data[i][j], end='')

            print("")
    
    def id(self) -> 'Matrix':
        """ identity matrix """
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.data[i][j] = i == j
        return self

    def zero(self) -> 'Matrix':
        """ zeros all elements"""
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.data[i][j] = 0
        return self

    def row(self, i: int):
        """ returns i-th row """ 
        assert i < self.nrows
        return self.data[i]

    def col(self, i: int):
        """ returns i-th col """ 
        assert i < self.ncols
        return [self.data[j][i] for j in range(self.nrows)]


    def random(self) -> 'Matrix':
        """ generates a random matrix """
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.data[i][j] = random.randint(0, self.q - 1)
        return self

    def random_row_with_weight(self, row: int, w: int) -> 'Matrix':
        """ generates a random weight w row """
        assert w > 0 and w < self.ncols
        self.zero()
        for i in range(w):
            self.data[row][i] = random.randint(1, self.q)

        # and now just simple apply a random permutation
        for i in range(self.ncols):
            pos = random.randint(0, self.ncols - i - 1)
            tmp = self.data[row][i]
            self.data[row][i] = self.data[row][i + pos]
            self.data[row][i + pos] = tmp
        return self

    def gauß(self, max_rank: Union[int, None] = None) -> int:
        """ simple Gaussian elimination. Is an inplace operation
        :return the rank of the matrix
        """
        if max_rank is None:
            max_rank = self.nrows
        
        assert isinstance(max_rank, int)
        row = 0
        for col in range(self.ncols):
            if row >= min(max_rank, self.nrows): break

            # find pivot
            sel = -1
            for i in range(row, self.nrows):
                if self.data[i][col] == 1:
                    sel = i 
                    break

            if sel == -1:
                return row

            self.__swap_rows(sel, row)

            # solve remaining coordinates
            for i in range(self.nrows):
                if i == row: continue 
                if self.data[i][col] == 0: continue

                for j in range(self.ncols):
                    self.data[i][j] += self.data[row][j]
                    self.data[i][j] %= self.q

            row += 1
        
        return row

    def mul(self, B: 'Matrix') -> 'Matrix':
        """ simple multiplication """
        B_r, B_c = B.nrows, B.ncols
        assert self.q == B.q and self.ncols == B_r
        C = Matrix(self.nrows, B_c, self.q)
        # each column in B
        for i in range(B_c):
            # each row in A
            for j in range(self.nrows):  
                sum = 0
                # each element in a row in A
                for k in range(self.ncols):  
                    sum += self[j, k] * B[k, i]

                C.data[j][i] = sum % self.q
        return C

    def add(self, B: 'Matrix') -> 'Matrix':
        """ simple inplace additions """
        B_r, B_c = B.nrows, B.ncols
        assert self.q == B.q and self.ncols == B_c and self.nrows == B_r
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.data[i][j] += B[i, j]
                self.data[i][j] %= self.q
        return self

    def sub(self, B: 'Matrix') -> 'Matrix':
        """ simple inplace subtraction """
        B_r, B_c = B.nrows, B.ncols
        assert self.q == B.q and self.ncols == B_c and self.nrows == B_r
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.data[i][j] = self.data[i][j] + (self.q - B[i, j])
                self.data[i][j] %= self.q
        return self

    def eq(self, B: 'Matrix'):
        """ simple comparsion """
        B_r, B_c = B.nrows, B.ncols
        assert self.q == B.q and self.ncols == B_c and self.nrows == B_r
        for i in range(self.nrows):
            for j in range(self.ncols):
                if self.data[i][j] != B[i, j]:
                    return False
        return True

    def transpose(self) -> 'Matrix':
        """ simple transpose, not inplace """
        T = Matrix(self.ncols, self.nrows, q=self.q)
        
        for i in range(self.nrows):
            for j in range(self.ncols):
                T.data[j][i] = self.data[i][j]
        return T

    def popcnt_row(self, row: int) -> int:
        """ computes the hamming weight of a row"""
        assert row < self.nrows
        return sum([r != 0 for r in self.data[row]])
        
    def popcnt_col(self, col: int) -> int:
        """ computes the hamming weight of a column"""
        assert col < self.ncols
        t = 0
        for j in range(self.nrows):
            t += self.data[j][col] != 0
        return t

    def __swap_rows(self, i: int, j: int) -> None:
        """ swap the rows i and j """
        assert i < self.nrows and j < self.nrows
        if i == j: return
        for k in range(self.ncols):
            tmp = self.data[i][k]
            self.data[i][k] = self.data[j][k]
            self.data[j][k] = tmp

    def __swap_cols(self, i: int, j: int) -> None:
        """ swap the cols i and j """
        assert i < self.ncols and j < self.ncols
        if i == j: return
        for k in range(self.nrows):
            tmp = self.data[k][i]
            self.data[k][i] = self.data[k][j]
            self.data[k][j] = tmp

    def __add__(self, B: 'Matrix'):
        return self.add(B)

    def __sub__(self, B: 'Matrix'):
        return self.sub(B)

    def __mul__(self, B: 'Matrix'):
        return self.mul(B)
   
    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        ret = ""
        for i in range(self.nrows):
            for j in range(self.ncols):
                ret += str(self.data[i][j]).rjust(3, ' ') + " "
            ret += "\n"
        return ret

if __name__ == "__main__":
    nc, nr, q, w = 10, 5, 2, 2
    A = Matrix(nr, nc, q)
    A.print()
    print()

    A.random()
    A.print()
    print()

    rank = A.gauß()
    A.print()
    print("rank", rank)

    e = Matrix(1, nc, q)
    e.random_row_with_weight(0, w)
    e.print()
    print()

    eT = e.transpose()
    eT.print()
    print()

    C = A.mul(eT)
    C.print()

