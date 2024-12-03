#!/usr/bin/env python3
""" super simple permutation implementation. 
The only goal is to have zero dependencies """

import random
from matrix import Matrix

class Permutation:
    """this implementation tries to mimic the sage `Permutation` implementation """
    def __init__(self, n: int) -> None:
        assert n > 0
        self.n = n 
        self.__data = [i for i in range(n)]
    
    def random(self) -> 'Permutation':
        for i in range(self.n):
            pos = random.randint(0, self.n - i - 1)
            tmp = self.__data[i]
            self.__data[i] = self.__data[i + pos]
            self.__data[i + pos] = tmp

        return self
    
    def to_matrix(self, q=127) -> Matrix:
        m = Matrix(self.n, self.n, q).zero()
        for i in range(self.n):
            j = self.__data[i]
            m[j, i] = 1
        return m

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        ret = ""
        for i in range(self.n):
            ret += str(self.__data[i]) + " "
        return ret

if __name__ == "__main__":
    n = 3 
    p = Permutation(n).random()
    print(p)
    m = p.to_matrix()
    print(m)
