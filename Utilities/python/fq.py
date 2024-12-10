#!/usr/bin/env python3
""" super simple fq implementation. The only goal is to have zero dependencies """

import random


class Fq():
    def __init__(self, value: int = 0, q: int = 127) -> None:
        self.q = q
        self.__value = value % self.q
    
    def get(self) -> int:
        """ returns the internal value as an int """
        return self.__value

    def set(self, value: int):
        """ sets the internal value to an int """
        self.__value = value % self.q

    def add(self, value: int) -> 'Fq':
        """  """
        ret = Fq((self.__value + value) % self.q, self.q)
        return ret

    def sub(self, value: int) -> 'Fq':
        """  """
        ret = Fq((self.__value + (self.q - value)) % self.q, self.q)
        return ret

    def mul(self, value: int) -> 'Fq':
        """  """
        ret = Fq((self.__value * value) % self.q, self.q)
        return ret

    def random(self, lower: int = 0, upper: int = 0) -> 'Fq': 
        """ generates a random element """
        if upper == 0:
            upper = self.q-1
        assert (lower <= upper)
        self.__value = random.randint(lower, upper)
        return self

    def __add__(self, fq: 'Fq') -> 'Fq':
        return self.add(fq.__value)

    def __sub__(self, fq: 'Fq') -> 'Fq':
        return self.sub(fq.__value)

    def __mul__(self, fq: 'Fq') -> 'Fq':
        return self.mul(fq.__value)

    def __eq__(self, fq) -> bool:
        return self.__value == fq.__value
    
    def __str__(self) -> str:
        return str(self.__value)

if __name__ == "__main__":
    # just some simple test for q==127
    one = Fq(1)
    t = Fq(126)
    two = one + one
    assert two.get() == 2

    zero = t + one
    assert zero.get() == 0

    for i in range(0, 126):
        t = Fq(i) 
        x = t * one 
        assert t.get() == x.get()
