#!/usr/bin/env python3
""" super simple fq implementation. The only goal is to have zero dependencies """

from typing import Union
import random
import math


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

    def __add__(self, fq: Union['Fq', int]) -> 'Fq':
        if isinstance(fq, Fq): fq = fq.__value
        return self.add(fq)

    def __sub__(self, fq: Union['Fq', int]) -> 'Fq':
        if isinstance(fq, Fq): fq = fq.__value
        return self.sub(fq)

    def __mul__(self, fq: Union['Fq', int]) -> 'Fq':
        if isinstance(fq, Fq): fq = fq.__value
        return self.mul(fq)

    def __eq__(self, fq) -> bool:
        if isinstance(fq, Fq): fq = fq.__value
        return self.__value == fq
    
    def __lt__(self, fq) -> bool:
        if isinstance(fq, Fq): fq = fq.__value
        return self.__value < fq
    
    def __le__(self, fq) -> bool:
        if isinstance(fq, Fq): fq = fq.__value
        return self.__value <= fq
    
    def __gt__(self, fq) -> bool:
        if isinstance(fq, Fq): fq = fq.__value
        return self.__value > fq
    
    def __ge__(self, fq) -> bool:
        if isinstance(fq, Fq): fq = fq.__value
        return self.__value >= fq
    
    def __str__(self) -> str:
        return str(self.__value)


def rand_elements(rng, 
                  q, 
                  buffer: list,
                  l: int, 
                  MIN_VALUE=0, 
                  MAX_VALUE=0):
    """ generate a element mod Q
    :param rng: shake state 
    :param Q: prime 
    """
    if MAX_VALUE == 0:
        MAX_VALUE = q-1

    SPAN = MAX_VALUE - MIN_VALUE
    REQ_BITS = math.ceil(math.log2(SPAN))
    EL_MASK = (1 << REQ_BITS) - 1
    counter = 0
    while True:
        word = int(rng.read(8).hex(), 16)
        for _ in range(64//REQ_BITS):
            rnd_value = word & EL_MASK
            if (rnd_value <= SPAN):
                buffer[counter] = Fq(rnd_value + MIN_VALUE, q)
                counter += 1 
            if(counter >= l):
                return
            word >>= REQ_BITS

if __name__ == "__main__":
    from Crypto.Hash import SHAKE256
    l = [0 for _ in range(128)]
    seed = bytearray(32)
    rng = SHAKE256.new(seed)
    rand_elements(rng, 127, l, len(l))
    print(l)

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
