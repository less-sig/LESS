#!/usr/bin/env python3
""" super simple fq implementation. The only goal is to have zero dependencies """

import random


class Fq():
    q = 127
    def __init__(self, value: int) -> None:
        self.__value = value % Fq.q
    
    def get(self) -> int:
        return self.__value 

    def add(self, value: int):
        self.__value = (self.__value + value) % Fq.q
        return self

    def sub(self, value: int):
        self.__value = (self.__value + (Fq.q - value)) % Fq.q
        return self

    def mul(self, value: int):
        self.__value = (self.__value * value) % Fq.q
        return self

    def random(self) -> 'Fq': 
        """ generates a random element """
        self.__value = random.randint(0, Fq.q-1)
        return self

    def __add__(self, fq: 'Fq') -> 'Fq':
        return self.add(fq.__value)

    def __sub__(self, fq: 'Fq') -> 'Fq':
        return self.sub(fq.__value)

    def __mul__(self, fq: 'Fq') -> 'Fq':
        return self.mul(fq.__value)

    def __eq__(self, fq) -> bool:
        return self.__value == fq.__value


if __name__ == "__main__":
    # just some simple test
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
