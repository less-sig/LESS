#!/usr/bin/env python3
""" super simple fq implementation. The only goal is to have zero dependencies """

class Fq():
    q = 127
    def __init__(self, value: int) -> None:
        self.__value = value % Fq.q

    def add(self, value: int):
        self.__value = (self.__value + value) % Fq.q
        return self

    def sub(self, value: int):
        self.__value = (self.__value + (Fq.q - value)) % Fq.q
        return self

    def mul(self, value: int):
        self.__value = (self.__value * value) % Fq.q
        return self

    def __add__(self, fq: 'Fq'):
        return self.add(fq.__value)

    def __sub__(self, fq: 'Fq'):
        return self.sub(fq.__value)

    def __mul__(self, fq: 'Fq'):
        return self.mul(fq.__value)
