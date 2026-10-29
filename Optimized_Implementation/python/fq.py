#!/usr/bin/env python3
""" super simple fq implementation. The only goal is to have zero dependencies """
from config import Q


class Fq(int):
    def __new__(cls, value: int = 0):
        return super().__new__(cls, value % Q)

    def get(self) -> int:
        """ returns the internal value as an int """
        return int(self)

    def __add__(self, b) -> 'Fq':
        return ELEMENTS[int.__add__(self, b) % Q]

    def __sub__(self, b) -> 'Fq':
        return ELEMENTS[int.__sub__(self, b) % Q]

    def __rsub__(self, b) -> 'Fq':
        return ELEMENTS[int.__rsub__(self, b) % Q]

    def __mul__(self, b) -> 'Fq':
        return ELEMENTS[int.__mul__(self, b) % Q]

    def __truediv__(self, b) -> 'Fq':
        return self * Fq(b).inv()

    def __rtruediv__(self, b) -> 'Fq':
        return Fq(b) * self.inv()

    def __neg__(self) -> 'Fq':
        return ELEMENTS[-int(self) % Q]

    __radd__ = __add__
    __rmul__ = __mul__

    def inv(self) -> 'Fq':
        """ a^{-1} mod q """
        return FQ127_INV_TABLE[self]


ELEMENTS = [Fq(i) for i in range(Q)]
# NOTE: we define 0^{-1} := 0.
FQ127_INV_TABLE = [Fq(pow(i, Q - 2, Q)) for i in range(Q + 1)]


if __name__ == "__main__":
    from Crypto.Hash import SHAKE256
    from prng import rand_range_q_state_elements
    l = [0 for _ in range(128)]
    seed = bytearray(32)
    rng = SHAKE256.new(seed)
    l = rand_range_q_state_elements(rng, len(l))
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

    for i in range(1, 127):
        assert (Fq(i) * Fq(i).inv()).get() == 1
        assert Fq(i) / Fq(i) == 1 and Fq(i) - Fq(i) == 0 and -Fq(i) + Fq(i) == 0
        assert 1 - Fq(i) == Fq(1 - i) and isinstance(1 - Fq(i), Fq)
