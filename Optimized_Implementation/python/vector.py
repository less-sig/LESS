#!/usr/bin/env python3
from fq import Fq
from prng import rand_range_q_state_elements, fq_star_rnd_state_elements


class Vector(list):
    @staticmethod
    def default(n: int) -> 'Vector':
        """ Zero Initialised """
        return Vector([Fq(0)] * n)

    @staticmethod
    def rand(state, n: int) -> 'Vector':
        """ sample a random vector """
        return Vector(rand_range_q_state_elements(state, n))

    @staticmethod
    def rand_star(state, n: int) -> 'Vector':
        """ sample a random vector """
        return Vector(fq_star_rnd_state_elements(state, n))

    def add(self, b: 'Vector') -> 'Vector':
        """ c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..N """
        return Vector(x + y for x, y in zip(self, b))

    def sub(self, b: 'Vector') -> 'Vector':
        """ c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..N """
        return Vector(x - y for x, y in zip(self, b))

    def mul(self, b: 'Vector') -> 'Vector':
        """ c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..N """
        return Vector(x * y for x, y in zip(self, b))

    def div(self, b: 'Vector') -> 'Vector':
        """ c[i] = a[i]*b[i]^{-1} mod q,  a[i],b[i] < 127, for i in 0..N """
        return Vector(x * y.inv() for x, y in zip(self, b))

    def scalar_mul(self, a: Fq) -> 'Vector':
        """ c[i] = c[i]*a mod q,  a[i],b[i] < 127, for i in 0..N """
        return Vector(x * a for x in self)

    def __add__(self, b: 'Vector') -> 'Vector':
        return self.add(b)

    def __sub__(self, b: 'Vector') -> 'Vector':
        return self.sub(b)

    def __mul__(self, b) -> 'Vector':
        return self.scalar_mul(b) if isinstance(b, int) else self.mul(b)

    def __neg__(self) -> 'Vector':
        return Vector(Fq(0) - x for x in self)

    __rmul__ = __mul__
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__

    def acc(self) -> Fq:
        """ computes: sum_{i=0}^{N} row[i] mod q,  row[i] < 127, for i in 0..N """
        return Fq(sum(self, Fq(0)))

    def inv_non_ct(self) -> 'Vector':
        """ c[i] = a[i]^{-1} mod q, a[i] < 127 for i in 0..N """
        return Vector(x.inv() for x in self)

    def acc_inv(self) -> Fq:
        """ computes: sum_{i=0}^{N} a[i]^{-1} mod q, a[i] < 127 for i in 0..N """
        return self.inv_non_ct().acc()

    def all_same(self) -> bool:
        """ c[i] == c[j] for all i < j < n """
        return all(x == self[0] for x in self)

    def count_zero(self) -> int:
        """ number of zeros in row """
        return self.count(Fq(0))


def test():
    v = lambda l: Vector(Fq(x) for x in l)
    assert v([1, 2, 3]).mul(v([3, 2, 1])) == [3, 4, 3]
    assert v([1, 126]).scalar_mul(Fq(2)) == [2, 125]
    assert v([126, 2]).acc() == 1
    assert v([1, 2]).acc_inv() == (1 + 64)
    assert v([0, 1, 2]).inv_non_ct() == [0, 1, 64]
    assert v([3, 3, 3]).all_same() and not v([3, 3, 2]).all_same()
    assert v([0, 1, 0]).count_zero() == 2
    assert v([1, 126]) + v([1, 2]) == [2, 1]
    assert v([1, 0]) - v([2, 1]) == [126, 126]
    assert v([2, 3]) * v([3, 2]) == [6, 6] and v([2, 3]) * Fq(2) == [4, 6]
    assert v([6, 6]) / v([3, 2]) == [2, 3] and v([4, 6]) / Fq(2) == [2, 3]
    assert -v([1, 0]) == [126, 0]
    w = v([1, 2])
    w += v([1, 1])
    assert w == [2, 3] and isinstance(w, Vector)
    from prng import csprng_initialize
    assert all(isinstance(x, Fq) for x in Vector.rand_star(csprng_initialize(b""), 8))


if __name__ == "__main__":
    test()
