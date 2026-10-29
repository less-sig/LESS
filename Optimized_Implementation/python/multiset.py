#!/usr/bin/env python3
from config import Q_PAD
from vector import Vector


class Multiset(list):
    """ Multiset is actually a histogram """
    @staticmethod
    def default() -> 'Multiset':
        """ Zero Initialised """
        return Multiset([0] * Q_PAD)

    @staticmethod
    def from_row(row: Vector) -> 'Multiset':
        ret = Multiset.default()
        for v in row:
            ret[v] += 1
        return ret

    def partial_cmp(self, b: 'Multiset') -> int:
        """
        return: 0 if multiset(row1) == multiset(row2)
                 x if multiset(row1) > multiset(row2)
                -x if multiset(row1) < multiset(row2)
        """
        i = 0
        while i + 1 < len(self) and self[i] == b[i]:
            i += 1
        return b[i] - self[i]


def test():
    from fq import Fq
    c = Multiset.from_row(Vector([Fq(1)] * Q_PAD))
    assert c[0] == 0 and c[1] == Q_PAD and sum(c) == Q_PAD

    a, b = Multiset.default(), Multiset.default()
    assert a.partial_cmp(b) == 0
    b[0] = 1
    assert b.partial_cmp(a) < 0
    assert a.partial_cmp(b) > 0


if __name__ == "__main__":
    test()
