#!/usr/bin/env python3
from config import K, N, N8
from prng import csprng_initialize, merge_exchange
from vector import Vector


class Monomial:
    def __init__(self, coeffs: Vector, perms: list[int]) -> None:
        self.coeffs = coeffs
        self.perms = perms

    @staticmethod
    def default(n: int) -> 'Monomial':
        """ creates a new Monomial, initialized with as id """
        return Monomial(Vector.default(n), list(range(n)))

    @staticmethod
    def rand(state, n: int) -> 'Monomial':
        """
        NOTE: a monomial matrix has to be invertible, so the coefficients are
        drawn from F_q^*.
        """
        coeffs = Vector.rand_star(state, n)
        return Monomial(coeffs, merge_exchange(list(range(n)), state))

    @staticmethod
    def rand_from_seed(seed: bytes, n: int) -> 'Monomial':
        return Monomial.rand(csprng_initialize(seed), n)

    def __mul__(self, r: 'Monomial') -> 'Monomial':
        """ l = l*r """
        # `perms[i]` is the destination of column i, so applying `l` and
        # then `r` sends column i to `r.perms[l.perms[i]]`.
        return Monomial(Vector(c * r.coeffs[p] for c, p in zip(self.coeffs, self.perms)),
                        [r.perms[p] for p in self.perms])

    def inv(self) -> 'Monomial':
        """ l = r**-1 """
        ret = Monomial.default(len(self.perms))
        for i, t in enumerate(self.perms):
            ret.perms[t] = i
            ret.coeffs[t] = self.coeffs[i].inv()
        return ret

    def __eq__(self, b) -> bool:
        return self.perms == b.perms and self.coeffs == b.coeffs

    def __str__(self) -> str:
        return f"perms: {self.perms}\ncoeffs: {self.coeffs}"


class Permutation:
    def __init__(self, perms: list[int]) -> None:
        self.perms = perms

    @staticmethod
    def from_information_set_composition(mu_tilde: Monomial, is_pivot_column: list[int]) -> 'Permutation':
        """
        translation of monomial_compose_action_information_set
        NOTE: Only the permutation is computed, as this is the only thing we need
        since the adaption of canonical forms.
        :param mu_tilde[in]: input monomial matrix.
        :param is_pivot_column[in]: information set, where a 1 in the array
                symbolizes that the corresponding columns is a pivot column.
        :return pi_tilde: mu_tilde^{-1} * information set
        """
        row_idx = mu_tilde.inv().perms
        return Permutation([row_idx[c] for c in range(len(row_idx)) if is_pivot_column[c] == 1])

    @staticmethod
    def from_composition_action(Q_in: Monomial, in_: 'Permutation') -> 'Permutation':
        """
        NOTE: translation of 'monomial_compose_action'
        composes a compactly stored action of a monomial on an IS with a regular
        monomial.
        NOTE: Only the permutation is computed, as this is the only thing we need
        since the adaption of canonical forms.
        """
        reverse_Q = Q_in.inv().perms
        return Permutation([reverse_Q[p] for p in in_.perms])

    def bytes(self) -> bytes:
        """ compress a permutation to a byte string """
        ret = bytearray(N8)
        for p in self.perms:
            ret[p // 8] ^= 1 << (p % 8)
        return bytes(ret)

    @staticmethod
    def check_canonical_action(b: bytes) -> bool:
        """
        checks if the given (N+7)/8 bytes are a valid canonical form action,
        e.g. this function checks if exactly `K` ones are in the input.
        """
        w = int.from_bytes(b, "little")
        return w.bit_length() <= N and w.bit_count() == K

    def __str__(self) -> str:
        return f"perms: {self.perms}\n"

def test():
    n = 32
    s = csprng_initialize(bytes([1] * 32))
    a = Monomial.rand(s, n)
    assert a.inv().inv() == a
    assert sorted(a.perms) == list(range(n)) and 0 not in a.coeffs
    assert Monomial.rand_from_seed(bytes([9] * 32), n) == Monomial.rand(csprng_initialize(bytes([9] * 32)), n)

    # m * m^-1 must be the identity monomial
    t = a * a.inv()
    assert t.perms == list(range(n)) and t.coeffs == [1] * n

    # pi_tilde must select exactly the K pivot columns
    mu_tilde = Monomial.rand_from_seed(bytes([2] * 32), 2 * K)
    pi = Permutation.from_information_set_composition(mu_tilde, [1] * K + [0] * K)
    assert all(mu_tilde.perms[pi.perms[j]] == j for j in range(K))

    # composing with the identity monomial must leave the action unchanged
    out = Permutation.from_composition_action(Monomial.default(2 * K), Permutation(list(range(K))))
    assert out.perms == list(range(K))
    assert Permutation.check_canonical_action(out.bytes())


# just some simple tests
if __name__ == "__main__":
    from matrix import Matrix
    test()
    n, r, q = 10, 5, 2
    m1 = Monomial.rand_from_seed(bytes(32), n)
    print(m1)

    G1 = Matrix(r, n, q).id()
    G2 = G1 * m1
    print(G2)
