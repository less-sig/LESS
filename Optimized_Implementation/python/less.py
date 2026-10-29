#!/usr/bin/env python3
from Crypto.Hash import SHAKE256
from typing import Tuple

from .drbg import NIST_KAT_DRBG
from monomial import Monomial
from matrix import Matrix
from fq import Fq
from cf import CF


class less:
    """
    """
    N = 252
    K = N//2
    NK = N - K
    N8 = N // 8
    Q = 127
    T = 247
    W = 30
    NUM_KEYPAIRS = 2 
    
    HASH_DIGEST_LENGTH = 32
    SEED_LENGTH_BYTES = 16
    PRIVATE_KEY_SEED_LENGTH_BYTES = 2*SEED_LENGTH_BYTES
    SEED_TREE_MAX_PUBLISHED_BYTES = 1488

    class pk:
        def __init__(self) -> None:
            """
            """
            self.G_0_seed = bytearray(less.SEED_LENGTH_BYTES)
            self.SF_G = bytearray(less.K*less.NK + less.N)

    class sk:
        def __init__(self) -> None:
            """
            """
            self.compressed_sk = bytearray(less.PRIVATE_KEY_SEED_LENGTH_BYTES)
            self.G_0_seed = bytearray(less.SEED_LENGTH_BYTES)

    class sig:
        def __init__(self) -> None:
            """
            """
            self.cf_monom_actions = [bytearray(less.N8)] * less.W
            self.tree_salt = bytearray(less.HASH_DIGEST_LENGTH)
            self.digest = bytearray(less.HASH_DIGEST_LENGTH)
            self.seed_storage = bytearray(less.SEED_TREE_MAX_PUBLISHED_BYTES + 1)


    def __init__(self, rbg=NIST_KAT_DRBG) -> None:
        """
        """
        self.rbg = rbg
    
    def set_random(self, rbg):
        """ Set the key material RBG."""
        self.rbg        =   rbg 

    def shake256(self, x, l):
        """ shake256s(x, l): Internal hook."""
        return bytearray(SHAKE256.new(x).read(l))

    def keygen(self, seed_sk=None):
        sk = less.sk()
        pk = less.pk()

        if seed_sk == None:
            seed_sk = self.rbg(less.PRIVATE_KEY_SEED_LENGTH_BYTES)

        sk.compressed_sk = bytearray(seed_sk.random_bytes(less.PRIVATE_KEY_SEED_LENGTH_BYTES))
        sk_shake_state = SHAKE256.new(sk.compressed_sk)
        
        private_monomial_seeds = []*(less.NUM_KEYPAIRS - 1)
        for i in range(less.NUM_KEYPAIRS - 1):
            private_monomial_seeds[i] = sk_shake_state.read(less.PRIVATE_KEY_SEED_LENGTH_BYTES)

        # public stuff
        sk.G_0_seed = bytearray(seed_sk.random_bytes(less.SEED_LENGTH_BYTES))
        pk.G_0_seed = sk.G_0_seed

        G0_rref = Matrix(less.K, less.NK, less.Q).random_from_seed(sk.G_0_seed)
        G0_rred_pivots = [less.K+i for i in range(less.NK)]

        tmp_full_G = Matrix(less.K, less.N, less.Q)
        self.generator_rref_expand(tmp_full_G, G0_rref, G0_rred_pivots)
        
        for i in range(less.NUM_KEYPAIRS - 1):
            private_Q_inv = Monomial(self.N, self.Q).random_from_seed(private_monomial_seeds[i])
            private_Q = private_Q_inv.inv()
            result_G = private_Q.apply(tmp_full_G)
            pivot_cols = [0 for _ in range(self.N)]
            result_G.rref(pivot_cols)
            self.compress_rref(pk.SF_G, result_G, pivot_cols)
        return sk, pk

    def sign(self, 
             sk, 
             m: bytearray,
             mlen: int):
        
        ephermal_monomials_seed = self.rbg(self.SEED_LENGTH_BYTES)
        cf_seed = self.rbg(self.SEED_LENGTH_BYTES)
        for i in range(self.T):
            G0 = mu_tilde(full_G0)
        pass

    def open(self):
        pass

    def compress_rref(self, 
                      compressed: bytearray, 
                      G: Matrix,
                      pivot_cols: list[int]):
        for i in range(self.N//8):
            t = (pivot_cols[8*i + 0] << 0) | \
                (pivot_cols[8*i + 1] << 1) | \
                (pivot_cols[8*i + 2] << 2) | \
                (pivot_cols[8*i + 3] << 3) | \
                (pivot_cols[8*i + 4] << 4) | \
                (pivot_cols[8*i + 5] << 5) | \
                (pivot_cols[8*i + 6] << 6) | \
                (pivot_cols[8*i + 7] << 7)
            compressed[i] = t

        # TODO this is only correct for CAT I
        compressed[self.N//8] = (pivot_cols[self.N-4]<<0) |\
                                (pivot_cols[self.N-3]<<1) |\
                                (pivot_cols[self.N-2]<<2) |\
                                (pivot_cols[self.N-1]<<3)
        i = self.N//8 + 1
        state = 0
        for r in range(self.K):
            for c in range(self.N):
                if not pivot_cols[c]:
                    if state == 0:
                        compressed[i] = G[r][c].get()
                    else:
                        compressed[i] |= (G[r][c].get() << (8-state))
                        i = i + 1 
                        if state < 7:
                            compressed[i] |= (G[r][c].get() << state)
                    if state < 7:
                        state += 1 
                    else:
                        state = 0

    def generator_rref_expand(self, 
                              full: Matrix, 
                              compact: Matrix,
                              compact_pivots: list):
        placed_dense_cols = 0 
        for col_idx in range(less.N):
            if (placed_dense_cols < less.NK) and (col_idx == compact_pivots[placed_dense_cols]):
                for row_idx in range(less.K):
                    full[row_idx, col_idx] = compact[row_idx, placed_dense_cols]
            else:
                for row_idx in range(less.K):
                    full[row_idx, col_idx] = Fq(row_idx == (col_idx - placed_dense_cols), less.Q)
                
        pass

    def build_full_generator_matrix(self, A: Matrix) -> Matrix:
        """
        """
        G = Matrix(self.K, self.N, self.Q).zero()
        for i in range(self.K):
            G[i, i] = 1

        for i in range(self.K):
            for j in range(self.NK):
                G[i, j+self.K] = A[i, j]
        return G

    def compress(self, A: Monomial) -> bytearray:
        """
        compress a canonical action
        """
        ret = bytearray(self.N)
        for i in range(self.N):
            ret[A.perm[i]] = 1
        return ret

    def cf(self, G: Matrix) -> Tuple[bool, Matrix]:
        """
        :param G
        :return canonical form of a matrix
        """
        return CF(G)

    def blind(self, 
              prng,
              G: Matrix):
        """ randomly chooses two monomial matrices
        :param prng
        :param G non-IS part of a generator matrix 
        :return left * G * right
        """
        left = Monomial(self.NK, self.Q).random_from_rng(prng)
        right = Monomial(self.NK, self.Q).random_from_rng(prng)
        G = right.apply_right(G)
        return left.apply_left(G)

