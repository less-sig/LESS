#!/usr/bin/env python3
from Crypto.Hash import SHAKE256

from .drbg import NIST_KAT_DRBG
from monomial import Monomial
from matrix import Matrix



class less:
    N = 256
    K = 128
    NK = N - K
    N8 = N // 8
    Q = 127
    T = 247
    W = 30
    NUM_KEYPAIRS = 2 
    
    HASH_DIGEST_LENGTH = 32
    SEED_LENGTH_BYTES = 16
    PRIVATE_KEY_SEED_LENGTH_BYTES = 2*SEED_LENGTH_BYTES
    SEED_TREE_MAX_PUBLISHED_BYTES = 1472

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
        
        for i in range(less.NUM_KEYPAIRS):
            Monomial 

    def sign(self):
        pass

    def open(self):
        pass

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
        """
        ret = bytearray(self.N)
        for i in range(self.N):
            ret[A.perm[i]] = 1
        return ret
