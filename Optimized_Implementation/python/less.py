#!/usr/bin/env python3
from config import (
    Hash, N, K, T, W, N8, NUM_KEYPAIRS, RREF_MAT_PACKEDBYTES, HASH_DIGEST_LENGTH,
    SEED_LENGTH_BYTES, PRIVATE_KEY_SEED_LENGTH_BYTES, SEED_TREE_MAX_PUBLISHED_BYTES,
    MAX_PUBLISHED_SEEDS, SIGN_PIVOT_REUSE_LIMIT, LESS_VERIFY_PIVOT_REUSE_LIMIT,
)
from matrix import Matrix, MatrixNormalized
from monomial import Monomial, Permutation
from prng import randombytes, csprng_initialize, csprng_randombytes
from tree import build_ggm, ggm_path, rebuild_ggm, seed_leaves
from utils import sample_challenge, verify


class PrivateKey:
    def __init__(self, compressed_sk: bytes) -> None:
        self.compressed_sk = compressed_sk


class PublicKey:
    """
    Public key: the first gen. matrix is shrunk to just a seed, all the
    others are stored in RREF form
    """
    def __init__(self, g0_seed: bytes, sf_g: list[bytes]) -> None:
        self.g0_seed = g0_seed
        self.sf_g = sf_g

    def bytes(self) -> bytes:
        return self.g0_seed + b"".join(self.sf_g)

    @staticmethod
    def from_bytes(pk: bytes) -> 'PublicKey':
        return PublicKey(pk[:SEED_LENGTH_BYTES],
                         [pk[SEED_LENGTH_BYTES + i * RREF_MAT_PACKEDBYTES:][:RREF_MAT_PACKEDBYTES]
                          for i in range(NUM_KEYPAIRS - 1)])


class Signature:
    def __init__(self, digest: bytes, salt: bytes, cf_monom_actions: list[bytes], seed_storage: bytes) -> None:
        self.digest = digest
        self.salt = salt
        self.cf_monom_actions = cf_monom_actions
        self.seed_storage = seed_storage


def expand_sk(sk: PrivateKey):
    """ expanding it onto private seeds """
    sk_shake_state = csprng_initialize(sk.compressed_sk)

    # Generating public code G_0
    g0_seed = csprng_randombytes(sk_shake_state, SEED_LENGTH_BYTES)

    # The first private key monomial is an ID matrix, no need for random
    # generation, hence NUM_KEYPAIRS-1
    private_monomial_seeds = [csprng_randombytes(sk_shake_state, PRIVATE_KEY_SEED_LENGTH_BYTES)
                              for _ in range(NUM_KEYPAIRS - 1)]
    return sk_shake_state, g0_seed, private_monomial_seeds


def expand_g0(g0_seed: bytes):
    g0_rref = Matrix.rand_from_seed(g0_seed, K, N - K)
    g0_column_pos = [i + K for i in range(N - K)]
    # generator_get_pivot_flags: everything is a pivot but the stored columns
    g0_initial_pivot_flags = [1] * N
    for i in range(K):
        g0_initial_pivot_flags[g0_column_pos[i]] = 0
    return Matrix.from_compressed_rref(g0_column_pos, g0_rref), g0_initial_pivot_flags


def ephemeral_monomial(round_seed: bytes, salt: bytes, i: int) -> Monomial:
    """
    the ephemeral monomial is drawn from the round seed, the salt and the
    round index, see `monomial_sample_salt`
    """
    mu_input = round_seed + salt + i.to_bytes(2, "little")
    return Monomial.rand(csprng_initialize(mu_input), N)


def hash_rows(hasher, a: MatrixNormalized):
    for row in a.rows:
        hasher.update(bytes(row))


def less_keygen() -> tuple[PrivateKey, PublicKey]:
    """ LESS_keygen """
    # generating private key from a single seed
    sk = PrivateKey(randombytes(PRIVATE_KEY_SEED_LENGTH_BYTES))
    _, g0_seed, private_monomial_seeds = expand_sk(sk)
    tmp_full_g, _ = expand_g0(g0_seed)

    sf_g = []
    was_pivot_column = [0] * N
    # note that the first "keypair" is just the public generator G_0, stored
    # as a seed and the identity matrix (not stored)
    for seed in private_monomial_seeds:
        is_pivot_column = [0] * N
        # expand inverse monomial from seed
        private_q_inv = Monomial.rand_from_seed(seed, N)
        # NOTE: the inversion is not implemented in constant time.
        private_q = private_q_inv.inv()
        result_g = tmp_full_g * private_q
        result_g.rref(is_pivot_column, was_pivot_column, K)
        sf_g.append(result_g.compress(is_pivot_column))
    return sk, PublicKey(g0_seed, sf_g)


def less_sign(m: bytes, sk: PrivateKey) -> tuple[Signature, int]:
    # Private key expansion
    sk_shake_state, g0_seed, private_monomial_seeds = expand_sk(sk)

    # generate the salt from a TRNG
    salt = randombytes(HASH_DIGEST_LENGTH)

    # Ephemeral monomial generation
    ephem_monomials_seed = csprng_randombytes(sk_shake_state, SEED_LENGTH_BYTES)

    # create the prng for the "blinding" monomials for the canonical form computation
    cf_seed = csprng_randombytes(sk_shake_state, SEED_LENGTH_BYTES)
    cf_shake_state = csprng_initialize(cf_seed)

    seed_tree = build_ggm(ephem_monomials_seed, salt)
    linearized_rounds_seeds = seed_leaves(seed_tree)

    # Public G_0 expansion
    full_g0, g0_initial_pivot_flags = expand_g0(g0_seed)

    hasher = Hash.new()
    pi_tilde = []
    i = 0
    while i < T:
        mu_tilde = ephemeral_monomial(linearized_rounds_seeds[i], salt, i)
        G0 = full_g0 * mu_tilde

        is_pivot_column = [0] * N
        permuted_pivot_flag = [0] * N
        for t in range(N):
            permuted_pivot_flag[mu_tilde.perms[t]] = g0_initial_pivot_flags[t]

        if not G0.rref(is_pivot_column, permuted_pivot_flag, SIGN_PIVOT_REUSE_LIMIT):
            return None, 0
        # just copy the non IS
        a_j = MatrixNormalized.from_information_set(G0, is_pivot_column)
        pi = Permutation.from_information_set_composition(mu_tilde, is_pivot_column)

        a_j.blind(cf_shake_state)
        if a_j.cf() == 0:
            # this gives us 2**64 retries. That should be more than enough.
            s = linearized_rounds_seeds[i]
            v = (int.from_bytes(s[:8], "little") + 1) % 2**64
            linearized_rounds_seeds[i] = v.to_bytes(8, "little") + s[8:]
            continue

        pi_tilde.append(pi)
        hash_rows(hasher, a_j)
        i += 1

    hasher.update(m)
    hasher.update(salt)
    digest = hasher.digest()

    fixed_weight_string = sample_challenge(digest)
    indices_to_publish = [int(f != 0) for f in fixed_weight_string]
    seeds = ggm_path(seed_tree, indices_to_publish)

    cf_monom_actions = []
    for i in range(T):
        if fixed_weight_string[i] != 0:
            j = fixed_weight_string[i]
            Q_to_multiply = Monomial.rand_from_seed(private_monomial_seeds[j - 1], N)
            mono_action = Permutation.from_composition_action(Q_to_multiply, pi_tilde[i])
            cf_monom_actions.append(mono_action.bytes())

    seed_storage = seeds.ljust(SEED_TREE_MAX_PUBLISHED_BYTES, b"\x00")
    return Signature(digest, salt, cf_monom_actions, seed_storage), len(seeds) // SEED_LENGTH_BYTES


def less_verify(pk: PublicKey, m: bytes, sig: Signature) -> bool:
    """
    NOTE: non-constant time
    :param pk[in]: public key
    :param m[in]: message for which a signature was computed
    :param sig[in]: signature
    :return false: on failure
             true: on success
    """
    fixed_weight_string = sample_challenge(sig.digest)
    published_seed_indexes = [int(f != 0) for f in fixed_weight_string]

    seed_tree = rebuild_ggm(published_seed_indexes, sig.seed_storage, sig.salt)
    if seed_tree is None:
        return False
    linearized_rounds_seeds = seed_leaves(seed_tree)

    g0_full, g0_initial_pivot_flags = expand_g0(pk.g0_seed)
    gs = [Matrix.expand(sf_g, [0] * N) for sf_g in pk.sf_g]

    hasher = Hash.new()
    employed_monoms = 0
    for i in range(T):
        is_pivot_column = [0] * N
        if fixed_weight_string[i] == 0:
            # same ephemeral monomial derivation as in `less_sign`
            mu_tilde = ephemeral_monomial(linearized_rounds_seeds[i], sig.salt, i)
            g_prime = g0_full * mu_tilde

            permuted_pivot_flags = [0] * N
            for t in range(N):
                permuted_pivot_flags[mu_tilde.perms[t]] = g0_initial_pivot_flags[t]
        else:
            gi_initial_pivot_flags = [0] * N
            key_pos = fixed_weight_string[i] - 1
            action = sig.cf_monom_actions[employed_monoms]

            if not Permutation.check_canonical_action(action):
                return False

            permuted_pivot_flags = [0] * N
            g_prime = Matrix.unpack_coset_rep(gs[key_pos], action, gi_initial_pivot_flags,
                                              permuted_pivot_flags)
            employed_monoms += 1

        if not g_prime.rref(is_pivot_column, permuted_pivot_flags, LESS_VERIFY_PIVOT_REUSE_LIMIT):
            return False

        # copy the non IS
        a_i = MatrixNormalized.from_information_set(g_prime, is_pivot_column)

        # compute the canonical form and exit on any failure
        if a_i.cf() == 0:
            return False

        # NOTE: as we increase the size of the `normalized_IS_t`
        # we need to hash the values row by row.
        hash_rows(hasher, a_i)

    hasher.update(m)
    hasher.update(sig.salt)
    recomputed_digest = hasher.digest()

    # NOTE: does not have to be constant time.
    return verify(recomputed_digest, sig.digest) == 0


def signature_size(leaves: int) -> int:
    """ LESS_SIGNATURE_SIZE """
    return HASH_DIGEST_LENGTH * 2 + N8 * W + leaves * SEED_LENGTH_BYTES + 1


def crypto_sign_keypair() -> tuple[bytes, bytes]:
    sk, pk = less_keygen()
    return pk.bytes(), sk.compressed_sk


def crypto_sign(m: bytes, sk: bytes) -> bytes:
    sig, leaves = less_sign(m, PrivateKey(sk))
    sm = m + sig.digest + sig.salt + b"".join(sig.cf_monom_actions) + sig.seed_storage
    return sm[:len(m) + signature_size(leaves) - 1] + bytes([leaves])


def crypto_sign_open(sm: bytes, pk: bytes):
    """
    verifying a signed message sm[0],sm[1],...,sm[smlen-1]
    under public key pk[0],pk[1],...
    and producing original message m[0],m[1],...,m[*mlen-1]
    """
    num_seeds_published = sm[-1]
    if num_seeds_published > MAX_PUBLISHED_SEEDS:
        return None
    mlen = len(sm) - signature_size(num_seeds_published)
    if mlen < 0:
        return None

    m, sig = sm[:mlen], sm[mlen:]
    actions = sig[2 * HASH_DIGEST_LENGTH:2 * HASH_DIGEST_LENGTH + N8 * W]
    sig = Signature(sig[:HASH_DIGEST_LENGTH],
                    sig[HASH_DIGEST_LENGTH:2 * HASH_DIGEST_LENGTH],
                    [actions[i * N8:(i + 1) * N8] for i in range(W)],
                    # only the published seeds are stored, the rest of the buffer stays zero
                    sig[2 * HASH_DIGEST_LENGTH + N8 * W:-1].ljust(SEED_TREE_MAX_PUBLISHED_BYTES, b"\x00"))
    return m if less_verify(PublicKey.from_bytes(pk), m, sig) else None


def test():
    from prng import init_randombytes

    def keypair_and_signature(msg: bytes):
        init_randombytes(bytes(range(48)))
        sk, pk = less_keygen()
        sig, leaves = less_sign(msg, sk)
        assert 0 < leaves <= MAX_PUBLISHED_SEEDS, f"{leaves} published leaves is out of range"
        return pk, sig

    msg = b"abc"
    pk, sig = keypair_and_signature(msg)
    assert less_verify(pk, msg, sig), "a valid signature was rejected"

    # a different message must not verify
    assert not less_verify(pk, b"ghi", sig)

    sig.digest = bytes([sig.digest[0] ^ 1]) + sig.digest[1:]
    assert not less_verify(pk, msg, sig), "a flipped digest bit was accepted"


if __name__ == "__main__":
    test()
