#!/usr/bin/env python3
from config import T, W, NUM_KEYPAIRS, HASH_DIGEST_LENGTH
from prng import csprng_initialize, csprng_randombytes


def verify(a: bytes, b: bytes) -> int:
    """Compare two arrays for equality in constant time.
    :returns 0 if the byte arrays are equal, 1 otherwise
    """
    r = 0
    for x, y in zip(a, b):
        r |= x ^ y
    return int(r != 0)


def sample_challenge(digest: bytes) -> list[int]:
    """ Expands a digest into a fixed weight string with elements in Z_{NUM_KEYPAIRS}.
    :return fixed_weight_string: array of length T, with hamming weight W,
        where values > 0 are between [1, NUM_KEYPAIRS-1)
    """
    # Initialize CSPRNG from the digest
    shake_state = csprng_initialize(digest)

    keypair_bits = (NUM_KEYPAIRS - 1).bit_length()
    position_bits = (T - 1).bit_length()

    # Clear the first T-W entries to zero
    fixed_weight_string = [0] * T

    # NOTE: `rnd_buf` and `c` are shared between the two loops below, the
    # shuffle keeps consuming the bits the value loop left over.
    rnd_buf, c = 0, 0

    def next_bits(bits: int) -> int:
        nonlocal rnd_buf, c
        if c == 0:
            rnd_buf = int.from_bytes(csprng_randombytes(shake_state, 8), "little")
            c = 64 // bits
        v = rnd_buf & ((1 << bits) - 1)
        rnd_buf >>= bits
        c -= 1
        return v

    if NUM_KEYPAIRS != 2:
        # Rejection sampling for values in [1, NUM_KEYPAIRS-1)
        for i in range(T - W, T):
            value = next_bits(keypair_bits)
            while value >= NUM_KEYPAIRS - 1:
                value = next_bits(keypair_bits)
            fixed_weight_string[i] = value + 1
    else:
        # NUM_KEYPAIRS == 2: the last W values are 1
        fixed_weight_string[T - W:] = [1] * W

    # TODO: 
    # Fisher-Yates-like shuffle into the first T-W zero positions
    for p in range(T - W, T):
        pos = next_bits(position_bits)
        while pos > p:
            pos = next_bits(position_bits)

        # Swap
        fixed_weight_string[p], fixed_weight_string[pos] = fixed_weight_string[pos], fixed_weight_string[p]
    return fixed_weight_string


def test():
    fws = sample_challenge(bytes([0x5a] * HASH_DIGEST_LENGTH))
    assert sum(v != 0 for v in fws) == W, "the challenge must have hamming weight W"
    assert all(v < NUM_KEYPAIRS for v in fws)
    assert fws == sample_challenge(bytes([0x5a] * HASH_DIGEST_LENGTH))
    assert fws != sample_challenge(bytes([0x01] * HASH_DIGEST_LENGTH))
    assert verify(b"abc", b"abc") == 0 and verify(b"abc", b"abd") == 1


if __name__ == "__main__":
    test()
