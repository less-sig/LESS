#!/usr/bin/env python3
from config import Xof
from fq import Fq

# NOTE: this holds the squeezing state, not the absorbing one. Consecutive
# calls to `randombytes` have to return consecutive bytes of the same
# stream, so the state must not be re-finalized on every call.
PLATFORM_CSPRNG_STATE = Xof.new()


def init_randombytes(buf: bytes):
    global PLATFORM_CSPRNG_STATE
    PLATFORM_CSPRNG_STATE = Xof.new(bytes(buf))


def randombytes(n: int) -> bytes:
    return PLATFORM_CSPRNG_STATE.read(n)


def csprng_initialize(seed: bytes):
    return Xof.new(bytes(seed))


def csprng_initialize_ds(seed: bytes, domain_seperator: int):
    """
    \\param seed[in]: seed which is feed into the prng
    \\param domain_seperator[in]: domain seperator as defined in the spec,
            which is fed into the prng after the seed
    """
    return Xof.new(bytes(seed) + domain_seperator.to_bytes(2, "little"))


def csprng_randombytes(state, n: int) -> bytes:
    return state.read(n)


def rand_range_impl(state, n: int, min_v: int) -> list[Fq]:
    span = 127 - min_v
    req_bits = 7
    el_mask = (1 << req_bits) - 1
    buffer = []
    while True:
        word = int.from_bytes(state.read(8), "little")
        for _ in range(64 // req_bits):
            rnd_value = word & el_mask
            if rnd_value < span:
                buffer.append(Fq(min_v + rnd_value))
                if len(buffer) >= n:
                    return buffer
            word >>= req_bits


def rand_range_q_state_elements(state, n: int) -> list[Fq]:
    return rand_range_impl(state, n, 0)


def fq_star_rnd_state_elements(state, n: int) -> list[Fq]:
    return rand_range_impl(state, n, 1)


def merge_exchange(buffer: list[int], state) -> list[int]:
    n = len(buffer)
    rand, ctr = 0, 0

    # Find the largest power of 2 less than n
    t = 1
    while t < n - t:
        t += t

    p_step = t
    while p_step > 0:
        q, r, d = t << 1, 0, p_step
        while True:
            q >>= 1
            for i in range(n - d):
                if (i & p_step) == r:
                    # Refill random bits if exhausted
                    if ctr == 0:
                        rand = int.from_bytes(state.read(8), "little")
                        ctr = 64

                    # Conditional swap
                    if rand & 1:
                        buffer[i], buffer[i + d] = buffer[i + d], buffer[i]

                    # Drain random bit and counter
                    rand >>= 1
                    ctr -= 1
            d = q - p_step
            r = p_step
            if q == p_step:
                break
        p_step >>= 1
    return buffer


def test():
    init_randombytes(bytes(32))
    a, b = randombytes(32), randombytes(32)
    assert a != b, "consecutive calls must not repeat the same bytes"
    init_randombytes(bytes(32))
    assert randombytes(64) == a + b

    # the domain separator is appended as two little endian bytes
    seed = bytes([7] * 16)
    a = csprng_randombytes(csprng_initialize_ds(seed, 0x0102), 32)
    assert a == csprng_randombytes(csprng_initialize(seed + b"\x02\x01"), 32)

    buffer = rand_range_q_state_elements(csprng_initialize(b"seed"), 256)
    assert all(v < 127 for v in buffer)
    buffer = fq_star_rnd_state_elements(csprng_initialize(b"seed"), 256)
    assert all(1 <= v < 127 for v in buffer)

    buffer = merge_exchange(list(range(128)), csprng_initialize(b"123"))
    assert sorted(buffer) == list(range(128))
    assert sum(buffer[i] == i for i in range(128)) < 5


if __name__ == "__main__":
    test()
