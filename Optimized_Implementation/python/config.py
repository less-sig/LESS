#!/usr/bin/env python3
import os
from Crypto.Hash import SHAKE128, SHAKE256, SHA3_256, SHA3_384, SHA3_512

Q = 127
Q_PAD = 128

CATEGORIES = {
    "cat1": (252, 126, 16, 25, SHAKE128, SHA3_256),
    "cat3": (400, 200, 24, 51, SHAKE256, SHA3_384),
    "cat5": (548, 274, 32, 79, SHAKE256, SHA3_512),
}

TARGETS = {
    "cat1_192": (2, 192, 36, [0, 0, 0, 0, 0, 0, 0, 0, 128], [1, 2, 4, 8, 16, 32, 64, 128, 128],
                 [0, 0, 0, 0, 0, 0, 0, 64, 128], [255, 191], [128, 64], 87),
    "cat1_68":  (4, 68, 42, [0, 0, 0, 0, 8, 8, 8, 8], [1, 2, 4, 8, 8, 16, 32, 64],
                 [0, 0, 0, 4, 0, 0, 0, 64], [71, 11], [64, 4], 26),
    "cat1_45":  (8, 45, 34, [0, 0, 0, 0, 2, 2, 26], [1, 2, 4, 8, 14, 28, 32],
                 [0, 0, 0, 1, 0, 12, 32], [57, 45, 14], [32, 12, 1], 11),
    "cat3_220": (2, 220, 68, [0, 0, 0, 0, 0, 0, 0, 8, 56], [1, 2, 4, 8, 16, 32, 64, 120, 192],
                 [0, 0, 0, 0, 0, 0, 4, 24, 192], [247, 223, 123], [192, 24, 4], 119),
    "cat3_102": (4, 102, 61, [0, 0, 0, 0, 0, 4, 12, 12], [1, 2, 4, 8, 16, 28, 48, 96],
                 [0, 0, 0, 0, 2, 4, 0, 96], [107, 55, 29], [96, 4, 2], 41),
    "cat5_345": (2, 345, 75, [0, 0, 0, 0, 0, 2, 2, 2, 50, 178], [1, 2, 4, 8, 16, 30, 60, 120, 192, 256],
                 [0, 0, 0, 0, 1, 0, 0, 24, 64, 256], [433, 369, 217, 30], [256, 64, 24, 1], 169),
    "cat5_137": (4, 137, 79, [0, 0, 0, 2, 2, 2, 18, 18, 18], [1, 2, 4, 6, 12, 24, 32, 64, 128],
                 [0, 0, 1, 0, 0, 8, 0, 0, 128], [145, 41, 6], [128, 8, 1], 58),
}

PARAMS = os.environ.get("LESS_PARAMS", "cat1_192")
N, K, SEED_LENGTH_BYTES, SIGN_PIVOT_REUSE_LIMIT, Xof, Hash = CATEGORIES[PARAMS[:4]]
(NUM_KEYPAIRS, T, W, TREE_OFFSETS, TREE_NODES_PER_LEVEL, TREE_LEAVES_PER_LEVEL,
 TREE_LEAVES_START_INDICES, TREE_CONSECUTIVE_LEAVES, MAX_PUBLISHED_SEEDS) = TARGETS[PARAMS]
T_LOG2 = len(TREE_OFFSETS) - 1
TREE_SUBROOTS = len(TREE_LEAVES_START_INDICES)

N8 = (N + 7) // 8
K8 = (K + 7) // 8

# BITS_TO_REPRESENT(Q) is 7 for Q = 127
RREF_MAT_PACKEDBYTES = (7 * (N - K) * K + 7) // 8 + (N + 7) // 8
HASH_DIGEST_LENGTH = 2 * SEED_LENGTH_BYTES
SALT_LENGTH_BYTES = HASH_DIGEST_LENGTH
NUM_LEAVES_SEED_TREE = T
NUM_NODES_SEED_TREE = 2 * NUM_LEAVES_SEED_TREE - 1
SEED_TREE_MAX_PUBLISHED_BYTES = MAX_PUBLISHED_SEEDS * SEED_LENGTH_BYTES + 1
# length of the private key seed is doubled to avoid multikey attacks
PRIVATE_KEY_SEED_LENGTH_BYTES = 2 * SEED_LENGTH_BYTES

# number of columns to reuse in vrfy
LESS_VERIFY_PIVOT_REUSE_LIMIT = K

# size of `PublicKey` in bytes, this is what `CRYPTO_PUBLICKEYBYTES` of
# `api.h` resolves to and what names the KAT files
LESS_CRYPTO_PUBLICKEYBYTES = SEED_LENGTH_BYTES + (NUM_KEYPAIRS - 1) * RREF_MAT_PACKEDBYTES
# size of `PrivateKey` in bytes
LESS_CRYPTO_SECRETKEYBYTES = PRIVATE_KEY_SEED_LENGTH_BYTES
# `LESS_SIGNATURE_SIZE(MAX_PUBLISHED_SEEDS)`
LESS_SIGNATURE_MAX_SIZE = HASH_DIGEST_LENGTH * 2 + N8 * W + MAX_PUBLISHED_SEEDS * SEED_LENGTH_BYTES + 1
