//! The parameter set is selected via a cargo feature:
//!
//! | feature     | N   | K   | T   | W  | NUM_KEYPAIRS |
//! |-------------|-----|-----|-----|----|--------------|
//! | `cat1_192`  | 252 | 126 | 192 | 36 | 2            |
//! | `cat1_68`   | 252 | 126 |  68 | 42 | 4            |
//! | `cat1_45`   | 252 | 126 |  45 | 34 | 8            |
//! | `cat3_220`  | 400 | 200 | 220 | 68 | 2            |
//! | `cat3_102`  | 400 | 200 | 102 | 61 | 4            |
//! | `cat5_345`  | 548 | 274 | 345 | 75 | 2            |
//! | `cat5_137`  | 548 | 274 | 137 | 79 | 4            |
//!
//! Exactly one of them must be enabled, a non default set is selected with
//! `--no-default-features --features <name>`.

#![allow(non_upper_case_globals)]

#[cfg(not(any(feature = "cat1", feature = "cat3", feature = "cat5")))]
compile_error!("no parameter set selected: enable one of cat1_192, cat1_68, cat1_45, cat3_220, cat3_102, cat5_345, cat5_137");

/***************************** Category 1 *********************************/
#[cfg(feature = "cat1")]
mod category {
    pub const N: usize = 252;
    pub const K: usize = 126;
    pub const SEED_LENGTH_BYTES: usize = 16;
    /// Ensures probability of non-CT operation is < 2^-64
    pub const SIGN_PIVOT_REUSE_LIMIT: usize = 25;

    /// the XOF, SHAKE128 for this category
    pub type Xof = sha3::Shake128;
    /// the hash, HASH_DIGEST_LENGTH*8 == 256 bits
    pub type Hash = sha3::Sha3_256;
}

/***************************** Category 3 *********************************/
#[cfg(feature = "cat3")]
mod category {
    pub const N: usize = 400;
    pub const K: usize = 200;
    pub const SEED_LENGTH_BYTES: usize = 24;
    /// Ensures probability of non-CT operation is < 2^-64
    pub const SIGN_PIVOT_REUSE_LIMIT: usize = 51;

    /// the XOF, SHAKE256 for this category
    pub type Xof = sha3::Shake256;
    /// the hash, HASH_DIGEST_LENGTH*8 == 384 bits
    pub type Hash = sha3::Sha3_384;
}

/***************************** Category 5 *********************************/
#[cfg(feature = "cat5")]
mod category {
    pub const N: usize = 548;
    pub const K: usize = 274;
    pub const SEED_LENGTH_BYTES: usize = 32;
    /// Ensures probability of non-CT operation is < 2^-64
    pub const SIGN_PIVOT_REUSE_LIMIT: usize = 79;

    /// the XOF, SHAKE256 for this category
    pub type Xof = sha3::Shake256;
    /// the hash, HASH_DIGEST_LENGTH*8 == 512 bits
    pub type Hash = sha3::Sha3_512;
}

/************************* Optimization corner ****************************/
#[cfg(feature = "cat1_192")]
mod target {
    pub const NUM_KEYPAIRS: usize = 2;
    pub const T: usize = 192;
    pub const W: usize = 36;
    pub const T_LOG2: usize = 8;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 0, 0, 0, 128];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 16, 32, 64, 128, 128];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 0, 0, 64, 128];
    pub const TREE_SUBROOTS: usize = 2;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [255, 191];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [128, 64];
    pub const MAX_PUBLISHED_SEEDS: usize = 87;
}

#[cfg(feature = "cat1_68")]
mod target {
    pub const NUM_KEYPAIRS: usize = 4;
    pub const T: usize = 68;
    pub const W: usize = 42;
    pub const T_LOG2: usize = 7;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 8, 8, 8, 8];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 8, 16, 32, 64];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 4, 0, 0, 0, 64];
    pub const TREE_SUBROOTS: usize = 2;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [71, 11];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [64, 4];
    pub const MAX_PUBLISHED_SEEDS: usize = 26;
}

#[cfg(feature = "cat1_45")]
mod target {
    pub const NUM_KEYPAIRS: usize = 8;
    pub const T: usize = 45;
    pub const W: usize = 34;
    pub const T_LOG2: usize = 6;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 2, 2, 26];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 14, 28, 32];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 1, 0, 12, 32];
    pub const TREE_SUBROOTS: usize = 3;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [57, 45, 14];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [32, 12, 1];
    pub const MAX_PUBLISHED_SEEDS: usize = 11;
}

#[cfg(feature = "cat3_220")]
mod target {
    pub const NUM_KEYPAIRS: usize = 2;
    pub const T: usize = 220;
    pub const W: usize = 68;
    pub const T_LOG2: usize = 8;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 0, 0, 8, 56];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 16, 32, 64, 120, 192];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 0, 4, 24, 192];
    pub const TREE_SUBROOTS: usize = 3;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [247, 223, 123];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [192, 24, 4];
    pub const MAX_PUBLISHED_SEEDS: usize = 119;
}

#[cfg(feature = "cat3_102")]
mod target {
    pub const NUM_KEYPAIRS: usize = 4;
    pub const T: usize = 102;
    pub const W: usize = 61;
    pub const T_LOG2: usize = 7;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 4, 12, 12];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 16, 28, 48, 96];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 2, 4, 0, 96];
    pub const TREE_SUBROOTS: usize = 3;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [107, 55, 29];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [96, 4, 2];
    pub const MAX_PUBLISHED_SEEDS: usize = 41;
}

#[cfg(feature = "cat5_345")]
mod target {
    pub const NUM_KEYPAIRS: usize = 2;
    pub const T: usize = 345;
    pub const W: usize = 75;
    pub const T_LOG2: usize = 9;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 0, 2, 2, 2, 50, 178];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 8, 16, 30, 60, 120, 192, 256];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 0, 0, 1, 0, 0, 24, 64, 256];
    pub const TREE_SUBROOTS: usize = 4;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [433, 369, 217, 30];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [256, 64, 24, 1];
    pub const MAX_PUBLISHED_SEEDS: usize = 169;
}

#[cfg(feature = "cat5_137")]
mod target {
    pub const NUM_KEYPAIRS: usize = 4;
    pub const T: usize = 137;
    pub const W: usize = 79;
    pub const T_LOG2: usize = 8;
    pub const TREE_OFFSETS: [usize; T_LOG2 + 1] = [0, 0, 0, 2, 2, 2, 18, 18, 18];
    pub const TREE_NODES_PER_LEVEL: [usize; T_LOG2 + 1] = [1, 2, 4, 6, 12, 24, 32, 64, 128];
    pub const TREE_LEAVES_PER_LEVEL: [usize; T_LOG2 + 1] = [0, 0, 1, 0, 0, 8, 0, 0, 128];
    pub const TREE_SUBROOTS: usize = 3;
    pub const TREE_LEAVES_START_INDICES: [usize; TREE_SUBROOTS] = [145, 41, 6];
    pub const TREE_CONSECUTIVE_LEAVES: [usize; TREE_SUBROOTS] = [128, 8, 1];
    pub const MAX_PUBLISHED_SEEDS: usize = 58;
}

pub use category::{Hash, Xof, K, N, SEED_LENGTH_BYTES, SIGN_PIVOT_REUSE_LIMIT};
pub use target::{
    MAX_PUBLISHED_SEEDS, NUM_KEYPAIRS, T, T_LOG2, W,
    TREE_CONSECUTIVE_LEAVES, TREE_LEAVES_PER_LEVEL, TREE_LEAVES_START_INDICES,
    TREE_NODES_PER_LEVEL, TREE_OFFSETS, TREE_SUBROOTS,
};

/// the squeezing state of the XOF
pub type XofState = <Xof as sha3::digest::ExtendableOutput>::Reader;

pub const N8: usize = (N+7)/8;
pub const K8: usize = (K+7)/8;

/// rounds x to the next multiple of n
const fn next_multiple(x: usize, n: usize) -> usize {
    ((x + (n - 1)) / n) * n
}

/// In case of the optimized implementation, we need that all vectors
/// are of a length, which is a multiple of 32
#[cfg(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon"))]
pub const N_pad: usize = next_multiple(N, 32);
#[cfg(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon"))]
pub const K_pad: usize = next_multiple(K, 32);
#[cfg(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon"))]
pub const N_K_pad: usize = next_multiple(N - K, 32);

/// in case of the reference implementation, we do not need this behaviour.
#[cfg(not(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon")))]
pub const N_pad: usize = N;
#[cfg(not(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon")))]
pub const K_pad: usize = K;
#[cfg(not(any(target_feature = "avx2", target_feature = "avx512f", target_feature = "neon")))]
pub const N_K_pad: usize = N - K;


/// BITS_TO_REPRESENT(Q) is 7 for Q = 127
pub const RREF_MAT_PACKEDBYTES: usize = (7*(N-K)*K + 7)/8 + (N + 7)/8;
pub const HASH_DIGEST_LENGTH: usize = 2*SEED_LENGTH_BYTES;
pub const SALT_LENGTH_BYTES: usize = HASH_DIGEST_LENGTH;
pub const NUM_LEAVES_SEED_TREE: usize = T;
pub const NUM_NODES_SEED_TREE: usize = (2*NUM_LEAVES_SEED_TREE) - 1;
pub const SEED_TREE_MAX_PUBLISHED_BYTES: usize = MAX_PUBLISHED_SEEDS*SEED_LENGTH_BYTES + 1;
/// length of the private key seed is doubled to avoid multikey attacks
pub const PRIVATE_KEY_SEED_LENGTH_BYTES: usize = 2*SEED_LENGTH_BYTES;

/// number of columns to reuse in vrfy
pub const LESS_VERIFY_PIVOT_REUSE_LIMIT: usize = K;

/// size of `PublicKey` in bytes, this is what `CRYPTO_PUBLICKEYBYTES` of
/// `api.h` resolves to and what names the KAT files
pub const LESS_CRYPTO_PUBLICKEYBYTES: usize = SEED_LENGTH_BYTES + (NUM_KEYPAIRS-1)*RREF_MAT_PACKEDBYTES;
/// size of `PrivateKey` in bytes
pub const LESS_CRYPTO_SECRETKEYBYTES: usize = PRIVATE_KEY_SEED_LENGTH_BYTES;
/// `LESS_SIGNATURE_SIZE(MAX_PUBLISHED_SEEDS)`
pub const LESS_SIGNATURE_MAX_SIZE: usize =
    HASH_DIGEST_LENGTH*2 + N8*W + MAX_PUBLISHED_SEEDS*SEED_LENGTH_BYTES + 1;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_constants() {
        assert_eq!(N, 2*K);
        assert_eq!(N8, (N+7)/8);
        assert_eq!(K8, (K+7)/8);
        assert_eq!(HASH_DIGEST_LENGTH, 2*SEED_LENGTH_BYTES);
        assert_eq!(SALT_LENGTH_BYTES, HASH_DIGEST_LENGTH);
        assert_eq!(PRIVATE_KEY_SEED_LENGTH_BYTES, 2*SEED_LENGTH_BYTES);
        assert_eq!(NUM_NODES_SEED_TREE, 2*T - 1);
        assert_eq!(N_pad % 32, 0);
        assert_eq!(K_pad % 32, 0);
        assert_eq!(N_K_pad % 32, 0);
        assert!(N_pad >= N && N_pad - N < 32);
        assert!(K_pad >= K && K_pad - K < 32);
        assert!(N_K_pad >= N - K && N_K_pad - (N - K) < 32);
        assert!(W < T);
    }

    /// the digest of `Hash` must be exactly HASH_DIGEST_LENGTH bytes
    #[test]
    fn hash_output_length() {
        use sha3::Digest;
        assert_eq!(<Hash as Digest>::output_size(), HASH_DIGEST_LENGTH);
    }

    /// every leaf is accounted for exactly once, and the number of nodes per
    /// level sums up to the size of the linearized tree
    #[test]
    fn tree_constants() {
        let mut leaves = 0;
        for i in 0..TREE_SUBROOTS {
            leaves += TREE_CONSECUTIVE_LEAVES[i];
        }
        assert_eq!(leaves, T);

        let mut nodes = 0;
        for i in 0..=T_LOG2 {
            nodes += TREE_NODES_PER_LEVEL[i];
        }
        assert_eq!(nodes, NUM_NODES_SEED_TREE);

        for i in 0..TREE_SUBROOTS {
            assert!(TREE_LEAVES_START_INDICES[i] + TREE_CONSECUTIVE_LEAVES[i] <= NUM_NODES_SEED_TREE);
        }
        for i in 0..=T_LOG2 {
            assert!(TREE_LEAVES_PER_LEVEL[i] <= TREE_NODES_PER_LEVEL[i]);
        }
    }

    /// T_LOG2 must be LOG2(T)
    #[test]
    fn t_log2() {
        let bits = |n: usize| -> usize { usize::BITS as usize - n.leading_zeros() as usize };
        let log2 = if bits(T) > bits(T-1) { bits(T-1) } else { bits(T) };
        assert_eq!(T_LOG2, log2);
    }

    #[cfg(feature = "cat1_192")]
    #[test]
    fn cat1_192() {
        assert_eq!((N, K, T, W, NUM_KEYPAIRS), (252, 126, 192, 36, 2));
        assert_eq!(SEED_LENGTH_BYTES, 16);
        assert_eq!(HASH_DIGEST_LENGTH, 32);
        assert_eq!(RREF_MAT_PACKEDBYTES, 13924);
        assert_eq!(NUM_NODES_SEED_TREE, 383);
        assert_eq!(SEED_TREE_MAX_PUBLISHED_BYTES, 1393);
        assert_eq!(SIGN_PIVOT_REUSE_LIMIT, 25);
        assert_eq!(N8, 32);
    }

    #[cfg(feature = "cat3_220")]
    #[test]
    fn cat3_220() {
        assert_eq!((N, K, T, W, NUM_KEYPAIRS), (400, 200, 220, 68, 2));
        assert_eq!(SEED_LENGTH_BYTES, 24);
        assert_eq!(HASH_DIGEST_LENGTH, 48);
        assert_eq!(NUM_NODES_SEED_TREE, 439);
        assert_eq!(N8, 50);
    }

    #[cfg(feature = "cat5_345")]
    #[test]
    fn cat5_345() {
        assert_eq!((N, K, T, W, NUM_KEYPAIRS), (548, 274, 345, 75, 2));
        assert_eq!(SEED_LENGTH_BYTES, 32);
        assert_eq!(HASH_DIGEST_LENGTH, 64);
        assert_eq!(NUM_NODES_SEED_TREE, 689);
        assert_eq!(N8, 69);
    }
}
