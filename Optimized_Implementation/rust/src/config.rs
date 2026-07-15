
pub const K: usize = 126;
pub const N: usize = 252;
pub const T: usize = 192;
pub const T_LOG2: usize = 8;
pub const W: usize = 32;
pub const N8: usize = (N+7)/8;

pub const N_pad: usize = 256;
pub const K_pad: usize = 128;
pub const N_K_pad: usize = 128;


pub const SEED_LENGTH_BYTES: usize  = 32;
pub const RREF_MAT_PACKEDBYTES: usize = 1000;
pub const NUM_KEYPAIRS: usize = 2;
pub const HASH_DIGEST_LENGTH: usize = 2*SEED_LENGTH_BYTES;
pub const SALT_LENGTH_BYTES: usize = HASH_DIGEST_LENGTH;
pub const NUM_LEAVES_SEED_TREE: usize = T;
pub const NUM_NODES_SEED_TREE: usize = ((2*NUM_LEAVES_SEED_TREE) - 1);
pub const SEED_TREE_MAX_PUBLISHED_BYTES: usize = MAX_PUBLISHED_SEEDS*SEED_LENGTH_BYTES + 1;
/// length of the private key seed is doubled to avoid multikey attacks
pub const PRIVATE_KEY_SEED_LENGTH_BYTES: usize = (2*SEED_LENGTH_BYTES);

pub const SIGN_PIVOT_REUSE_LIMIT: usize = 51;

pub const TREE_OFFSETS: [usize; 9] = [0, 0, 0, 0, 0, 0, 0, 0, 128];
pub const TREE_NODES_PER_LEVEL: [usize; 9] = [1, 2, 4, 8, 16, 32, 64, 128, 128];
pub const TREE_LEAVES_PER_LEVEL: [usize; 9] = [0, 0, 0, 0, 0, 0, 0, 64, 128];
pub const TREE_SUBROOTS: usize = 2;
pub const TREE_LEAVES_START_INDICES: [usize; 2] = [255, 191];
pub const TREE_CONSECUTIVE_LEAVES: [usize; 2] = [128, 64];
pub const MAX_PUBLISHED_SEEDS: usize = 87;