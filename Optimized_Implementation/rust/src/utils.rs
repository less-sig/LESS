use sha3::{
    Shake128
};

use crate::config::{
    T, W,
    HASH_DIGEST_LENGTH,
    NUM_KEYPAIRS,
};

use crate::prng::{
    csprng_randombytes,
    csprng_initialize
};
/// Expands a digest into a fixed weight string with elements in Z_{NUM_KEYPAIRS}.
/// \param fixed_weight_string[out]: array of length T, with hamming weight W,
///     where values > 0 are between [1, NUM_KEYPAIRS-1)
pub fn sample_challenge(
    fixed_weight_string: &mut [u8; T],
    digest: &[u8; HASH_DIGEST_LENGTH],
) {
    // Initialize CSPRNG from the digest
    let mut shake_state = Shake128::default();
    csprng_initialize(&mut shake_state, digest);

    let max_keypair_idx: u32 = NUM_KEYPAIRS as u32 - 1;
    let keypair_bits: u32 = (max_keypair_idx.ilog2() + 1) as u32;

    let t_minus_1: u32 = T as u32 - 1;
    let position_bits: u32 = (t_minus_1.ilog2() + 1) as u32;

    // Clear the first T-W entries to zero
    for entry in fixed_weight_string[..T - W].iter_mut() {
        *entry = 0;
    }

    if NUM_KEYPAIRS != 2 {
        // Rejection sampling for values in [1, NUM_KEYPAIRS-1)
        let keypair_mask: u64 = ((1u64 << keypair_bits) - 1) & !(0xFF << keypair_bits);
        let limit: u8 = NUM_KEYPAIRS as u8 - 1;

        let mut rnd_buf: u64 = 0;
        let mut c: u32 = 64 / keypair_bits;

        for i in (T - W)..T {
            loop {
                if c == 0 {
                    let mut buf: [u8; 8] = [0u8; 8];
                    // TODO csprng_randombytes(&mut buf, &mut shake_state);
                    rnd_buf = u64::from_le_bytes(buf);
                    c = 64 / keypair_bits;
                }

                let value: u8 = (rnd_buf & keypair_mask) as u8;
                rnd_buf >>= keypair_bits;
                c -= 1;

                if value < limit {
                    break;
                }
                fixed_weight_string[i] = value + 1;
            }
        }
    } else {
        // NUM_KEYPAIRS == 2: all values are 1
        for entry in fixed_weight_string.iter_mut() {
            *entry = 1;
        }
    }

    // Fisher-Yates-like shuffle into the first T-W zero positions
    let mut rnd_buf: u64 = 0;
    let mut c: u32 = 0;
    let position_mask: u64 = ((1u64 << position_bits) - 1) & !(0xFF << position_bits);

    for p in (T - W)..T {
        let mut pos: usize;
        loop {
            if c == 0 {
                let mut buf: [u8; 8] = [0u8; 8];
                // TODO csprng_randombytes(&mut buf, &mut shake_state);
                rnd_buf = u64::from_le_bytes(buf);
                c = 64 / position_bits;
            }

            pos = (rnd_buf & position_mask) as usize;
            rnd_buf >>= position_bits;
            c -= 1;

            if pos <= p {
                break;
            }
        }

        // Swap
        let tmp = fixed_weight_string[p];
        fixed_weight_string[p] = fixed_weight_string[pos];
        fixed_weight_string[pos] = tmp;
    }
}