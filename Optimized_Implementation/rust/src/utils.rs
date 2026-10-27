use crate::config::{
    T, W, Xof,
    HASH_DIGEST_LENGTH,
    NUM_KEYPAIRS,
};

use crate::prng::{
    csprng_randombytes,
    csprng_initialize
};

/// Compare two arrays for equality in constant time.
/// \param a[in]: pointer to the first byte array
/// \param b[in]: pointer to the second byte array
/// \returns 0 if the byte arrays are equal, 1 otherwise
pub fn verify(a: &[u8], b: &[u8]) -> u32 {
    let mut r = 0u8;

    for i in 0..a.len() {
        r |= a[i] ^ b[i];
    }

    (((r as u64).wrapping_neg()) >> 63) as u32
}

/// Expands a digest into a fixed weight string with elements in Z_{NUM_KEYPAIRS}.
/// \param fixed_weight_string[out]: array of length T, with hamming weight W,
///     where values > 0 are between [1, NUM_KEYPAIRS-1)
pub fn sample_challenge(
    fixed_weight_string: &mut [u8; T],
    digest: &[u8; HASH_DIGEST_LENGTH],
) {
    // Initialize CSPRNG from the digest
    let mut shake_state = csprng_initialize::<Xof>(digest);

    let max_keypair_idx: u32 = NUM_KEYPAIRS as u32 - 1;
    let keypair_bits: u32 = (max_keypair_idx.ilog2() + 1) as u32;

    let t_minus_1: u32 = T as u32 - 1;
    let position_bits: u32 = (t_minus_1.ilog2() + 1) as u32;

    // Clear the first T-W entries to zero
    for entry in fixed_weight_string[..T - W].iter_mut() {
        *entry = 0;
    }

    // NOTE: `rnd_buf` and `c` are shared between the two loops below, the
    // shuffle keeps consuming the bits the value loop left over.
    let mut rnd_buf: u64 = 0;
    let mut c: u32 = 0;

    if NUM_KEYPAIRS != 2 {
        // Rejection sampling for values in [1, NUM_KEYPAIRS-1)
        let keypair_mask: u64 = ((1u64 << keypair_bits) - 1) & !(0xFF << keypair_bits);
        let limit: u8 = NUM_KEYPAIRS as u8 - 1;

        for i in (T - W)..T {
            let value: u8;
            loop {
                if c == 0 {
                    let buf: [u8; 8] = csprng_randombytes(&mut shake_state, 8).try_into().unwrap();
                    rnd_buf = u64::from_le_bytes(buf);
                    c = 64 / keypair_bits;
                }

                let v: u8 = (rnd_buf & keypair_mask) as u8;
                rnd_buf >>= keypair_bits;
                c -= 1;

                if v < limit {
                    value = v;
                    break;
                }
            }
            fixed_weight_string[i] = value + 1;
        }
    } else {
        // NUM_KEYPAIRS == 2: the last W values are 1
        for entry in fixed_weight_string[T - W..].iter_mut() {
            *entry = 1;
        }
    }

    // Fisher-Yates-like shuffle into the first T-W zero positions
    let position_mask: u64 = ((1u64 << position_bits) - 1) & !(0xFF << position_bits);

    for p in (T - W)..T {
        let mut pos: usize;
        loop {
            if c == 0 {
                let buf: [u8; 8] = csprng_randombytes(&mut shake_state, 8).try_into().unwrap();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_challenge_weight() {
        let digest = [0x5au8; HASH_DIGEST_LENGTH];
        let mut fws = [0u8; T];
        sample_challenge(&mut fws, &digest);

        let weight = fws.iter().filter(|&&v| v != 0).count();
        assert_eq!(weight, W, "the challenge must have hamming weight W");
    }

    #[test]
    fn sample_challenge_range() {
        let digest = [0x11u8; HASH_DIGEST_LENGTH];
        let mut fws = [0u8; T];
        sample_challenge(&mut fws, &digest);

        for i in 0..T {
            assert!((fws[i] as usize) < NUM_KEYPAIRS);
        }
    }

    #[test]
    fn sample_challenge_deterministic() {
        let digest = [0x7fu8; HASH_DIGEST_LENGTH];
        let mut a = [0u8; T];
        let mut b = [0u8; T];
        sample_challenge(&mut a, &digest);
        sample_challenge(&mut b, &digest);
        assert_eq!(a, b);
    }

    #[test]
    fn sample_challenge_depends_on_digest() {
        let mut a = [0u8; T];
        let mut b = [0u8; T];
        sample_challenge(&mut a, &[0x00u8; HASH_DIGEST_LENGTH]);
        sample_challenge(&mut b, &[0x01u8; HASH_DIGEST_LENGTH]);
        assert_ne!(a, b);
    }

    /// the non zero entries must not all sit in the last W positions, the
    /// shuffle has to spread them over the full length
    #[test]
    fn sample_challenge_is_shuffled() {
        let mut spread = false;
        for d in 0..8u8 {
            let mut fws = [0u8; T];
            sample_challenge(&mut fws, &[d; HASH_DIGEST_LENGTH]);
            if fws[..T - W].iter().any(|&v| v != 0) {
                spread = true;
                break;
            }
        }
        assert!(spread);
    }
}