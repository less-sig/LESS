#![allow(unused_imports)]

use sha3::{
    Digest, Sha3_256, Sha3_512, Shake128, Shake256, Shake128ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
    digest::core_api::XofReaderCoreWrapper,
};

use crate::fq::Fq;

//use crate::KyberError;
//use rand_core::*;

/// Fills buffer x with len bytes, RNG must satisfy the
/// RngCore trait and CryptoRng marker trait requirements
//pub fn randombytes<R>(x: &mut [u8], len: usize, rng: &mut R) -> Result<(), KyberError>
//where
//    R: RngCore + CryptoRng,
//{
//    match rng.try_fill_bytes(&mut x[..len]) {
//        Ok(_) => Ok(()),
//        Err(_) => Err(KyberError::RandomBytesGeneration),
//    }
//}

fn rand_range_impl<const MIN_V: u8, R>(
    mut state: R,
    buffer: &mut [Fq],
)
where 
    R: XofReader
{
    let span = 127 - MIN_V;
    let req_bits: usize = 7;
    let el_mask: u8 = (1u8 << req_bits) - 1;

    let mut count: usize = 0;
    let mut buf = [0u8; 8];

    while count < buffer.len() {
        state.read(&mut buf);
        let mut word: u64 = u64::from_ne_bytes(buf);
        let iters = (size_of::<u64>() * 8) / req_bits;
        for _ in 0..iters {
            let rnd_value = word as u8 & el_mask;
            if rnd_value < span {
                buffer[count] = Fq(MIN_V + rnd_value);
                count += 1;
                if count >= buffer.len() {
                    return;
                }
            }
            word >>= req_bits;
        }
    }
}

pub fn rand_range_q_state_elements<R>(
    state: R,
    buffer: &mut [Fq],
)
where 
    R: XofReader
{
    rand_range_impl::<0, R>(state, buffer);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple() {
        let mut hasher = Shake128::default();
        hasher.update(b"123");
        let reader = hasher.finalize_xof();
        let mut buffer = [Fq(0); 32];
        rand_range_impl::<1, XofReaderCoreWrapper<Shake128ReaderCore>>(reader, &mut buffer);
        for i in 0..32 {
            assert_ne!(buffer[i], Fq(0));
        }
    }
}
