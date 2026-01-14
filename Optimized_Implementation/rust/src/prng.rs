#![allow(unused_imports)]



//use crate::KyberError;
//use rand_core::*;
//
///// Fills buffer x with len bytes, RNG must satisfy the
///// RngCore trait and CryptoRng marker trait requirements
//pub fn randombytes<R>(x: &mut [u8], len: usize, rng: &mut R) -> Result<(), KyberError>
//where
//    R: RngCore + CryptoRng,
//{
//    match rng.try_fill_bytes(&mut x[..len]) {
//        Ok(_) => Ok(()),
//        Err(_) => Err(KyberError::RandomBytesGeneration),
//    }
//}
//// how to use it correclty: https://github.com/RustCrypto/KEMs/blob/master/ml-kem/src/crypto.rs
//fn rand_range_impl<EL, R>(
//    mut fill_word: R,
//    buffer: &mut [EL],
//    min_value: EL,
//    max_value: EL,
//)
//where
//    EL: Copy
//        + From<u8>
//        + core::ops::Add<Output = EL>
//        + core::ops::Sub<Output = EL>
//        + PartialOrd
//        + Into<u64>,
//    R: FnMut(&mut u64),
//{
//    let span: u64 = (max_value - min_value).into();
//    let req_bits: usize = bits_to_represent(span);
//    let el_mask: u64 = (1u64 << req_bits) - 1;
//
//    let mut count: usize = 0;
//    let mut word: u64 = 0;
//
//    while count < buffer.len() {
//        fill_word(&mut word);
//
//        let iters = (size_of::<u64>() * 8) / req_bits;
//        for _ in 0..iters {
//            let rnd_value = word & el_mask;
//            if rnd_value <= span {
//                buffer[count] = EL::from(0u8) + min_value + EL::from(rnd_value as u8);
//                count += 1;
//                if count >= buffer.len() {
//                    return;
//                }
//            }
//            word >>= req_bits;
//        }
//    }
//}
//pub fn rand_with_state<EL>(
//    shake_state: &mut ShakeState,
//    buffer: &mut [EL],
//    min_value: EL,
//    max_value: EL,
//)
//where
//    EL: Copy
//        + From<u8>
//        + core::ops::Add<Output = EL>
//        + core::ops::Sub<Output = EL>
//        + PartialOrd
//        + Into<u64>,
//{
//    rand_range_impl(
//        |word| {
//            csprng_randombytes(
//                word as *mut u64 as *mut u8,
//                core::mem::size_of::<u64>(),
//                shake_state,
//            );
//        },
//        buffer,
//        min_value,
//        max_value,
//    );
//}
//
//
//#[cfg(test)]
//mod tests {
//    use super::*;
//
//    #[test]
//    fn it_works() {
//
//    }
//}
