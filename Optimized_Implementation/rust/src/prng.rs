#![allow(unused_imports)]

use sha3::{
    Digest, Sha3_256, Sha3_512, Shake128, Shake256, Shake128Core, Shake128ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
    digest::core_api::XofReaderCoreWrapper,
};

use std::sync::{LazyLock, Mutex};
use crate::fq::Fq;

use crate::error::LESSError;
use sha3::digest::core_api::CoreWrapper;

static PLATFORM_CSPRNG_STATE: LazyLock<Mutex<CoreWrapper<Shake128Core>>> = LazyLock::new(||
    Mutex::new(Shake128::default())
);

pub fn init_randombytes(buf: &[u8]) {
    let mut hasher = PLATFORM_CSPRNG_STATE.lock().unwrap();
    hasher.update(buf);

}
pub fn randombytes(buf: &mut [u8]) {
    let hasher = PLATFORM_CSPRNG_STATE.lock().unwrap();
    let mut reader = hasher.clone().finalize_xof();
    reader.read(buf);
}


pub fn csprng_initialize<R>(state: &mut R, buf: &[u8]) where R: ExtendableOutput {
    state.update(buf);
    // TODO state.finalize_xof();
}
pub fn csprng_randombytes<R>(buf: &mut [u8], state: &mut  R) where R: XofReader {
    state.read(buf);
}
fn rand_range_impl<const MIN_V: u8, R>(
    buffer: &mut [Fq],
    state: &mut R,
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
    buffer: &mut [Fq],
    state: &mut R,
)
where 
    R: XofReader
{
    rand_range_impl::<0, R>(buffer, state);
}



pub fn merge_exchange<R>(
    buffer: &mut [u16],
    state: &mut R,
)
where
    R: XofReader
{
    let n = buffer.len() as u32;

    let mut rand: u64 = 0;
    let mut ctr: u32 = 0;

    // Find the largest power of 2 less than n
    let mut t: u32 = 1;
    while t < n - t {
        t += t;
    }

    for p_step in (0..).map(|i| t >> i).take_while(|&x| x > 0) {
        let mut q = t << 1u32;
        let mut r: u32 = 0;
        let mut d: u32 = p_step;
        loop {
            q >>= 1;
            for i in 0..(n - d) {
                if (i & p_step) == r {
                    // Refill random bits if exhausted
                    if ctr == 0 {
                        let mut buf: [u8; 8] = [0u8; 8];
                        state.read(&mut buf);
                        rand = u64::from_le_bytes(buf);
                        ctr = 64;
                    }

                    // Constant-time mask: 0xFFFF if (rand & 1) == 1, else 0x0000
                    let mask: u16 = if (rand & 1) == 1 { 0xFFFF } else { 0 };

                    // Conditional swap using masked XOR
                    let a = buffer[i as usize];
                    let b = buffer[(i + d) as usize];
                    buffer[i as usize] = a ^ (mask & (a ^ b));
                    buffer[(i + d) as usize] = b ^ (mask & (a ^ b));

                    // Drain random bit and counter
                    rand >>= 1;
                    ctr -= 1;
                }
            }
            d = q - p_step;
            r = p_step;
            if q == p_step {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init_randombytes() {
        let seed = [0u8; 32];
        init_randombytes(&seed);

        let mut out_buf = [0u8; 32];
        randombytes(&mut out_buf);
        // TODO not really correct
        for i in 0..32 {
            assert_ne!(out_buf[i], 0u8);
        }
    }

    #[test]
    fn test_rand_range_impl() {
        let mut hasher = Shake128::default();
        hasher.update(b"123");
        let mut reader = hasher.finalize_xof();
        let mut buffer = [Fq(0); 32];
        rand_range_impl::<1, XofReaderCoreWrapper<Shake128ReaderCore>>(&mut buffer, &mut reader);
        for i in 0..32 {
            assert_ne!(buffer[i], Fq(0));
        }
    }

    #[test]
    fn test_merge_exchange() {
        let mut hasher = Shake128::default();
        hasher.update(b"123");
        let mut reader = hasher.finalize_xof();
        let mut buffer: [u16; 128] = (0..128).collect::<Vec<_>>().try_into().unwrap();
        merge_exchange(&mut buffer, &mut reader);

        let mut cnt = 0u32;
        for i in 0..128 {
            cnt += (buffer[i] == i as u16) as u32;
        }

        assert!(cnt < 5);
    }
}
