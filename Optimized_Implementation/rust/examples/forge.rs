//! Brute forces false signatures and checks that `less_verify` rejects every
//! one of them.
//!
//! Each strategy mutates a valid signature, or builds one out of random bytes,
//! and the result is tallied as rejected, accepted (a forgery) or panicked (a
//! denial of service). The example exits non zero if any forgery is accepted.
//!
//! Use `--release`, the debug build is orders of magnitude slower.

use less::config::*;
use less::less::*;
use less::prng::init_randombytes;
use std::panic::{self, AssertUnwindSafe};

/// deterministic xorshift, so that a failing trial can be reproduced
fn xorshift(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// a uniform value below `n`
fn below(state: &mut u64, n: usize) -> usize {
    (xorshift(state) as usize) % n
}

/// Number of bytes a serialized signature really carries, the way
/// `crypto_sign_open` splits it: only `leaves*SEED_LENGTH_BYTES` of
/// `seed_storage` are transmitted, the rest of that field is never sent and
/// mutating it cannot change the outcome.
fn sig_bytes(leaves: usize) -> usize {
    HASH_DIGEST_LENGTH * 2 + N8 * W + leaves * SEED_LENGTH_BYTES
}

/// addresses one byte of a signature as if the four fields were contiguous
fn sig_byte(sig: &mut Signature, idx: usize) -> &mut u8 {
    let mut i = idx;
    if i < HASH_DIGEST_LENGTH {
        return &mut sig.digest[i];
    }
    i -= HASH_DIGEST_LENGTH;
    if i < HASH_DIGEST_LENGTH {
        return &mut sig.salt[i];
    }
    i -= HASH_DIGEST_LENGTH;
    if i < N8 * W {
        return &mut sig.cf_monom_actions[i / N8][i % N8];
    }
    i -= N8 * W;
    &mut sig.seed_storage[i]
}

/// a copy of a signature, so that a single field can be tampered with
fn copy_sig(sig: &Signature) -> Signature {
    Signature {
        digest: sig.digest,
        salt: sig.salt,
        cf_monom_actions: sig.cf_monom_actions,
        seed_storage: sig.seed_storage,
    }
}

/// a copy of a public key, so that a single byte can be tampered with
fn copy_pk(pk: &PublicKey) -> PublicKey {
    PublicKey { g0_seed: pk.g0_seed, sf_g: pk.sf_g }
}

/// keygen and sign once, from a fixed seed
fn keypair_and_signature(msg: &[u8]) -> (PublicKey, Signature, usize) {
    let seed: [u8; 48] = core::array::from_fn(|i| i as u8);
    init_randombytes(&seed);

    let mut sk = PrivateKey { compressed_sk: [0u8; PRIVATE_KEY_SEED_LENGTH_BYTES] };
    let mut pk = PublicKey {
        g0_seed: [0u8; SEED_LENGTH_BYTES],
        sf_g: [[0u8; RREF_MAT_PACKEDBYTES]; NUM_KEYPAIRS - 1],
    };
    less_keygen(&mut sk, &mut pk);

    let mut sig = Signature {
        digest: [0u8; HASH_DIGEST_LENGTH],
        salt: [0u8; HASH_DIGEST_LENGTH],
        cf_monom_actions: [[0u8; N8]; W],
        seed_storage: [0u8; SEED_TREE_MAX_PUBLISHED_BYTES],
    };
    let leaves = less_sign(&mut sig, msg, &sk, &pk);

    (pk, sig, leaves)
}

#[derive(Default)]
struct Tally {
    rejected: u64,
    accepted: u64,
    panicked: u64,
    first_panic: Option<String>,
}

const STRATEGIES: &[&str] = &[
    "1 bit flipped in digest",
    "1 bit flipped in salt",
    "1 bit flipped in cf_monom_actions",
    "1 bit flipped in seed_storage",
    "1 bit flipped anywhere",
    "8 bits flipped anywhere",
    "digest fully random",
    "salt fully random",
    "cf_monom_actions fully random",
    "seed_storage fully random",
    "whole signature random",
    "1 bit flipped in the message",
    "1 bit flipped in the public key",
];

fn main() {
    let trials: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(256);

    let msg = b"LESS rust port forgery search";
    let (pk, sig, leaves) = keypair_and_signature(msg);

    assert!(leaves > 0 && leaves <= MAX_PUBLISHED_SEEDS,
            "{leaves} published leaves is out of range");
    assert!(less_verify(&pk, msg, &sig), "the honest signature was rejected");

    let sig_bytes = sig_bytes(leaves);
    println!("N={N} K={K} T={T} W={W}, signature {sig_bytes} bytes, {leaves} published leaves");
    println!("{trials} trials per strategy\n");

    /* the verifier's own panic messages are tallied, not printed */
    let captured: std::sync::Arc<std::sync::Mutex<Option<String>>> =
        std::sync::Arc::new(std::sync::Mutex::new(None));
    {
        let captured = std::sync::Arc::clone(&captured);
        panic::set_hook(Box::new(move |info| {
            *captured.lock().unwrap() = Some(format!("{info}"));
        }));
    }

    let mut state: u64 = 0x243f_6a88_85a3_08d3;
    let mut flagged: Vec<(&str, Tally)> = Vec::new();

    for (s, name) in STRATEGIES.iter().enumerate() {
        let mut t = Tally::default();

        for _ in 0..trials {
            let mut bad = copy_sig(&sig);
            let mut bad_pk = copy_pk(&pk);
            let mut bad_msg = *msg;

            match s {
                0 => {
                    let i = below(&mut state, HASH_DIGEST_LENGTH);
                    bad.digest[i] ^= 1 << below(&mut state, 8);
                }
                1 => {
                    let i = below(&mut state, HASH_DIGEST_LENGTH);
                    bad.salt[i] ^= 1 << below(&mut state, 8);
                }
                2 => {
                    let i = below(&mut state, W);
                    let j = below(&mut state, N8);
                    bad.cf_monom_actions[i][j] ^= 1 << below(&mut state, 8);
                }
                3 => {
                    let i = below(&mut state, leaves * SEED_LENGTH_BYTES);
                    bad.seed_storage[i] ^= 1 << below(&mut state, 8);
                }
                4 => {
                    let i = below(&mut state, sig_bytes);
                    *sig_byte(&mut bad, i) ^= 1 << below(&mut state, 8);
                }
                5 => {
                    for _ in 0..8 {
                        let i = below(&mut state, sig_bytes);
                        *sig_byte(&mut bad, i) ^= 1 << below(&mut state, 8);
                    }
                }
                6 => {
                    for v in bad.digest.iter_mut() {
                        *v = xorshift(&mut state) as u8;
                    }
                }
                7 => {
                    for v in bad.salt.iter_mut() {
                        *v = xorshift(&mut state) as u8;
                    }
                }
                8 => {
                    for row in bad.cf_monom_actions.iter_mut() {
                        for v in row.iter_mut() {
                            *v = xorshift(&mut state) as u8;
                        }
                    }
                }
                9 => {
                    for v in bad.seed_storage[..leaves * SEED_LENGTH_BYTES].iter_mut() {
                        *v = xorshift(&mut state) as u8;
                    }
                }
                10 => {
                    for i in 0..sig_bytes {
                        *sig_byte(&mut bad, i) = xorshift(&mut state) as u8;
                    }
                }
                11 => {
                    let i = below(&mut state, bad_msg.len());
                    bad_msg[i] ^= 1 << below(&mut state, 8);
                }
                _ => {
                    if below(&mut state, 2) == 0 {
                        let i = below(&mut state, SEED_LENGTH_BYTES);
                        bad_pk.g0_seed[i] ^= 1 << below(&mut state, 8);
                    } else {
                        let i = below(&mut state, NUM_KEYPAIRS - 1);
                        let j = below(&mut state, RREF_MAT_PACKEDBYTES);
                        bad_pk.sf_g[i][j] ^= 1 << below(&mut state, 8);
                    }
                }
            }

            let r = panic::catch_unwind(AssertUnwindSafe(|| {
                less_verify(&bad_pk, &bad_msg, &bad)
            }));

            match r {
                Ok(true) => t.accepted += 1,
                Ok(false) => t.rejected += 1,
                Err(_) => {
                    t.panicked += 1;
                    if t.first_panic.is_none() {
                        t.first_panic = captured.lock().unwrap().take();
                    }
                }
            }
        }

        println!("{name:<36} rejected {:>6}  accepted {:>4}  panicked {:>5}",
                 t.rejected, t.accepted, t.panicked);

        if t.accepted > 0 || t.panicked > 0 {
            flagged.push((name, t));
        }
    }

    let _ = panic::take_hook();
    println!();

    if flagged.is_empty() {
        println!("every forged signature was cleanly rejected, no panics");
        return;
    }

    let mut forged = false;
    for (name, t) in &flagged {
        if t.accepted > 0 {
            println!("FORGERY ACCEPTED: {name}, {} times", t.accepted);
            forged = true;
        }
        if t.panicked > 0 {
            println!("VERIFIER PANICKED: {name}, {} of {trials}", t.panicked);
            if let Some(p) = &t.first_panic {
                println!("    {p}");
            }
        }
    }

    if forged {
        std::process::exit(1);
    }
}
