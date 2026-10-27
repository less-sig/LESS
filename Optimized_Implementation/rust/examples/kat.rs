//! Reads the NIST KAT files of `Utilities/KAT_Generation/KAT` and verifies
//! every signed message in them. This is the verification half of
//! `Reference_Implementation/lib/nist/PQCgenKAT_sign.c`
//!
//! usage:
//!     cargo run --release --example kat
//!     cargo run --release --no-default-features --features cat3_220 --example kat
//!     cargo run --release --example kat -- <path/to.rsp> [max_entries]

use std::env;
use std::fs;
use std::process::ExitCode;

use less::config::{
    HASH_DIGEST_LENGTH, LESS_CRYPTO_PUBLICKEYBYTES, N8, NUM_KEYPAIRS, RREF_MAT_PACKEDBYTES,
    SEED_LENGTH_BYTES, SEED_TREE_MAX_PUBLISHED_BYTES, W,
};
use less::less::{less_verify, PublicKey, Signature};

/// one `count = ...` block of a `.rsp` file
struct Entry {
    count: usize,
    mlen: usize,
    pk: Vec<u8>,
    sm: Vec<u8>,
    smlen: usize,
}

fn from_hex(s: &str) -> Result<Vec<u8>, String> {
    let s = s.trim();
    if s.len() % 2 != 0 {
        return Err(format!("odd length hex string of {} chars", s.len()));
    }
    (0..s.len() / 2)
        .map(|i| {
            u8::from_str_radix(&s[2 * i..2 * i + 2], 16)
                .map_err(|_| format!("`{}` is not hex", &s[2 * i..2 * i + 2]))
        })
        .collect()
}

/// splits a `.rsp` file into its `count = ...` blocks
fn parse(text: &str) -> Result<Vec<Entry>, String> {
    let mut out: Vec<Entry> = Vec::new();

    let mut count = 0usize;
    let mut mlen = 0usize;
    let mut smlen = 0usize;
    let mut pk: Vec<u8> = Vec::new();
    let mut sm: Vec<u8> = Vec::new();
    let mut open = false;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let (key, value) = match line.split_once('=') {
            Some((k, v)) => (k.trim(), v.trim()),
            None => continue,
        };

        match key {
            "count" => {
                if open {
                    out.push(Entry { count, mlen, pk: std::mem::take(&mut pk),
                                     sm: std::mem::take(&mut sm), smlen });
                }
                count = value.parse().map_err(|_| format!("bad count `{}`", value))?;
                open = true;
            }
            "mlen" => mlen = value.parse().map_err(|_| format!("bad mlen `{}`", value))?,
            "smlen" => smlen = value.parse().map_err(|_| format!("bad smlen `{}`", value))?,
            "pk" => pk = from_hex(value)?,
            "sm" => sm = from_hex(value)?,
            // `seed`, `msg` and `sk` are only needed to regenerate the KAT
            _ => {}
        }
    }

    if open {
        out.push(Entry { count, mlen, pk, sm, smlen });
    }

    Ok(out)
}

/// splits `sm` the way `crypto_sign_open` does and verifies it
fn verify_entry(e: &Entry) -> Result<(), String> {
    if e.sm.len() != e.smlen {
        return Err(format!("sm is {} bytes, smlen says {}", e.sm.len(), e.smlen));
    }
    if e.pk.len() != LESS_CRYPTO_PUBLICKEYBYTES {
        return Err(format!("pk is {} bytes, expected {}", e.pk.len(),
                           LESS_CRYPTO_PUBLICKEYBYTES));
    }

    /* the number of published seeds is the very last byte, see `crypto_sign` */
    let leaves = e.sm[e.smlen - 1] as usize;
    let published = leaves * SEED_LENGTH_BYTES;

    /* LESS_SIGNATURE_SIZE(leaves) */
    let sig_len = HASH_DIGEST_LENGTH * 2 + N8 * W + published + 1;
    if sig_len > e.smlen {
        return Err(format!("{} published leaves do not fit into {} bytes", leaves, e.smlen));
    }
    if published > SEED_TREE_MAX_PUBLISHED_BYTES {
        return Err(format!("{} published leaves exceed the maximum", leaves));
    }

    let mlen = e.smlen - sig_len;
    if mlen != e.mlen {
        return Err(format!("recovered message is {} bytes, mlen says {}", mlen, e.mlen));
    }

    let msg = &e.sm[..mlen];

    let mut pk = PublicKey {
        g0_seed: [0u8; SEED_LENGTH_BYTES],
        sf_g: [[0u8; RREF_MAT_PACKEDBYTES]; NUM_KEYPAIRS - 1],
    };
    pk.g0_seed.copy_from_slice(&e.pk[..SEED_LENGTH_BYTES]);
    for i in 0..(NUM_KEYPAIRS - 1) {
        let off = SEED_LENGTH_BYTES + i * RREF_MAT_PACKEDBYTES;
        pk.sf_g[i].copy_from_slice(&e.pk[off..off + RREF_MAT_PACKEDBYTES]);
    }

    let mut sig = Signature {
        digest: [0u8; HASH_DIGEST_LENGTH],
        salt: [0u8; HASH_DIGEST_LENGTH],
        cf_monom_actions: [[0u8; N8]; W],
        seed_storage: [0u8; SEED_TREE_MAX_PUBLISHED_BYTES],
    };

    let mut off = mlen;
    sig.digest.copy_from_slice(&e.sm[off..off + HASH_DIGEST_LENGTH]);
    off += HASH_DIGEST_LENGTH;
    sig.salt.copy_from_slice(&e.sm[off..off + HASH_DIGEST_LENGTH]);
    off += HASH_DIGEST_LENGTH;
    for i in 0..W {
        sig.cf_monom_actions[i].copy_from_slice(&e.sm[off..off + N8]);
        off += N8;
    }
    /* only the published seeds are stored, the rest of the buffer stays zero */
    sig.seed_storage[..published].copy_from_slice(&e.sm[off..off + published]);

    if !less_verify(&pk, msg, &sig) {
        return Err("less_verify rejected the signature".to_string());
    }

    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();

    let path = if args.len() > 1 {
        args[1].clone()
    } else {
        format!("../../Utilities/KAT_Generation/KAT/PQCsignKAT_{}.rsp",
                LESS_CRYPTO_PUBLICKEYBYTES)
    };
    let limit: usize = if args.len() > 2 {
        args[2].parse().unwrap_or(usize::MAX)
    } else {
        usize::MAX
    };

    let text = match fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("could not open <{}>: {}", path, e);
            return ExitCode::from(1);
        }
    };

    let entries = match parse(&text) {
        Ok(e) => e,
        Err(e) => {
            eprintln!("could not parse <{}>: {}", path, e);
            return ExitCode::from(3);
        }
    };

    println!("{}", path);
    println!("public key bytes = {}", LESS_CRYPTO_PUBLICKEYBYTES);
    println!("entries = {}", entries.len());

    let mut ok = 0usize;
    let mut failed = 0usize;

    for e in entries.iter().take(limit) {
        match verify_entry(e) {
            Ok(()) => ok += 1,
            Err(msg) => {
                failed += 1;
                println!("count = {:<4} FAILED: {}", e.count, msg);
            }
        }
    }

    println!("verified = {}, failed = {}", ok, failed);

    if failed == 0 && ok > 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(4)
    }
}
