TODO List:
- rewrite 
```
    let mut fixed_weight_string = [0u8; T];
    sample_challenge(&mut fixed_weight_string, &sig.digest);
```
to
```
    let fixed_weight_string =sample_challenge(&mut fixed_weight_string, &sig.digest, T);
```
- simd gaus
- simd transpose
- TODOs in the code
- less_sign/verify should return proper errors
- `monomoal_from_salt` function is missing
- CI






## KAT

`examples/kat.rs` reads the NIST KAT files from `Utilities/KAT_Generation/KAT`
and verifies every signed message in them, like the `crypto_sign_open` half of
`Reference_Implementation/lib/nist/PQCgenKAT_sign.c`.

Build and run it for the default parameter set (`cat1_192`):

```sh
cargo run --release --example kat
```

The KAT file is picked by the size of the public key, so the selected
parameter set decides which `PQCsignKAT_<bytes>.rsp` is read:

```sh
cargo run --release --no-default-features --features cat1_192  --example kat
cargo run --release --no-default-features --features cat1_68  --example kat
cargo run --release --no-default-features --features cat1_45  --example kat
cargo run --release --no-default-features --features cat3_220 --example kat
cargo run --release --no-default-features --features cat3_102 --example kat
cargo run --release --no-default-features --features cat5_345 --example kat
cargo run --release --no-default-features --features cat5_137 --example kat
```

To build without running it:

```sh
cargo build --release --example kat
```

A different file, and a cap on the number of entries, can be passed after `--`:

```sh
cargo run --release --example kat -- ../../Utilities/KAT_Generation/KAT/PQCsignKAT_13940.rsp 10
```

The example exits non zero if any entry fails to verify. Use `--release`, the
debug build is orders of magnitude slower.
