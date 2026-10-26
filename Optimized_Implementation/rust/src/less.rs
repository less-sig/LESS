#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables)]


/// TODO:
/// type alias: MatrixRREF = Matrix<K, N_pad>
/// type alias: XOF = Xof

use sha3::{
    Digest,
};
use crate::config::{Xof, Hash, N, K, T, N_pad, N_K_pad, W, N8, PRIVATE_KEY_SEED_LENGTH_BYTES, RREF_MAT_PACKEDBYTES, NUM_KEYPAIRS, HASH_DIGEST_LENGTH, SEED_TREE_MAX_PUBLISHED_BYTES, SEED_LENGTH_BYTES, NUM_NODES_SEED_TREE, SIGN_PIVOT_REUSE_LIMIT};
use crate::matrix::{Matrix, MatrixNormalized};
use crate::monomial::{Monomial, Permutation};
use crate::prng::{
    randombytes,
    csprng_randombytes, csprng_initialize
};
use crate::tree::{build_ggm, ggm_path, seed_leaves};
use crate::utils::sample_challenge;

pub struct PrivateKey {
    compressed_sk: [u8; PRIVATE_KEY_SEED_LENGTH_BYTES]
}

/// Public key: the first gen. matrix is shrunk to just a seed, all the
 /// others are stored in RREF form
pub struct PublicKey {
    g0_seed: [u8; SEED_LENGTH_BYTES],
    sf_g: [[u8; RREF_MAT_PACKEDBYTES]; NUM_KEYPAIRS - 1],
}

pub struct Signature {
    digest: [u8; HASH_DIGEST_LENGTH],
    salt: [u8; HASH_DIGEST_LENGTH],
    cf_monom_actions: [[u8; N8]; W],
    seed_storage: [u8; SEED_TREE_MAX_PUBLISHED_BYTES],
}


// LESS_keygen
// \param SK[out]: pointer to an uninitialized secret key data structure
// \param PK[out]: pointer to an uninitialized public key data structure
pub fn less_keygen(sk: &mut PrivateKey, pk: &mut PublicKey) {
    /* generating private key from a single seed */
    randombytes(&mut sk.compressed_sk);

    /* expanding it onto private seeds */
    let mut sk_shake_state = csprng_initialize::<Xof>(&sk.compressed_sk);

    /* Generating public code G_0 */
    csprng_randombytes(&mut pk.g0_seed, &mut sk_shake_state);

    let g0_rref = Matrix::<K, {N - K}, N_K_pad>::rand_from_seed::<Xof>(&pk.g0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);

    let tmp_full_g = Matrix::<K, N, N_pad>::from_compressed_rref(&g0_column_pos, &g0_rref);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    let mut private_monomial_seeds: [[u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1] =
        [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
    for i in 0..(NUM_KEYPAIRS - 1) {
        csprng_randombytes(&mut private_monomial_seeds[i], &mut sk_shake_state);
    }
    let mut was_pivot_column = [0u8; N];
    /* note that the first "keypair" is just the public generator G_0, stored
     * as a seed and the identity matrix (not stored) */
    for i in 0..(NUM_KEYPAIRS - 1) {
        let mut is_pivot_column = [0u8; N];
        /* expand inverse monomial from seed */
        let mut private_q = Monomial::<N, N_pad>::default();
        let private_q_inv = Monomial::<N, N_pad>::rand_from_seed::<Xof>(&private_monomial_seeds[i]);
        // NOTE: the inversion is not implemented in constant time.
        Monomial::<N, N_pad>::inv_non_ct(&mut private_q, &private_q_inv);

        let mut result_g = Matrix::<K, N, N_pad>::default();
        Matrix::monomial_mul(&mut result_g, &tmp_full_g, &private_q); // TODO replace with `*` operator

        let _ = result_g.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, K);
        // note that the result is stored at i-1 as the first
        // public key element is just a seed
        result_g.compress(&mut pk.sf_g[i], &is_pivot_column);
    }
} /* end LESS_keygen */

pub fn less_sign(sig: &mut Signature, m: &[u8], sk: &PrivateKey, pk: &PublicKey) -> usize {
    let mut g0_initial_pivot_flags = [0u8; N];
    let mut is_pivot_column = [0u8; N];

    /*         Private key expansion        */
    /* expand sequence of seeds for private inverse-monomial matrices */
    let mut sk_shake_state = csprng_initialize::<Xof>(&sk.compressed_sk);

    /* Generating seed for public code G_0 */
    let mut G_0_seed = [0u8; SEED_LENGTH_BYTES];
    csprng_randombytes(&mut G_0_seed, &mut sk_shake_state);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    let mut private_monomial_seeds = [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
    for i in 0..(NUM_KEYPAIRS - 1) {
        csprng_randombytes(&mut private_monomial_seeds[i], &mut sk_shake_state);
    }

    // generate the salt from a TRNG
    randombytes(&mut sig.salt);

    /*         Ephemeral monomial generation        */
    let mut ephem_monomials_seed = [0u8; SEED_LENGTH_BYTES];
    csprng_randombytes(&mut ephem_monomials_seed, &mut sk_shake_state);

    /* create the prng for the "blinding" monomials for the canonical form computation */
    let mut cf_seed = [0u8; SEED_LENGTH_BYTES];
    csprng_randombytes(&mut cf_seed, &mut sk_shake_state);
    let mut cf_shake_state = csprng_initialize::<Xof>(&cf_seed);

    let mut seed_tree = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
    build_ggm(&mut seed_tree, &ephem_monomials_seed, &sig.salt);

    let mut linearized_rounds_seeds = [0u8; T*SEED_LENGTH_BYTES];
    seed_leaves(&mut linearized_rounds_seeds, &seed_tree);

    // /*         Public G_0 expansion                  */
    let g0_rref = Matrix::<K, {N - K}, N_K_pad>::rand_from_seed::<Xof>(&G_0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);
    /* generator_get_pivot_flags: everything is a pivot but the stored columns */
    for i in 0..N {
        g0_initial_pivot_flags[i] = 1;
    }
    for i in 0..K {
        g0_initial_pivot_flags[g0_column_pos[i] as usize] = 0;
    }
    let full_g0 = Matrix::<K, N, N_pad>::from_compressed_rref::<{N - K}, N_K_pad>(&g0_column_pos, &g0_rref);

    let mut G0 = Matrix::<K, N, N_pad>::default();
    let mut hasher = Hash::new();
    let mut pi_tilde: [Permutation::<K>; T] = core::array::from_fn(|_| Permutation::<K>::default());

    /* the ephemeral monomial is drawn from the round seed, the salt and the
     * round index, see `monomial_sample_salt` */
    let mut mu_input = [0u8; SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH + 2];
    mu_input[SEED_LENGTH_BYTES..SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH]
        .copy_from_slice(&sig.salt);

    let mut i = 0usize;
    while i < T {
        mu_input[..SEED_LENGTH_BYTES].copy_from_slice(
            &linearized_rounds_seeds[i*SEED_LENGTH_BYTES..(i+1)*SEED_LENGTH_BYTES]);
        mu_input[SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH..]
            .copy_from_slice(&(i as u16).to_le_bytes());

        let mut mu_state = csprng_initialize::<Xof>(&mu_input);
        let mu_tilde = Monomial::<N, N_pad>::rand(&mut mu_state);
        Matrix::monomial_mul(&mut G0, &full_g0, &mu_tilde);

        is_pivot_column = [0u8; N];

        let mut permuted_pivot_flag: [u8; N] = [0u8; N];
        for t in 0..N {
            permuted_pivot_flag[mu_tilde.perms[t] as usize] = g0_initial_pivot_flags[t];
        }

        if G0.rref::<true>(&mut is_pivot_column, &mut permuted_pivot_flag, SIGN_PIVOT_REUSE_LIMIT) == false {
            return 0;
        }
        // just copy the non IS
        let mut a_j = MatrixNormalized::<{N - K}, N_K_pad>::from_information_set::<K, N, N_pad>(&G0, &is_pivot_column);
        pi_tilde[i] = Permutation::<K>::from_information_set_composition::<N, N_pad, N>(&mu_tilde, &is_pivot_column);


        a_j.blind(&mut cf_shake_state);
        let t = a_j.cf();

        if t == 0 {
            // this gives us 2**64 retries. That should be more than enough.
            let s = &mut linearized_rounds_seeds[i*SEED_LENGTH_BYTES..i*SEED_LENGTH_BYTES + 8];
            let v = u64::from_le_bytes(s.try_into().unwrap()).wrapping_add(1);
            s.copy_from_slice(&v.to_le_bytes());
            continue;
        }

        // NOTE: as we increase the size of the `normalized_IS_t`
        // we need to hash the values row by row.
        let mut row_bytes = [0u8; N - K];
        for sl in 0..K {
            for col in 0..(N - K) {
                row_bytes[col] = a_j[sl][col].0;
            }
            hasher.update(&row_bytes);
        }

        i += 1;
    }
    hasher.update(m);
    hasher.update(sig.salt);
    sig.digest = *hasher.finalize().as_array().unwrap();


    let mut fixed_weight_string = [0u8; T];
    sample_challenge(&mut fixed_weight_string, &sig.digest);

    let mut indices_to_publish = [0u8; T];
    for i in 0..T {
        indices_to_publish[i] = (fixed_weight_string[i] != 0) as u8;
    }

    sig.seed_storage = [0u8; SEED_TREE_MAX_PUBLISHED_BYTES];

    let num_seeds_publised = ggm_path(&seed_tree, &indices_to_publish, &mut sig.seed_storage);
    let mut emitted_monoms = 0;
    for i in 0..T {
        if fixed_weight_string[i] != 0 {
            let j = fixed_weight_string[i] as usize;
            let Q_to_multiply = Monomial::<N, N_pad>::rand_from_seed::<Xof>(&private_monomial_seeds[j - 1]);
            let mono_action = Permutation::<K>::from_composition_action(&Q_to_multiply, &pi_tilde[i]);
            sig.cf_monom_actions[emitted_monoms] = mono_action.bytes();

            emitted_monoms += 1;
        }
    }

    num_seeds_publised
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::prng::{init_randombytes, PRNG_TEST_LOCK};

    /// The vectors the C implementation produces for the same seed, generated
    /// by `kat/run.sh`. Both sides drive the SHAKE based `platform_csprng_state`
    /// from the seed below, so every byte has to agree.
    #[cfg(feature = "cat1_192")]
    const REFERENCE: &str = include_str!("../kat/reference_cat1_192.txt");

    #[cfg(feature = "cat1_192")]
    fn field(name: &str) -> &'static str {
        for line in REFERENCE.lines() {
            if let Some((k, v)) = line.split_once(" = ") {
                if k == name {
                    return v.trim();
                }
            }
        }
        panic!("no field `{}` in the reference vectors", name);
    }

    #[cfg(feature = "cat1_192")]
    fn hex(bytes: &[u8]) -> String {
        let mut s = String::new();
        for b in bytes {
            s.push_str(&format!("{:02x}", b));
        }
        s
    }

    /// the parameters the reference vectors were generated with
    #[cfg(feature = "cat1_192")]
    #[test]
    fn reference_parameters() {
        assert_eq!(field("N"), N.to_string());
        assert_eq!(field("K"), K.to_string());
        assert_eq!(field("T"), T.to_string());
        assert_eq!(field("W"), W.to_string());
        assert_eq!(field("SEED_LENGTH_BYTES"), SEED_LENGTH_BYTES.to_string());
        assert_eq!(field("RREF_MAT_PACKEDBYTES"), RREF_MAT_PACKEDBYTES.to_string());
    }

    #[cfg(feature = "cat1_192")]
    #[test]
    fn keygen_and_sign_match_the_reference() {
        let _guard = PRNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());

        let seed: [u8; 48] = core::array::from_fn(|i| i as u8);
        init_randombytes(&seed);

        let mut sk = PrivateKey { compressed_sk: [0u8; PRIVATE_KEY_SEED_LENGTH_BYTES] };
        let mut pk = PublicKey {
            g0_seed: [0u8; SEED_LENGTH_BYTES],
            sf_g: [[0u8; RREF_MAT_PACKEDBYTES]; NUM_KEYPAIRS - 1],
        };
        less_keygen(&mut sk, &mut pk);

        assert_eq!(hex(&sk.compressed_sk), field("sk"), "secret key");

        // `pubkey_t` is a packed struct: G_0_seed followed by SF_G
        let mut pk_bytes: Vec<u8> = Vec::new();
        pk_bytes.extend_from_slice(&pk.g0_seed);
        for i in 0..(NUM_KEYPAIRS - 1) {
            pk_bytes.extend_from_slice(&pk.sf_g[i]);
        }
        assert_eq!(hex(&pk_bytes), field("pk"), "public key");

        let msg = b"LESS rust port reference vector";
        assert_eq!(hex(msg), field("msg"), "message");

        let mut sig = Signature {
            digest: [0u8; HASH_DIGEST_LENGTH],
            salt: [0u8; HASH_DIGEST_LENGTH],
            cf_monom_actions: [[0u8; N8]; W],
            seed_storage: [0u8; SEED_TREE_MAX_PUBLISHED_BYTES],
        };
        let leaves = less_sign(&mut sig, msg, &sk, &pk);

        assert_eq!(hex(&sig.salt), field("salt"), "salt");
        assert_eq!(hex(&sig.digest), field("digest"), "digest");

        let mut actions: Vec<u8> = Vec::new();
        for i in 0..W {
            actions.extend_from_slice(&sig.cf_monom_actions[i]);
        }
        assert_eq!(hex(&actions), field("cf_monom_actions"), "cf_monom_actions");
        assert_eq!(hex(&sig.seed_storage), field("seed_storage"), "seed_storage");
        assert_eq!(leaves.to_string(), field("leaves"), "published leaves");
    }
}
