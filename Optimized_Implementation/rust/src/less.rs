#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables)]


/// TODO:
/// type alias: MatrixRREF = Matrix<K, N_pad>

use sha3::{
    Digest,
};
use crate::config::{
    Xof, Hash, N, K, T, N_pad, N_K_pad, W, N8, 
    PRIVATE_KEY_SEED_LENGTH_BYTES, RREF_MAT_PACKEDBYTES, NUM_KEYPAIRS, 
    HASH_DIGEST_LENGTH, SEED_TREE_MAX_PUBLISHED_BYTES, SEED_LENGTH_BYTES, 
    NUM_NODES_SEED_TREE, SIGN_PIVOT_REUSE_LIMIT, LESS_VERIFY_PIVOT_REUSE_LIMIT
};
use crate::matrix::{Matrix, MatrixNormalized};
use crate::monomial::{Monomial, Permutation};
use crate::prng::{
    randombytes,
    csprng_randombytes, csprng_initialize
};
use crate::tree::{build_ggm, ggm_path, rebuild_ggm, seed_leaves};
use crate::utils::{sample_challenge, verify};

pub struct PrivateKey {
    pub compressed_sk: [u8; PRIVATE_KEY_SEED_LENGTH_BYTES]
}

/// Public key: the first gen. matrix is shrunk to just a seed, all the
 /// others are stored in RREF form
pub struct PublicKey {
    pub g0_seed: [u8; SEED_LENGTH_BYTES],
    pub sf_g: [[u8; RREF_MAT_PACKEDBYTES]; NUM_KEYPAIRS - 1],
}

pub struct Signature {
    pub digest: [u8; HASH_DIGEST_LENGTH],
    pub salt: [u8; HASH_DIGEST_LENGTH],
    pub cf_monom_actions: [[u8; N8]; W],
    pub seed_storage: [u8; SEED_TREE_MAX_PUBLISHED_BYTES],
}

type MatrixRREF = Matrix<K, {N - K}, N_K_pad>;
type MatrixFull = Matrix<K, N, N_pad>;
type MonomialFull = Monomial::<N, N_pad>;
type NonIS = MatrixNormalized::<{N - K}, N_K_pad>;


// LESS_keygen
// \param SK[out]: pointer to an uninitialized secret key data structure
// \param PK[out]: pointer to an uninitialized public key data structure
pub fn less_keygen(sk: &mut PrivateKey, pk: &mut PublicKey) {
    /* generating private key from a single seed */
    randombytes(&mut sk.compressed_sk);

    /* expanding it onto private seeds */
    let mut sk_shake_state = csprng_initialize::<Xof>(&sk.compressed_sk);

    /* Generating public code G_0 */
    pk.g0_seed
        .copy_from_slice(&csprng_randombytes(&mut sk_shake_state, SEED_LENGTH_BYTES));

    let g0_rref = MatrixRREF::rand_from_seed::<Xof>(&pk.g0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);

    let tmp_full_g = MatrixFull::from_compressed_rref(&g0_column_pos, &g0_rref);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    let mut private_monomial_seeds: [[u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1] =
        [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
    for i in 0..(NUM_KEYPAIRS - 1) {
        private_monomial_seeds[i].copy_from_slice(&csprng_randombytes(
            &mut sk_shake_state,
            PRIVATE_KEY_SEED_LENGTH_BYTES,
        ));
    }
    let mut was_pivot_column = [0u8; N];
    /* note that the first "keypair" is just the public generator G_0, stored
     * as a seed and the identity matrix (not stored) */
    for i in 0..(NUM_KEYPAIRS - 1) {
        let mut is_pivot_column = [0u8; N];
        /* expand inverse monomial from seed */
        let mut private_q = MonomialFull::default();
        let private_q_inv = MonomialFull::rand_from_seed::<Xof>(&private_monomial_seeds[i]);
        // NOTE: the inversion is not implemented in constant time.
        MonomialFull::inv_non_ct(&mut private_q, &private_q_inv);

        let mut result_g = MatrixFull::default();
        Matrix::monomial_mul(&mut result_g, &tmp_full_g, &private_q); // TODO replace with `*` operator

        let _ = result_g.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, K);
        // note that the result is stored at i-1 as the first
        // public key element is just a seed
        result_g.compress(&mut pk.sf_g[i], &is_pivot_column);
    }
} /* end LESS_keygen */

pub fn less_sign(sig: &mut Signature, m: &[u8], sk: &PrivateKey, pk: &PublicKey) -> usize {
    let mut g0_initial_pivot_flags = [0u8; N];
    let mut is_pivot_column;

    /*         Private key expansion        */
    /* expand sequence of seeds for private inverse-monomial matrices */
    let mut sk_shake_state = csprng_initialize::<Xof>(&sk.compressed_sk);

    /* Generating seed for public code G_0 */
    let G_0_seed = csprng_randombytes(&mut sk_shake_state, SEED_LENGTH_BYTES);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    let mut private_monomial_seeds = [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
    for i in 0..(NUM_KEYPAIRS - 1) {
        private_monomial_seeds[i].copy_from_slice(&csprng_randombytes(
            &mut sk_shake_state,
            PRIVATE_KEY_SEED_LENGTH_BYTES,
        ));
    }

    // generate the salt from a TRNG
    randombytes(&mut sig.salt);

    /*         Ephemeral monomial generation        */
    let ephem_monomials_seed = csprng_randombytes(&mut sk_shake_state, SEED_LENGTH_BYTES);

    /* create the prng for the "blinding" monomials for the canonical form computation */
    let cf_seed = csprng_randombytes(&mut sk_shake_state, SEED_LENGTH_BYTES);
    let mut cf_shake_state = csprng_initialize::<Xof>(&cf_seed);

    let mut seed_tree = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
    build_ggm(
        &mut seed_tree,
        ephem_monomials_seed.as_slice().try_into().unwrap(),
        &sig.salt,
    );

    let mut linearized_rounds_seeds = [0u8; T*SEED_LENGTH_BYTES];
    seed_leaves(&mut linearized_rounds_seeds, &seed_tree);

    // /*         Public G_0 expansion                  */
    let g0_rref = MatrixRREF::rand_from_seed::<Xof>(&G_0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);
    /* generator_get_pivot_flags: everything is a pivot but the stored columns */
    for i in 0..N {
        g0_initial_pivot_flags[i] = 1;
    }
    for i in 0..K {
        g0_initial_pivot_flags[g0_column_pos[i] as usize] = 0;
    }
    let full_g0 = MatrixFull::from_compressed_rref::<{N - K}, N_K_pad>(&g0_column_pos, &g0_rref);

    let mut G0 = MatrixFull::default();
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
        let mu_tilde = MonomialFull::rand(&mut mu_state);
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
        let mut a_j = NonIS::from_information_set::<K, N, N_pad>(&G0, &is_pivot_column);
        pi_tilde[i] = Permutation::<K>::from_information_set_composition::<N, N_pad>(&mu_tilde, &is_pivot_column);


        a_j.blind(&mut cf_shake_state);
        let t = a_j.cf();

        if t == 0 {
            // this gives us 2**64 retries. That should be more than enough.
            let s = &mut linearized_rounds_seeds[i*SEED_LENGTH_BYTES..i*SEED_LENGTH_BYTES + 8];
            let v = u64::from_le_bytes(s.try_into().unwrap()).wrapping_add(1);
            s.copy_from_slice(&v.to_le_bytes());
            continue;
        }

        // NOTE: as we added padding: we need to hash the values row by row.
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
            let Q_to_multiply = MonomialFull::rand_from_seed::<Xof>(&private_monomial_seeds[j - 1]);
            let mono_action = Permutation::<K>::from_composition_action(&Q_to_multiply, &pi_tilde[i]);
            sig.cf_monom_actions[emitted_monoms] = mono_action.bytes();

            emitted_monoms += 1;
        }
    }

    num_seeds_publised
}


/// NOTE: non-constant time
/// TODO: example
/// \param pk[in]: public key
/// \param m[in]: message for which a signature was computed
/// \param sig[in]: signature
/// \param num_seeds_published[in]: number of seeds in `sig.seed_storage`
/// \return false: on failure
///         true: on success
pub fn less_verify(pk: &PublicKey, m: &[u8], sig: &Signature) -> bool {
    let mut fixed_weight_string = [0u8; T];
    sample_challenge(&mut fixed_weight_string, &sig.digest);

    let mut published_seed_indexes = [0u8; T];
    for i in 0..T {
        published_seed_indexes[i] = (fixed_weight_string[i] != 0) as u8;
    }

    let mut seed_tree = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
    let rebuilding_seeds_went_fine =
        rebuild_ggm(&mut seed_tree, &published_seed_indexes, &sig.seed_storage, &sig.salt);
    if rebuilding_seeds_went_fine == 0 {
        return false;
    }

    let mut linearized_rounds_seeds = [0u8; T * SEED_LENGTH_BYTES];
    seed_leaves(&mut linearized_rounds_seeds, &seed_tree);

    let mut employed_monoms = 0usize;

    let g0_rref = MatrixRREF::rand_from_seed::<Xof>(&pk.g0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);

    let mut g0_initial_pivot_flags = [0u8; N];
    for i in 0..N {
        g0_initial_pivot_flags[i] = 1;
    }
    for i in 0..K {
        g0_initial_pivot_flags[g0_column_pos[i] as usize] = 0;
    }
    let g0_full = MatrixFull::from_compressed_rref::<{N - K}, N_K_pad>(
        &g0_column_pos, &g0_rref);

    let mut gs: Vec<Matrix<K, N, N_pad>> = Vec::with_capacity(NUM_KEYPAIRS - 1);
    let mut gi_initial_pivot_flags = [0u8; N];
    for i in 0..(NUM_KEYPAIRS - 1) {
        let mut g = MatrixFull::default();
        g.expand(&pk.sf_g[i], &mut gi_initial_pivot_flags);
        gs.push(g);
    }

    let mut g_prime = MatrixFull::default();
    let mut hasher = Hash::new();

    /* same ephemeral monomial derivation as in `less_sign` */
    let mut mu_input = [0u8; SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH + 2];
    mu_input[SEED_LENGTH_BYTES..SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH]
        .copy_from_slice(&sig.salt);

    for i in 0..T {
        let mut is_pivot_column = [0u8; N];

        if fixed_weight_string[i] == 0 {
            mu_input[..SEED_LENGTH_BYTES].copy_from_slice(
                &linearized_rounds_seeds[i*SEED_LENGTH_BYTES..(i+1)*SEED_LENGTH_BYTES]);
            mu_input[SEED_LENGTH_BYTES + HASH_DIGEST_LENGTH..]
                .copy_from_slice(&(i as u16).to_le_bytes());

            let mut mu_state = csprng_initialize::<Xof>(&mu_input);
            let mu_tilde = MonomialFull::rand(&mut mu_state);

            Matrix::monomial_mul(&mut g_prime, &g0_full, &mu_tilde);

            let mut permuted_pivot_flags = [0u8; N];
            for t in 0..N {
                permuted_pivot_flags[mu_tilde.perms[t] as usize] = g0_initial_pivot_flags[t];
            }

            if g_prime.rref::<false>(&mut is_pivot_column, &mut permuted_pivot_flags,
                                     LESS_VERIFY_PIVOT_REUSE_LIMIT) == false {
                return false;
            }
        } else {
            gi_initial_pivot_flags = [0u8; N];
            let key_pos = (fixed_weight_string[i] - 1) as usize;

            if !Permutation::<K>::check_canonical_action(&sig.cf_monom_actions[employed_monoms]) {
                return false;
            }

            let mut g0_permuted_pivot_flags = [0u8; N];
            Matrix::unpack_coset_rep(&mut g_prime, &gs[key_pos],
                                     &sig.cf_monom_actions[employed_monoms],
                                     &gi_initial_pivot_flags,
                                     &mut g0_permuted_pivot_flags);

            if g_prime.rref::<false>(&mut is_pivot_column, &mut g0_permuted_pivot_flags,
                                     LESS_VERIFY_PIVOT_REUSE_LIMIT) == false {
                return false;
            }

            employed_monoms += 1;
        }

        // copy the non IS
        let mut a_i = NonIS::from_information_set::<K, N, N_pad>(
            &g_prime, &is_pivot_column);

        // compute the canonical form and exit on any failure
        if a_i.cf() == 0 {
            return false;
        }

        // NOTE: as we increase the size of the `normalized_IS_t`
        // we need to hash the values row by row.
        let mut row_bytes = [0u8; N - K];
        for sl in 0..K {
            for col in 0..(N - K) {
                row_bytes[col] = a_i[sl][col].0;
            }
            hasher.update(&row_bytes);
        }
    }

    hasher.update(m);
    hasher.update(sig.salt);
    let recomputed_digest = hasher.finalize();

    // NOTE: does not have to be constant time.
    verify(&recomputed_digest, &sig.digest) == 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::MAX_PUBLISHED_SEEDS;
    use crate::prng::{init_randombytes, PRNG_TEST_LOCK};

    /// keygen and sign once, from a fixed seed
    fn keypair_and_signature(msg: &[u8]) -> (PublicKey, Signature) {
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
        assert!(leaves > 0 && leaves <= MAX_PUBLISHED_SEEDS,
                "{} published leaves is out of range", leaves);

        (pk, sig)
    }

    /// a copy of a signature, so that a single field can be tampered with
    fn copy(sig: &Signature) -> Signature {
        Signature {
            digest: sig.digest,
            salt: sig.salt,
            cf_monom_actions: sig.cf_monom_actions,
            seed_storage: sig.seed_storage,
        }
    }

    /// a signature this port produces has to verify
    #[test]
    fn sign_and_verify_round_trip() {
        let _guard = PRNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());

        let msg = b"abc";
        let (pk, sig) = keypair_and_signature(msg);

        assert!(less_verify(&pk, msg, &sig), "a valid signature was rejected");
    }

    /// signing is deterministic for a fixed prng state
    #[test]
    fn sign_is_deterministic() {
        let _guard = PRNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());

        let msg = b"abc";
        let (_, a) = keypair_and_signature(msg);
        let (_, b) = keypair_and_signature(msg);

        assert_eq!(a.digest, b.digest);
        assert_eq!(a.salt, b.salt);
        assert_eq!(a.cf_monom_actions, b.cf_monom_actions);
        assert_eq!(a.seed_storage, b.seed_storage);
    }

    /// every part of the signature matters
    #[test]
    fn verify_rejects_a_tampered_signature() {
        let _guard = PRNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());

        let msg = b"def";
        let (pk, sig) = keypair_and_signature(msg);

        // a different message must not verify
        assert!(!less_verify(&pk, b"ghi", &sig));

        let mut t = copy(&sig);
        t.digest[0] ^= 1;
        assert!(!less_verify(&pk, msg, &t), "a flipped digest bit was accepted");

        let mut t = copy(&sig);
        t.salt[0] ^= 1;
        assert!(!less_verify(&pk, msg, &t), "a flipped salt bit was accepted");

        let mut t = copy(&sig);
        t.cf_monom_actions[0][0] ^= 1;
        assert!(!less_verify(&pk, msg, &t), "a flipped action bit was accepted");

        let mut t = copy(&sig);
        t.seed_storage[0] ^= 1;
        assert!(!less_verify(&pk, msg, &t), "a flipped seed bit was accepted");
    }
}

