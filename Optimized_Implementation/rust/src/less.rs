#![allow(dead_code)]


/// TODO:
/// type alias: MatrixRREF = Matrix<K, N_pad>
/// type alias: XOF = Shake128

use sha3::{
    Digest, Shake128, Sha3_256,
};
use crate::config::{N, K, T, N_pad, K_pad, N_K_pad, W, N8, PRIVATE_KEY_SEED_LENGTH_BYTES, RREF_MAT_PACKEDBYTES, NUM_KEYPAIRS, HASH_DIGEST_LENGTH, SEED_TREE_MAX_PUBLISHED_BYTES, SEED_LENGTH_BYTES, NUM_NODES_SEED_TREE, SIGN_PIVOT_REUSE_LIMIT};
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
    g0_seed: [u8; PRIVATE_KEY_SEED_LENGTH_BYTES],
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
    let mut sk_shake_state = csprng_initialize::<Shake128>(&sk.compressed_sk);

    /* Generating public code G_0 */
    csprng_randombytes(&mut pk.g0_seed, &mut sk_shake_state);

    let g0_rref = Matrix::<K, N_K_pad>::rand_from_seed::<Shake128>(&pk.g0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);

    let tmp_full_g = Matrix::<K, N_pad>::from_compressed_rref(&g0_column_pos, &g0_rref);

    /* The first private key monomial is an ID matrix, no need for random
     * generation, hence NUM_KEYPAIRS-1 */
    let mut private_monomial_seeds: [[u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1] =
        [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
    for i in 0..(NUM_KEYPAIRS - 1) {
        csprng_randombytes(&mut private_monomial_seeds[i], &mut sk_shake_state);
    }
    let mut was_pivot_column = [0u8; N_pad];
    /* note that the first "keypair" is just the public generator G_0, stored
     * as a seed and the identity matrix (not stored) */
    for i in 0..(NUM_KEYPAIRS - 1) {
        let mut is_pivot_column = [0u8; N_pad];
        /* expand inverse monomial from seed */
        let mut private_q = Monomial::<N>::default();
        let private_q_inv = Monomial::<N>::rand_from_seed::<Shake128>(&private_monomial_seeds[i]);
        // NOTE: the inversion is not implemented in constant time.
        Monomial::<N>::inv_non_ct(&mut private_q, &private_q_inv);

        let mut result_g = Matrix::<K, N_pad>::default();
        Matrix::monomial_mul(&mut result_g, &tmp_full_g, &private_q); // TODO replace with `*` operator

        result_g.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, K);
        // note that the result is stored at i-1 as the first
        // public key element is just a seed
        result_g.compress(&mut pk.sf_g[i], &is_pivot_column);
    }
} /* end LESS_keygen */

pub fn less_sign(sig: &mut Signature, m: &[u8], sk: &PrivateKey, pk: &PublicKey) -> usize {
    let mut g0_initial_pivot_flags = [0u8; N_pad];
    let mut is_pivot_column = [0u8; N_pad];

    /*         Private key expansion        */
    /* expand sequence of seeds for private inverse-monomial matrices */
    let mut sk_shake_state = csprng_initialize::<Shake128>(&sk.compressed_sk);

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
    let mut cf_shake_state = csprng_initialize::<Shake128>(&cf_seed);

    let mut seed_tree = [0u8; NUM_NODES_SEED_TREE * SEED_LENGTH_BYTES];
    build_ggm(&mut seed_tree, &ephem_monomials_seed, &sig.salt);

    let mut linearized_rounds_seeds = [0u8; T*SEED_LENGTH_BYTES];
    seed_leaves(&mut linearized_rounds_seeds, &seed_tree);

    // /*         Public G_0 expansion                  */
    let mut g0_rref = Matrix::<K, N_K_pad>::rand_from_seed::<Shake128>(&G_0_seed);
    let g0_column_pos: [u16; N_K_pad] = std::array::from_fn(|i| (i + K) as u16);
    for i in 0..K {
        g0_initial_pivot_flags[i] = 0;
    }
    let full_g0 = Matrix::<K, N_pad>::from_compressed_rref::<N_K_pad>(&g0_column_pos, &g0_rref);

    let mut G0 = Matrix::<K, N_pad>::default();
    let mut hasher = Sha3_256::new();
    let mut pi_tilde: [Permutation::<K>; T] = core::array::from_fn(|_| Permutation::<K>::default());
    for i in 0..T {
        let mu_tilde = Monomial::<N>::rand_from_seed::<Shake128>(&linearized_rounds_seeds[i*SEED_LENGTH_BYTES..(i+1)*SEED_LENGTH_BYTES]);
        Matrix::monomial_mul(&mut G0, &full_g0, &mu_tilde);

        let mut permuted_pivot_flag: [u8; N_pad] = [0u8; N_pad];
        for t in 0..N {
            permuted_pivot_flag[mu_tilde.perms[t] as usize] = g0_initial_pivot_flags[t];
        }

        G0.rref::<true>(&mut is_pivot_column, &mut permuted_pivot_flag, SIGN_PIVOT_REUSE_LIMIT);
        let mut a_j = MatrixNormalized::<K_pad>::from_information_set::<K, N_pad>(&G0, &is_pivot_column);
        pi_tilde[i] = Permutation::<K>::from_information_set_composition::<N, N_pad>(&mu_tilde, &is_pivot_column);


        a_j.blind(&mut cf_shake_state);
        a_j.cf();

        for i in 0..K {
            // TODO without unsafe?
            let bytes: &[u8] = unsafe {
                std::slice::from_raw_parts(a_j[i].0.as_ptr() as *const u8, 126)
            };
            hasher.update(bytes);
        }

    }
    hasher.update(m);
    hasher.update(sig.salt);
    sig.digest = *hasher.finalize().as_array().unwrap();


    let mut fixed_weight_string = [0u8; T];
    sample_challenge(&mut fixed_weight_string, &sig.digest);

    let num_seeds_publised = ggm_path(&seed_tree, &fixed_weight_string, &mut sig.seed_storage);
    let mut emitted_monoms = 0;
    for i in 0..T {
        if fixed_weight_string[i] != 0 {
            let j = fixed_weight_string[i] as usize;
            let Q_to_multiply = Monomial::<N>::rand_from_seed::<Shake128>(&private_monomial_seeds[j]);
            let mono_action = Permutation::<K>::from_composition_action(&Q_to_multiply, &pi_tilde[i]);
            sig.cf_monom_actions[emitted_monoms] = mono_action.bytes();

            emitted_monoms += 1;
        }
    }

    num_seeds_publised
}


#[cfg(test)]
mod tests {
    #[test]
    fn keygen() {

    }
}
