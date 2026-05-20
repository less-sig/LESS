#![allow(dead_code)]

use crate::config::{
    N, K, T, W, N8,
    PRIVATE_KEY_SEED_LENGTH_BYTES,
    RREF_MAT_PACKEDBYTES, 
    NUM_KEYPAIRS,
    HASH_DIGEST_LENGTH,
    SEED_TREE_MAX_PUBLISHED_BYTES
};
use crate::monomial::Monomial;
use crate::matrix::Matrix;
use crate::prng::randombytes;



pub struct PrivateKey {
    compressed_sk: [u8; PRIVATE_KEY_SEED_LENGTH_BYTES]
}

/// Public key: the first gen. matrix is shrunk to just a seed, all the
 /// others are stored in RREF form
pub struct PublicKey {
    g_0_seed: [u8; PRIVATE_KEY_SEED_LENGTH_BYTES],
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
//pub fn less_keygen(sk: &mut PrivateKey, pk: &mut PublicKey) {
//    /* generating private key from a single seed */
//    randombytes(&mut sk.compressed_sk);
//
//    /* expanding it onto private seeds */
//    let mut sk_shake_state = ShakeState::default();
//    initialize_csprng(&mut sk_shake_state, &sk.compressed_sk);
//
//    /* Generating public code G_0 */
//    csprng_randombytes(&mut pk.g0_seed, &mut sk_shake_state);
//
//    let mut g0_rref = rref_generator_mat_t::default();
//    generator_sample(&mut g0_rref, &pk.g0_seed);
//
//    let mut tmp_full_g = Matrix::<K, N>::default();
//    generator_rref_expand(&mut tmp_full_g, &g0_rref);
//
//    /* The first private key monomial is an ID matrix, no need for random
//     * generation, hence NUM_KEYPAIRS-1 */
//    let mut private_monomial_seeds: [[u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1] =
//        [[0u8; PRIVATE_KEY_SEED_LENGTH_BYTES]; NUM_KEYPAIRS - 1];
//    for i in 0..(NUM_KEYPAIRS - 1) {
//        csprng_randombytes(&mut private_monomial_seeds[i], &mut sk_shake_state);
//    }
//    let mut was_pivot_column = [0u8; N_pad];
//    /* note that the first "keypair" is just the public generator G_0, stored
//     * as a seed and the identity matrix (not stored) */
//    for i in 0..(NUM_KEYPAIRS - 1) {
//        let mut is_pivot_column = [0u8; N_pad];
//        /* expand inverse monomial from seed */
//        let mut private_q = Monomial::<N>::default();
//        let mut private_q_inv = Monomial::<N>::default();
//        monomial_sample_prikey(&mut private_q_inv, &private_monomial_seeds[i]);
//        /// NOTE: the inversion is not implemented in constant time.
//        monomial_inv(&mut private_q, &private_q_inv);
//
//        let mut result_g = Matrix::<K, N>::default();
//        generator_monomial_mul(&mut result_g, &tmp_full_g, &private_q);
//
//        generator_RREF_pivot_reuse(&mut result_g, &mut is_pivot_column, &was_pivot_column, K);
//        /* note that the result is stored at i-1 as the first
//         * public key element is just a seed */
//        compress_rref(&mut pk.sf_g[i], &result_g, &is_pivot_column);
//    }
//} /* end LESS_keygen */