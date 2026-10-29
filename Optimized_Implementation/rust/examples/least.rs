use futures::{SinkExt, StreamExt};
use sha3::Digest;
use tokio::net::{TcpListener, TcpStream};
use tokio_util::bytes::Bytes;
use tokio_util::codec::{Framed, LengthDelimitedCodec};

use less::config::{
    Hash, Xof, K, N, N_K_pad, N_pad, T, HASH_DIGEST_LENGTH, PRIVATE_KEY_SEED_LENGTH_BYTES,
    SEED_LENGTH_BYTES,
};
use less::matrix::Matrix;
use less::monomial::Monomial;
use less::prng::{csprng_initialize, randombytes};
use less::utils::{sample_challenge, verify};

type MatrixFull = Matrix<K, N, N_pad>;
type MonomialFull = Monomial<N, N_pad>;
/// a generator matrix in RREF together with its pivot columns
type Rref = (MatrixFull, [u8; N]);

///sigma = (cmt, salt, {mu^(j)}_{j=1..s}), with s = T
struct Signature {
    cmt: [u8; HASH_DIGEST_LENGTH],
    salt: [u8; SEED_LENGTH_BYTES],
    mu: Vec<MonomialFull>,
}

impl Signature {
    /// the signature `least_sign` returns on failure, never verifies
    fn empty() -> Self {
        Signature { cmt: [0u8; HASH_DIGEST_LENGTH], salt: [0u8; SEED_LENGTH_BYTES], mu: Vec::new() }
    }
}

/// one party's end of the full mesh, `peers[i]` talks to party `i`,
/// `peers[id]` is `None`
struct Net {
    id: usize,
    peers: Vec<Option<Framed<TcpStream, LengthDelimitedCodec>>>,
}

impl Net {
    /// One step of the broadcast channel: party `from` evaluates `make_fn` and
    /// sends the result as one frame to every peer, every other party reads
    /// one frame from `from`. `make_fn` is only evaluated by the active party.
    ///
    /// # Parameters
    /// - `from`: the publishing party
    /// - `make_fn`: function producing the message, 
    ///         only called if `from == self.id`
    ///
    /// # Returns
    /// the published message
    async fn publish(&mut self, from: usize, make_fn: impl FnOnce() -> Vec<u8>) -> Vec<u8> {
        if from != self.id {
            let p = self.peers[from].as_mut().unwrap();
            // this blocks
            return p.next().await.unwrap().unwrap().to_vec();
        }
        let msg = make_fn();
        for p in self.peers.iter_mut().flatten() {
            p.send(Bytes::copy_from_slice(&msg)).await.unwrap();
        }
        msg
    }
}

/// Connects `n` parties pairwise over localhost, one TCP connection per pair.
///
/// # Parameters
/// - `n`: number of parties
///
/// # Returns
/// one `Net` per party, `nets[i].id == i`
async fn mesh(n: usize) -> Vec<Net> {
    let mut nets: Vec<Net> = (0..n)
        .map(|id| Net { id, peers: (0..n).map(|_| None).collect() })
        .collect();
    for i in 0..n {
        for j in i + 1..n {
            let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
            let a = TcpStream::connect(listener.local_addr().unwrap()).await.unwrap();
            let (b, _) = listener.accept().await.unwrap();
            nets[i].peers[j] = Some(Framed::new(a, LengthDelimitedCodec::new()));
            nets[j].peers[i] = Some(Framed::new(b, LengthDelimitedCodec::new()));
        }
    }
    nets
}

/// TODO move into Matrix
/// SampleGenerator(seed_pk): G_0 in systematic form, as in `less_sign`.
///
/// # Parameters
/// - `seed`: seed_pk
///
/// # Returns
/// G_0 and its pivot columns (the first K)
fn sample_generator(seed: &[u8]) -> Rref {
    let rref = Matrix::<K, { N - K }, N_K_pad>::rand_from_seed::<Xof>(seed);
    let column_pos: [u16; N_K_pad] = core::array::from_fn(|i| (i + K) as u16);
    (MatrixFull::from_compressed_rref(&column_pos, &rref), core::array::from_fn(|i| (i < K) as u8))
}

/// Turn: the parties sorted by the sha3 hash of their generator matrix.
///
/// # Parameters
/// - `pks`: `pks[i]` is the serialized generator matrix of party `i`
///
/// # Returns
/// the party ids in turn order, `turn[i]` is l_{i+1} of the pdf
fn turn_order(pks: &[Vec<u8>]) -> Vec<usize> {
    let mut turn: Vec<usize> = (0..pks.len()).collect();
    turn.sort_by_key(|&i| Hash::digest(&pks[i]).to_vec());
    turn
}

/// Key setup, run by every party: everyone publishes pk_i = RREF(G_0 * Q_i),
/// this fixes the turn order, along which G_i = RREF(G_{i-1} * Q_{l_i}) is
/// built. The party's keypair is a LESS keypair.
///
/// # Parameters
/// - `net`: this party's network
/// - `g0`: G_0
/// - `q_inv`: this party's secret key Q^-1
///
/// # Returns
/// the turn order and G_N
async fn setup(net: &mut Net, g0: &Rref, q_inv: &MonomialFull) -> (Vec<usize>, Rref) {
    let mut q = MonomialFull::zero();
    Monomial::inv_non_ct(&mut q, q_inv);

    let mut pks = Vec::new();
    for i in 0..net.peers.len() {
        pks.push(net.publish(i, || { let (g, p) = g0.0.rref_mul(&q).unwrap(); g.to_bytes(&p) }).await);
    }
    let turn = turn_order(&pks);

    let mut g = g0.clone();
    for &l in &turn {
        let b = net.publish(l, || { let (g, p) = g.0.rref_mul(&q).unwrap(); g.to_bytes(&p) }).await;
        g = MatrixFull::from_bytes(&b);
    }
    (turn, g)
}

/// Threshold sign, Alg. 1 and 2, run by every party. Every party holds
/// exactly one turn (t = N). A wrong commitment opening or a missing
/// canonical form yields `Signature::empty()`.
///
/// # Parameters
/// - `net`: this party's network
/// - `g0`: G_0
/// - `q_inv`: this party's secret key Q^-1
/// - `turn`: the turn order from `setup`
/// - `msg`: message to sign
///
/// # Returns
/// the signature, identical on every party
async fn least_sign(net: &mut Net, g0: &Rref, q_inv: &MonomialFull, turn: &[usize], msg: &[u8]) -> Signature {
    // NOTE: n = N from the spec, but N is also the code length
    let n = turn.len();

    // lines 2-3
    let mut salt_i_and_nonce_r = [0u8; 2 * SEED_LENGTH_BYTES];
    // TODO sample deterministically
    randombytes(&mut salt_i_and_nonce_r);
    let mut c = Vec::new();
    for i in 0..n {
        c.push(net.publish(i, || Hash::digest(&salt_i_and_nonce_r).to_vec()).await);
    }

    // lines 4-9, the ephemeral monomials come from a fresh seed
    let mut seed = [0u8; SEED_LENGTH_BYTES];
    // TODO sample deterministically
    randombytes(&mut seed);
    let mut state = csprng_initialize::<Xof>(&seed);
    let mut mu_tilde: Vec<MonomialFull> = Vec::new();
    let mut g: Vec<Rref> = (0..T).map(|_| g0.clone()).collect();
    for i in 0..n {
        // NOTE: the pdf states `s`, but T is the number of rounds.
        for j in 0..T {
            let b = net.publish(turn[i], || {
                // line 7
                let mu_i = MonomialFull::rand(&mut state);
                // line 8
                let (g_i, p) = g[j].0.rref_mul(&mu_i).unwrap();
                mu_tilde.push(mu_i);
                // line 9
                g_i.to_bytes(&p)
            }).await;
            g[j] = MatrixFull::from_bytes(&b);
        }
    }

    // Part/Algorithm 2
    // line 1, hashed straight away for line 4
    let mut hasher = Hash::new();
    for j in 0..T {
        // if any of the CF fails we quit completly.
        if !g[j].0.cf_hash(&g[j].1, &mut hasher) {
            return Signature::empty();
        }
    }

    // lines 2-3
    let mut salt = [0u8; SEED_LENGTH_BYTES];
    for i in 0..n {
        let open = net.publish(i, || salt_i_and_nonce_r.to_vec()).await;
        // all parties verify c_i
        if Hash::digest(&open).to_vec() != c[i] {
            return Signature::empty();
        }
        for k in 0..SEED_LENGTH_BYTES {
            salt[k] ^= open[k];
        }
    }

    // line 4
    hasher.update(&salt);
    hasher.update(msg);
    let cmt: [u8; HASH_DIGEST_LENGTH] = *hasher.finalize().as_array().unwrap();

    // line 5
    let mut b = [0u8; T];
    sample_challenge(&mut b, &cmt);
    // println!("{:?}", b);

    // lines 6-14
    // b_j = 0: every party reveals its mu~_j
    // b_j = 1: round-robin, mu^_i = sk^-1 * mu^_{i-1} * mu~_i (line 12)
    let mut mu = Vec::new();
    for j in 0..T {
        // line 9
        let mut mu_hat = MonomialFull::default();
        for i in 0..n {
            let m = net.publish(turn[i], || {
                if b[j] == 0 {
                    return mu_tilde[j].to_bytes();
                }
                let mut r = q_inv.clone();
                Monomial::mul_non_ct(&mut r, &mu_hat);
                Monomial::mul_non_ct(&mut r, &mu_tilde[j]);
                // line 13
                r.to_bytes()
            }).await;
            let m = MonomialFull::from_bytes(&m);
            if b[j] == 0 {
                // line 7: mu^(j) = mu~_1 * ... * mu~_N, in turn order
                Monomial::mul_non_ct(&mut mu_hat, &m);
            } else {
                mu_hat = m;
            }
        }
        // line 14
        mu.push(mu_hat);
    }
    // line 15
    Signature { cmt, salt, mu }  
}

/// Threshold verify, Alg. 3. Malformed signatures (wrong number of
/// responses, non invertible monomials, no canonical form) are rejected
/// instead of panicking.
/// NOTE: non-constant time
///
/// # Parameters
/// - `g0`: G_0
/// - `gn`: G_N, the public key of the group
/// - `msg`: signed message
/// - `sig`: signature
///
/// # Returns
/// true if the signature is valid
fn least_verify(g0: &Rref, gn: &Rref, msg: &[u8], sig: &Signature) -> bool {
    // check correct formatting:w
    if sig.mu.len() != T || !sig.mu.iter().all(|m| m.is_valid()) {
        return false;
    }

    let mut b = [0u8; T];
    sample_challenge(&mut b, &sig.cmt);

    let mut hasher = Hash::new();
    for j in 0..T {
        // lines 1-4
        let g = if b[j] == 0 { g0 } else { gn };
        let Some((gj, p)) = g.0.rref_mul(&sig.mu[j]) else {
            return false;
        };
        if !gj.cf_hash(&p, &mut hasher) {
            return false;
        }
    }
    hasher.update(&sig.salt);
    hasher.update(msg);

    // lines 5-8
    verify(&hasher.finalize(), &sig.cmt) == 0
}

/// Runs `n` parties in one process: keygen, `setup`, `least_sign` and
/// `least_verify` on every party, `n` is the first argument (default 3).
#[tokio::main]
async fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let msg = b"LEAST threshold signature";

    // seed_pk is common, every party draws its own sk
    let mut seed_pk = [0u8; SEED_LENGTH_BYTES];
    randombytes(&mut seed_pk);

    // every party is a tokio thread. Which itself is async, but in `publish`
    // every one is blocked except the active party
    let mut parties = Vec::new();
    for mut net in mesh(n).await {
        parties.push(tokio::spawn(async move {
            // basically keygen:
            let g0 = sample_generator(&seed_pk);
            let mut sk = [0u8; PRIVATE_KEY_SEED_LENGTH_BYTES];
            randombytes(&mut sk);

            // sk^-1, as in `less_keygen`
            let q_inv = MonomialFull::rand_from_seed::<Xof>(&sk);

            let (turn, gn) = setup(&mut net, &g0, &q_inv).await;
            let sig = least_sign(&mut net, &g0, &q_inv, &turn, msg).await;
            (net.id, turn, least_verify(&g0, &gn, msg, &sig))
        }));
    }
    for p in parties {
        let (id, turn, ok) = p.await.unwrap();
        println!("party {id}: turn order {turn:?}, signature valid: {ok}");
    }
}

/// run with `cargo test --release --example least`
#[cfg(test)]
mod tests {
    use super::*;

    /// party 1 runs `setup` with its key but signs with a key it made up,
    /// so every b_j = 1 response is wrong and nobody accepts the signature
    #[tokio::test(flavor = "multi_thread")]
    async fn cheating_party_is_rejected() {
        // number of parties
        let n = 3;
        let msg = b"LEAST threshold signature";
        let mut seed_pk = [0u8; SEED_LENGTH_BYTES];
        randombytes(&mut seed_pk);

        let cheating_party = 0;
        let mut parties = Vec::new();
        for mut net in mesh(n).await {
            parties.push(tokio::spawn(async move {
                let g0 = sample_generator(&seed_pk);
                let mut sk = [0u8; PRIVATE_KEY_SEED_LENGTH_BYTES];
                randombytes(&mut sk);
                let q_inv = MonomialFull::rand_from_seed::<Xof>(&sk);
                let (turn, gn) = setup(&mut net, &g0, &q_inv).await;

                let mut wrong = [0u8; PRIVATE_KEY_SEED_LENGTH_BYTES];
                randombytes(&mut wrong);
                let key = if net.id == cheating_party { MonomialFull::rand_from_seed::<Xof>(&wrong) } else { q_inv };

                let sig = least_sign(&mut net, &g0, &key, &turn, msg).await;
                least_verify(&g0, &gn, msg, &sig)
            }));
        }
        for p in parties {
            assert!(!p.await.unwrap(), "the signature of a cheater was accepted");
        }
    }
}
