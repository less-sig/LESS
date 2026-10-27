// TODO remove
#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(non_snake_case)]

use std::ops::{Index, IndexMut, Mul, MulAssign};
use sha3::digest::{ExtendableOutput, XofReader};
use crate::config::N8;
use crate::fq::Fq;
use crate::prng::{merge_exchange, rand_range_q_state_elements, randombytes};
use crate::vector::Vector;

/// N: actual size
/// N_pad: padded size
#[derive(Debug, PartialEq, Eq)]
pub struct Monomial<const N: usize, const N_PAD: usize> {
    pub coeffs: Vector<N, N_PAD>,
    pub perms: [u16; N],
}

#[derive(Debug, PartialEq, Eq)]
pub struct Permutation<const N: usize> {
    pub perms: [u16; N],
}

impl<const N: usize, const N_PAD: usize> Monomial<N, N_PAD> {
    const _ASSERT_PAD: () = assert!(N_PAD >= N, "N_PAD must be >= N");

    /// creates a new Monomial, initialized with as id
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// use less::fq::Fq;
    /// let t = Monomial::<100, 100>::default();
    /// assert_eq!(t.coeffs[0].0, 0);
    /// ```
    #[inline]
    pub fn default() -> Self {
        let () = Self::_ASSERT_PAD;

        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }

        Monomial {
            coeffs: Vector::<N, N_PAD>::default(),
            perms: t,
        }
    }

    /// NOTE: TODO doc test
    /// creates a fully zero monomial matrix.
    /// NOTE: technically that's not a matrix.
    #[inline]
    pub fn zero() -> Self {
        let t: [u16; N] = [0u16; N];
        Monomial {
            coeffs: Vector::<N, N_PAD>::default(),
            perms: t,
        }
    }


    /// NOTE: TODO doc test
    /// NOTE: a monomial matrix has to be invertible, so the coefficients are
    /// drawn from F_q^*.
    /// NOTE: constant time
    pub fn rand<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let coeffs = Vector::<N, N_PAD>::rand_star(state);

        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }
        merge_exchange(&mut t, state);

        Monomial {
            coeffs: coeffs,
            perms: t,
        }
    }

    /// NOTE: constant time
    /// sample a random monomial matrix. This means the Fq elements must be != 0
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// ```
    ///
    /// # Parameters
    /// - `seed`: 32/64 byts to seed the prng
    pub fn rand_from_seed<S>(seed: &[u8]) -> Self
    where
        S: ExtendableOutput + Default + Clone
    {
        // static_assert(seed.len() % 32 == 0);
        let mut hasher = S::default();
        hasher.update(seed);
        let mut state = hasher.finalize_xof();
        Self::rand(&mut state)
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let A: Monomial<100, 100> = Monomial::default();
    /// let B: Monomial<100, 100> = Monomial::default();
    /// A == B;
    /// A != B;
    /// ```
    ///
    /// # Returns
    ///  1 if a == b, else 0
    fn eq(a: &Self, b: &Self) -> bool {
        for i in 0..a.dimension() {
            if a.perms[i] != b.perms[i] {
                return false;
            }
            if a.coeffs[i] != b.coeffs[i] {
                return false;
            }
        }

        true
    }

    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let result: Monomial<100, 128> = Monomial::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// The dimension of the monomial
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }


    /// NOTE: non ct
    /// l = l*r
    ///
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let a: Monomial<100, 100> = Monomial::default();
    /// let mut b: Monomial<100, 100> = Monomial::default();
    /// Monomial::mul_non_ct(&mut b, &a);
    /// ```
    pub fn mul_non_ct(l: &mut Self, r: &Self) {
        for i in 0..N {
            // `perms[i]` is the destination of column i, so applying `l` and
            // then `r` sends column i to `r.perms[l.perms[i]]`.
            // NOTE: only index `i` of `l` is written and only the already
            // read index `t` is used, so this is safe in place.
            let t: usize = l.perms[i] as usize;
            l.perms[i] = r.perms[t];
            l.coeffs[i] = Fq::mul(l.coeffs[i], r.coeffs[t]);
        }
    }

    /// NOTE non ct
    /// l = r**-1
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let a: Monomial<100, 100> = Monomial::default();
    /// let mut b: Monomial<100, 100> = Monomial::default();
    /// Monomial::inv_non_ct(&mut b, &a);
    /// ```
    pub fn inv_non_ct(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t = r.perms[i] as usize;
            l.perms[t] = i as u16;
            l.coeffs[t] = Fq::inv(r.coeffs[i])
        } 
    }
}

impl <const N: usize, const N_PAD: usize> Mul for Monomial<N, N_PAD> {
    type Output = Monomial<N, N_PAD>;
    #[inline]
    fn mul(self, r: Monomial<N, N_PAD>) -> Monomial<N, N_PAD> {
        let mut ret = self;
        Monomial::<N, N_PAD>::mul_non_ct(&mut ret, &r);
        ret
    }
}
impl <const N: usize, const N_PAD: usize> MulAssign for Monomial<N, N_PAD> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        Monomial::<N, N_PAD>::mul_non_ct(self, &other);
    }
}



impl<const N: usize> Permutation<N> {
    /// creates a new Monomial, initialized with as id
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Permutation;
    /// use less::fq::Fq;
    /// let t = Permutation::<100>::default();
    /// assert_eq!(t.perms[0], 0);
    /// ```
    pub fn default() -> Self {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }
        Permutation {
            perms: t,
        }
    }

    /// translation of monomial_compose_action_information_set
    /// NOTE: constant time implementation
    /// NOTE: Only the permutation is computed, as this is the only thing we need
    /// since the adaption of canonical forms.
    /// \tparam M_PRIME[in]: number of elements in mu_tilde and is_pivot_cilumn
    /// \tparam M_PRIME_PAD[in]: padded number of elements
    /// \param pi_tilde[out]: mu_tilde^{-1} * information set
    /// \param mu_tilde[in]: input monomial matrix.
    /// \param is_pivot_column[in]: information set, where a 1 in the array
    ///         symbolizes that the corresponding columns is a pivot column.
    #[inline]
    pub fn from_information_set_composition
        <const M_PRIME: usize, const M_PRIME_PAD: usize>
        (mu_tilde: &Monomial<M_PRIME, M_PRIME_PAD>, is_pivot_column: &[u8; M_PRIME]) -> Self {
        let mut ret = Permutation::default();

        let mut piv_idx = 0usize;
        // the scan runs over all M_PRIME columns of the code, not over the
        // N columns of the information set
        for col_idx in 0..M_PRIME {
            let mut row_idx = 0u16;
            for t in 0..M_PRIME {
                // NOTE: do not break here.
                if mu_tilde.perms[t] == col_idx as u16 { // TODO COMPUTE_CT_MASK(j, pos);
                    row_idx = t as u16;
                }
            }

            // NOTE: this should not leak information.
            if is_pivot_column[col_idx] == 1 {
                ret.perms[piv_idx] = row_idx;
                piv_idx += 1;
            }
        }

        ret
    }

    /// NOTE: translation of 'monomial_compose_action'
    /// NOTE: constant time implementation
    /// composes a compactly stored action of a monomial on an IS with a regular
    /// monomial.
    /// NOTE: Only the permutation is computed, as this is the only thing we need
    /// since the adaption of canonical forms.
    /// \param out[out]: output monomial
    /// \param Q_in[in]: input monomial
    /// \param in[in]: input monomial
    #[inline]
    pub fn from_composition_action
        <const N_PRIME: usize, const N_PRIME_PAD: usize>(Q_in: &Monomial<N_PRIME, N_PRIME_PAD>, in_: &Permutation<N>) -> Self {
        let mut ret = Permutation::default();
        // to compose with monomial_action_IS_t, reverse the convention
        // for Q storage: store in permutation[i] the idx of the source column landing
        // as the i-th after the GQ product, and in coefficients[i] the coefficient
        // by which the column is multiplied upon landing
        let mut reverse_Q = Monomial::<N_PRIME, N_PRIME_PAD>::default();
        for i in 0..N_PRIME {
            let pos = Q_in.perms[i];
            // `reverse_Q` has N_PRIME entries, all of them have to be scanned
            for j in 0..N_PRIME {
                let mask: u16 = (-((j as u16 == pos) as i16)) as u16; // TODO COMPUTE_CT_MASK(j, pos);
                let not_mask: u16 = !mask;
                let value: u16 = (reverse_Q.perms[j] & not_mask) ^ ((i as u16) & mask);
                reverse_Q.perms[j] = value;
            }
        }

        for i in 0..N { 
            let pos: u16 = in_.perms[i];
            let mut value = 0u16;
            for j in 0..N_PRIME {
                let mask: u16 = (-((j as u16 == pos) as i16)) as u16; // TODO COMPUTE_CT_MASK(j, pos);
                value ^= reverse_Q.perms[j] & mask;
            }
            ret.perms[i] = value;
        }
        ret
    }

    /// TODO doc test
    /// compress a permutation to a byte string
    #[inline]
    pub fn bytes(&self) -> [u8; N8] {
        let mut ret = [0u8; N8];
        for i in 0..N {
            let limb: usize = (self.perms[i] / 8) as usize;
            let pos:  usize = (self.perms[i] % 8) as usize;
            ret[limb] ^= 1u8 << pos;
        }

        ret
    }

    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let result: Monomial<100, 100> = Monomial::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// The dimension of the monomial
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }

    /// checks if the given (N+7)/8 bytes are a valid canonical form action,
    /// e.g. this function checks if exactly `N` ones are in the input.
    /// \param b[in]: compressed canonical form output from `bytes`
    /// \return true: if the hamming weight is N
    ///         false: if the hamming weight is not N
    pub fn check_canonical_action(b: &[u8; N8]) -> bool {
        let mut w = 0u32;
        for i in 0..N8 {
            w += b[i].count_ones();
        }

        w == (N as u32)
    }

    /// translate form `CheckCanonicalAction`
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Permutation;
    /// let result: Permutation<100> = Permutation::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    pub fn from_bytes(b: [bool; N]) -> Self {
        let mut t: [u16; N] = [0u16; N];

        for i in 0..N {
            if b[i] == true {
                t[i] = i as u16;
            }
        }
        Permutation { perms: t }
    }

    /// inverse of `from_bytes`
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Permutation;
    /// let result: Permutation<100> = Permutation::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    pub fn compress(&self, k: u32) -> [bool; N] {
        let mut t = [false; N];
        for i in 0..k as usize {
            let p: usize = self.perms[i] as usize;
            t[p] = true;
        }

        t
    }

    /// NOTE: non ct
    /// l = l*r
    ///
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Permutation;
    /// let a: Permutation<100> = Permutation::default();
    /// let mut b: Permutation<100> = Permutation::default();
    /// Permutation::mul_non_ct(&mut b, &a);
    /// ```
    pub fn mul_non_ct(l: &mut Self, r: &Self) {
        for i in 0..N {
            // see `Monomial::mul_non_ct`
            let t: usize = l.perms[i] as usize;
            l.perms[i] = r.perms[t];
        }
    }

    /// NOTE non ct
    /// l = r**-1
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Permutation;
    /// let a: Permutation<100> = Permutation::default();
    /// let mut b: Permutation<100> = Permutation::default();
    /// Permutation::inv_non_ct(&mut b, &a);
    /// ```
    pub fn inv_non_ct(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t = r.perms[i] as usize;
            l.perms[t] = i as u16;
        } 
    }
}

impl<const N: usize> Index<usize> for Permutation<N> {
    type Output = u16;
    fn index(&self, index: usize) -> &Self::Output {
        &self.perms[index]
    }
}
impl<const N: usize> IndexMut<usize> for Permutation<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.perms[index]
    }
}

impl <const N: usize> Mul for Permutation<N> {
    type Output = Permutation<N>;
    #[inline]
    fn mul(self, r: Permutation<N>) -> Permutation<N> {
        let mut ret = self;
        Permutation::<N>::mul_non_ct(&mut ret, &r);
        ret
    }
}
impl <const N: usize> MulAssign for Permutation<N> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        Permutation::<N>::mul_non_ct(self, &other);
    }
}




#[cfg(test)]
mod tests_Monomial {
    use super::*;
    use crate::config::Xof;
    use crate::prng::csprng_initialize;
    const N: usize = 32;

    #[test]
    fn new() {
        let c = Monomial::<N, N>::default();
        for i in 0..c.dimension() {
            assert_eq!(c.coeffs[i].0, 0);
            assert_eq!(c.perms[i], i as u16);
        }
    }

    #[test]
    fn mul() {
        let a = Monomial::<N, N>::default();
        let mut b = Monomial::<N, N>::zero();
        Monomial::mul_non_ct(&mut b, &a);

        for i in 0..N {
            assert_eq!(b.coeffs[i].0, 0);
            assert_eq!(b.perms[i], 0u16);
        }

        b = Monomial::<N, N>::default();
        Monomial::mul_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn inv() {
        let mut s = csprng_initialize::<Xof>(&[1u8; 32]);
        let a = Monomial::<N, N>::rand(&mut s);
        let mut b = Monomial::<N, N>::rand(&mut s);
        let mut c = Monomial::<N, N>::rand(&mut s);
        Monomial::inv_non_ct(&mut b, &a);
        Monomial::inv_non_ct(&mut c, &b);
        assert_ne!(b, a);
        assert_eq!(c, a);
    }

    #[test]
    fn zero() {
        let c = Monomial::<N, N>::zero();
        for i in 0..N {
            assert_eq!(c.coeffs[i].0, 0);
            assert_eq!(c.perms[i], 0);
        }
    }

    #[test]
    fn dimension() {
        assert_eq!(Monomial::<N, N>::default().dimension(), N);
    }

    #[test]
    fn eq() {
        let a = Monomial::<N, N>::default();
        let b = Monomial::<N, N>::default();
        assert!(Monomial::eq(&a, &b));
        assert!(!Monomial::eq(&a, &Monomial::<N, N>::zero()));
    }

    /// a monomial matrix must be invertible, so no coefficient may be zero
    /// and the permutation must be a permutation
    #[test]
    fn rand_is_a_monomial() {
        for seed in 0..64u8 {
            let mut s = csprng_initialize::<Xof>(&[seed; 32]);
            let a = Monomial::<N, N>::rand(&mut s);

            let mut seen = [false; N];
            for i in 0..N {
                assert_ne!(a.coeffs[i].0, 0, "seed {}: coefficient {} is zero", seed, i);
                assert!(a.coeffs[i].0 < 127);
                assert!(!seen[a.perms[i] as usize], "permutation is not injective");
                seen[a.perms[i] as usize] = true;
            }
        }
    }

    #[test]
    fn rand_from_seed_is_a_monomial() {
        for seed in 0..64u8 {
            let a = Monomial::<N, N>::rand_from_seed::<Xof>(&[seed; 32]);

            let mut seen = [false; N];
            for i in 0..N {
                assert_ne!(a.coeffs[i].0, 0, "seed {}: coefficient {} is zero", seed, i);
                assert!(!seen[a.perms[i] as usize]);
                seen[a.perms[i] as usize] = true;
            }
        }

        let a = Monomial::<N, N>::rand_from_seed::<Xof>(&[3u8; 32]);
        let b = Monomial::<N, N>::rand_from_seed::<Xof>(&[3u8; 32]);
        assert_eq!(a, b);
        let c = Monomial::<N, N>::rand_from_seed::<Xof>(&[4u8; 32]);
        assert_ne!(a, c);
    }

    /// rand and rand_from_seed must request the randomness in the same order
    #[test]
    fn rand_matches_rand_from_seed() {
        let mut s = csprng_initialize::<Xof>(&[9u8; 32]);
        let a = Monomial::<N, N>::rand(&mut s);
        let b = Monomial::<N, N>::rand_from_seed::<Xof>(&[9u8; 32]);
        assert_eq!(a, b);
    }

    /// m * m^-1 must be the identity monomial
    #[test]
    fn mul_inv_is_identity() {
        let a = Monomial::<N, N>::rand_from_seed::<Xof>(&[5u8; 32]);
        let mut a_inv = Monomial::<N, N>::zero();
        Monomial::inv_non_ct(&mut a_inv, &a);

        let mut t = Monomial::<N, N>::zero();
        for i in 0..N {
            t.perms[i] = a.perms[i];
            t.coeffs[i] = a.coeffs[i];
        }
        Monomial::mul_non_ct(&mut t, &a_inv);

        for i in 0..N {
            assert_eq!(t.perms[i], i as u16, "permutation at {}", i);
            assert_eq!(t.coeffs[i].0, 1, "coefficient at {}", i);
        }
    }
}

#[cfg(test)]
mod tests_Permutation {
    use super::*;
    const N: usize = 32;

    #[test]
    fn new() {
        let c = Permutation::<N>::default();
        for i in 0..c.dimension() {
            assert_eq!(c.perms[i], i as u16);
        }
    }

    #[test]
    fn compress() {
        let c = Permutation::<N>::default();
        let d = c.compress((N/2) as u32);
        let e = Permutation::from_bytes(d);
        // compress marks the first k images, from_bytes recovers them
        for i in 0..N/2 {
            assert!(d[i]);
            assert_eq!(e.perms[i], i as u16);
        }
        for i in N/2..N {
            assert!(!d[i]);
        }
    }

    #[test]
    fn dimension() {
        assert_eq!(Permutation::<N>::default().dimension(), N);
    }

    #[test]
    fn index() {
        let mut c = Permutation::<N>::default();
        c[3] = 17;
        assert_eq!(c[3], 17);
    }

    /// CosetRep: exactly K bits are set, one per image of the permutation
    #[test]
    fn bytes() {
        use crate::config::{K, N8};
        let c = Permutation::<K>::default();
        let b = c.bytes();
        let weight: u32 = b.iter().map(|x| x.count_ones()).sum();
        assert_eq!(weight, K as u32);
        assert_eq!(b.len(), N8);
    }

    /// pi_tilde must select exactly the K pivot columns
    #[test]
    fn from_information_set_composition() {
        use crate::config::{K, N, N_pad, Xof};

        let mu_tilde = Monomial::<N, N_pad>::rand_from_seed::<Xof>(&[2u8; 32]);

        // the first K columns are the information set
        let mut is_pivot_column = [0u8; N];
        for i in 0..K {
            is_pivot_column[i] = 1;
        }

        let pi = Permutation::<K>::from_information_set_composition::<N, N_pad>(
            &mu_tilde, &is_pivot_column);

        // pi.perms[j] is the source column that mu_tilde maps onto pivot j
        let mut seen = [false; N];
        for j in 0..K {
            let src = pi.perms[j] as usize;
            assert!(src < N);
            assert!(!seen[src], "pi_tilde is not injective");
            seen[src] = true;
            assert_eq!(mu_tilde.perms[src] as usize, j);
        }
    }

    /// composing with the identity monomial must leave the action unchanged
    #[test]
    fn from_composition_action_identity() {
        use crate::config::{K, N, N_pad};

        let id = Monomial::<N, N_pad>::default();
        let input = Permutation::<K>::default();
        let out = Permutation::<K>::from_composition_action::<N, N_pad>(&id, &input);
        assert_eq!(out.perms, input.perms);
    }

    #[test]
    fn mul() {
        let a = Permutation::<N>::default();
        let mut b = Permutation::<N>::default();
        Permutation::mul_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn inv() {
        let a = Permutation::<N>::default();
        let mut b = Permutation::<N>::default();
        Permutation::inv_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }
}
