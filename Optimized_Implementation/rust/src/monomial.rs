// TODO remove
#![allow(dead_code)]
#![allow(unused_imports)]

use std::ops::{Index, IndexMut, Mul, MulAssign};
use sha3::digest::{ExtendableOutput, XofReader};
use crate::config::N8;
use crate::fq::Fq;
use crate::prng::{merge_exchange, rand_range_q_state_elements, randombytes};
use crate::vector::Vector;

#[derive(Debug, PartialEq, Eq)]
pub struct Monomial<const N: usize> {
    pub coeffs: Vector<N>,
    pub perms: [u16; N],
}

#[derive(Debug, PartialEq, Eq)]
pub struct Permutation<const N: usize> {
    pub perms: [u16; N],
}

impl<const N: usize> Monomial<N> {
    /// creates a new Monomial, initialized with as id
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// use less::fq::Fq;
    /// let t = Monomial::<100>::default();
    /// assert_eq!(t.coeffs[0].0, 0);
    /// ```
    #[inline]
    pub fn default() -> Self {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }

        Monomial {
            coeffs: Vector::<N>::default(),
            perms: t,
        }
    }

    /// creates a fully zero monomial matrix.
    /// NOTE: technically that's not a matrix.
    #[inline]
    pub fn zero() -> Self {
        let t: [u16; N] = [0u16; N];
        Monomial {
            coeffs: Vector::<N>::default(),
            perms: t,
        }
    }


    /// NOTE: TODO
    pub fn rand<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }
        merge_exchange(&mut t, state);
        Monomial {
            coeffs: Vector::<N>::rand(state),
            perms: t,
        }
    }

    /// sample a random monomial matrix. This means the Fq elements must be != 0
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// ```
    ///
    /// # Parameters
    /// - `seed`: NxM matrix
    pub fn rand_from_seed<S>(seed: &[u8]) -> Self
    where
        S: ExtendableOutput + Default + Clone
    {
        // static_assert(seed.len() % 32 == 0);
        let mut hasher = S::default();
        hasher.update(seed);
        let mut state = hasher.finalize_xof();
        // Self::rand(&mut reader)


        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }
        let mut ret = Monomial {
            coeffs: Vector::<N>::rand(&mut state),
            perms: t,
        };

        merge_exchange(&mut ret.perms, &mut state);
        ret
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let A: Monomial<100> = Monomial::default();
    /// let B: Monomial<100> = Monomial::default();
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
    /// let result: Monomial<100> = Monomial::default();
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
    /// let a: Monomial<100> = Monomial::default();
    /// let mut b: Monomial<100> = Monomial::default();
    /// Monomial::mul_non_ct(&mut b, &a);
    /// ```
    pub fn mul_non_ct(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t: usize = r.perms[i] as usize;
            l.perms[i] = l.perms[t] ;
            l.coeffs[i] = Fq::mul(l.coeffs[t], r.coeffs[i]);
        }
    }

    /// NOTE non ct
    /// l = r**-1
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let a: Monomial<100> = Monomial::default();
    /// let mut b: Monomial<100> = Monomial::default();
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

impl <const N: usize> Mul for Monomial<N> {
    type Output = Monomial<N>;
    #[inline]
    fn mul(self, r: Monomial<N>) -> Monomial<N> {
        let mut ret = self;
        Monomial::<N>::mul_non_ct(&mut ret, &r);
        ret
    }
}
impl <const N: usize> MulAssign for Monomial<N> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        Monomial::<N>::mul_non_ct(self, &other);
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
    /// \param pi_tilde[out]: mu_tilde^{-1} * information set
    /// \param mu_tilde[in]: input monomial matrix.
    /// \param is_pivot_column[in]: information set, where a 1 in the array
    ///         symbolizes that the corresponding columns is a pivot column.
    #[inline]
    pub fn from_information_set_composition<const M_PRIME: usize, const N_PRIME: usize>(mu_tilde: &Monomial<M_PRIME>, is_pivot_column: &[u8; N_PRIME]) -> Self {
        let mut ret = Permutation::default();

        let mut piv_idx = 0usize;
        for col_idx in 0..N {
            let mut row_idx = 0u16;
            for t in 0..N {
                // NOTE: do not break here.
                // NOTE: the compiler should not be able to optimize a break here. As it doesnt know that each
                //      value in `permutation` is unique, hence it "could" be that there are multiple position == col_idx
                if mu_tilde.perms[t] == col_idx as u16 {
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
    pub fn from_composition_action<const N_PRIME: usize>(Q_in: &Monomial<N_PRIME>, in_: &Permutation<N>) -> Self {
        let mut ret = Permutation::default();
        // to compose with monomial_action_IS_t, reverse the convention
        // for Q storage: store in permutation[i] the idx of the source column landing
        // as the i-th after the GQ product, and in coefficients[i] the coefficient
        // by which the column is multiplied upon landing
        let mut reverse_Q = Monomial::<N_PRIME>::default();
        for i in 0..N_PRIME {
            // we want to compute:
            // reverse_Q.permutation[Q_in->permutation[i]] = i;

            // NOTE: the type `uint16_t` is rather important. Otherwise the compiler is needed to emit
            // costly zero-extend mov instructions. At least on a ryzen 5 7600X the resulting code
            // is ~2x slower.
            let pos = Q_in.perms[i];
            for j in 0..N {
                let mask: u16 = (-((j as u16 == pos) as i16)) as u16; // TODO COMPUTE_CT_MASK(j, pos);
                let not_mask: u16 = !mask;
                let value: u16 = (reverse_Q.perms[j] & not_mask) ^ ((i as u16) & mask);
                reverse_Q.perms[j] = value;
            }
        }

        for i in 0..N { 
            // we want to compute:
            // out->permutation[i] = reverse_Q.permutation[in->permutation[i]];
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
    /// let result: Monomial<100> = Monomial::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// The dimension of the monomial
    #[inline]
    pub fn dimension(&self) -> usize {
        N
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
            let t: usize = r.perms[i] as usize;
            l.perms[i] = l.perms[t] ;
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
mod tests {
    use super::*;
    const N: usize = 32;

    #[test]
    fn new() {
        let c = Monomial::<N>::default();
        for i in 0..c.dimension() {
            assert_eq!(c.coeffs[i].0, 0);
            assert_eq!(c.perms[i], i as u16);
        }
    }

    #[test]
    fn mul() {
        let a = Monomial::<N>::default();
        let mut b = Monomial::<N>::zero();
        Monomial::mul_non_ct(&mut b, &a);

        for i in 0..N {
            assert_eq!(b.coeffs[i].0, 0);
            assert_eq!(b.perms[i], 0u16);
        }

        b = Monomial::<N>::default();
        Monomial::mul_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn inv() {
        let a = Monomial::<N>::default();
        let mut b = Monomial::<N>::zero();
        Monomial::inv_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }
}

#[cfg(test)]
mod tests_MonomialInformationSet {
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
        // NOTE only the first k assert_eq!(e, c);
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