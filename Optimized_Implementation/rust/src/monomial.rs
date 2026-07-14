use std::ops::{Index, IndexMut, Mul, MulAssign};
use sha3::digest::{ExtendableOutput, XofReader};
use crate::fq::Fq;
use crate::prng::{merge_exchange, rand_range_q_state_elements};
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
    /// let t = Monomial::<100>::new();
    /// assert_eq!(t.coeffs[0].0, 0);
    /// ```
    #[inline]
    pub fn new() -> Self {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }

        Monomial {
            coeffs: Vector::<N>::init(),
            perms: t,
        }
    }
    /// same as `new`
    #[inline]
    pub fn init() -> Self {
        Self::new()
    }
    /// creates a fully zero monomial matrix.
    /// NOTE: technically that's not a matrix.
    #[inline]
    pub fn zero() -> Self {
        let t: [u16; N] = [0u16; N];
        Monomial {
            coeffs: Vector::<N>::init(),
            perms: t,
        }
    }

    /// NOTE:
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

    /// sample a random matrix
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
    pub fn rand_from_seed<S>(seed: &[u8; 32]) -> Self
    where
        S: ExtendableOutput + Default + Clone
    {
        let mut hasher = S::default();
        hasher.update(seed.as_slice());
        let mut reader = hasher.finalize_xof();
        Self::rand(&mut reader)
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let A: Monomial<100> = Monomial::init();
    /// let B: Monomial<100> = Monomial::init();
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
    /// let result: Monomial<100> = Monomial::init();
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
    /// let a: Monomial<100> = Monomial::init();
    /// let mut b: Monomial<100> = Monomial::init();
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
    /// let a: Monomial<100> = Monomial::init();
    /// let mut b: Monomial<100> = Monomial::init();
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
    /// let t = Permutation::<100>::new();
    /// assert_eq!(t.perms[0], 0);
    /// ```
    pub fn new() -> Self {
        let mut t: [u16; N] = [0u16; N];
        for i in 0..N {
            t[i] = i as u16;
        }
        Permutation {
            perms: t,
        }
    }
    pub fn init() -> Self {
        Self::new()
    }

    /// # Examples
    ///
    /// ```
    /// use less::monomial::Monomial;
    /// let result: Monomial<100> = Monomial::init();
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
    /// let result: Permutation<100> = Permutation::init();
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
    /// let result: Permutation<100> = Permutation::init();
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
    /// let a: Permutation<100> = Permutation::init();
    /// let mut b: Permutation<100> = Permutation::init();
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
    /// let a: Permutation<100> = Permutation::init();
    /// let mut b: Permutation<100> = Permutation::init();
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
        let c = Monomial::<N>::new();
        for i in 0..c.dimension() {
            assert_eq!(c.coeffs[i].0, 0);
            assert_eq!(c.perms[i], i as u16);
        }
    }

    #[test]
    fn mul() {
        let a = Monomial::<N>::new();
        let mut b = Monomial::<N>::zero();
        Monomial::mul_non_ct(&mut b, &a);

        for i in 0..N {
            assert_eq!(b.coeffs[i].0, 0);
            assert_eq!(b.perms[i], 0u16);
        }

        b = Monomial::<N>::new();
        Monomial::mul_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn inv() {
        let a = Monomial::<N>::new();
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
        let c = Permutation::<N>::new();
        for i in 0..c.dimension() {
            assert_eq!(c.perms[i], i as u16);
        }
    }
    #[test]
    fn compress() {
        let c = Permutation::<N>::new();
        let d = c.compress((N/2) as u32);
        let e = Permutation::from_bytes(d);
        // NOTE only the first k assert_eq!(e, c);
    }

    #[test]
    fn mul() {
        let a = Permutation::<N>::new();
        let mut b = Permutation::<N>::new();
        Permutation::mul_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn inv() {
        let a = Permutation::<N>::new();
        let mut b = Permutation::<N>::new();
        Permutation::inv_non_ct(&mut b, &a);
        assert_eq!(b, a);
    }
}