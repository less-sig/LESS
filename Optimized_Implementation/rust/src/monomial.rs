use std::ops::{Index, IndexMut};
use sha3::digest::ExtendableOutput;
use crate::fq::Fq;
use crate::prng::rand_range_q_state_elements;
use crate::vector::Vector;

#[derive(Debug, PartialEq, Eq)]
pub struct Monomial<const N: usize> {
    pub coeffs: Vector<N>,
    pub perms: [u16; N],
}

#[derive(Debug, PartialEq, Eq)]
pub struct MonomialInformationSet<const N: usize> {
    pub perms: [u16; N],
}

impl<const N: usize> Monomial<N> {
    ///
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
        let mut t: [u16; N] = [0u16; N];
        Monomial {
            coeffs: Vector::<N>::init(),
            perms: t,
        }
    }

    /// NOTE: only for testing. Not really useful.
    pub fn rand() -> Self {
        let mut t: [u16; N] = [0u16; N];// TODO use merge exchange
        Monomial {
            coeffs: Vector::<N>::rand(),
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
        let mut a = Self::init();
        // TODO let mut hasher = Shake128::default();
        a
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

    /// Create an instance from a `Vector`
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let result: Multiset<100> = Multiset::init();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// The dimension of the vector
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }


    /// l = l*r
    pub fn mul(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t: usize = r.perms[i] as usize;
            l.perms[i] = l.perms[t] ;
            l.coeffs[i] = Fq::mul(l.coeffs[t], r.coeffs[i]);
        }
    }

    /// NOTE non ct
    /// l = r**-1
    pub fn inv(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t = r.perms[i] as usize;
            l.perms[t] = i as u16;
            l.coeffs[t] = Fq::inv(r.coeffs[i])
        } 
    }
}


impl<const N: usize> MonomialInformationSet<N> {
    /// initialize with zero
    pub fn new() -> Self {
        MonomialInformationSet {
            perms: [0u16; N],
        }
    }

    /// translate form `CheckCanonicakAction`
    pub fn from_bytes(b: [bool; N]) -> Self {
        let mut t: [u16; N] = [0u16; N];

        for i in 0..N {
            if b[i] == true {
                t[i] = i as u16;
            }
        }
        MonomialInformationSet { perms: t }
    }

    ///
    pub fn compress(&self, k: u32) -> [bool; N] {
        let mut t = [false; N];
        for i in 0..k as usize {
            let p: usize = self.perms[i] as usize;
            t[p] = true;
        }

        t
    }

    /// NOTE non ct
    /// l = r**-1
    pub fn inv(l: &mut Self, r: &Self) {
        for i in 0..N {
            let t = r.perms[i] as usize;
            l.perms[t] = i as u16;
        } 
    }
}

impl<const N: usize> Index<usize> for MonomialInformationSet<N> {
    type Output = u16;
    fn index(&self, index: usize) -> &Self::Output {
        &self.perms[index]
    }
}
impl<const N: usize> IndexMut<usize> for MonomialInformationSet<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.perms[index]
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
}