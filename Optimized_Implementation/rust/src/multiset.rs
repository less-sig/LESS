// TODO remove
#![allow(dead_code)]

use std::ops::{Index, IndexMut};
use crate::vector::Vector;
use crate::matrix::MatrixNormalized;

/// Multiset is actually a histogram
/// NOTE: N: is most likely Q_PAD
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Multiset<const N:usize> (pub [u8; N]);


impl<const N:usize> Multiset<N> {
    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let result: Multiset<100> = Multiset::init();
    /// assert_eq!(result[0], 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Parameters
    /// - `dimension`: Size of the multiset
    ///
    /// # Returns
    /// A new multiset element
    #[inline]
    pub fn init() -> Self {
        Self([0; N])
    }

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let t = [0u8; 100];
    /// let result: Multiset<100> = Multiset::from_vector(&t);
    /// assert_eq!(result[0], 100);
    /// ```
    ///
    /// # Returns
    /// A new multiset element
    #[inline]
    pub fn from_vector(coefficients: &[u8; N]) -> Self {
        let mut ret = Self::init();
        for i in 0..N {
            ret[coefficients[i] as usize] += 1;
        }
        ret
    }

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// use less::vector::Vector;
    /// let t = Vector::<100>::from_u8(1);
    /// let mut c = Multiset::<100>::init();
    /// Multiset::<100>::from_row(&mut c, &t);
    /// assert_eq!(c[0], 0);
    /// assert_eq!(c[1], 100);
    /// ```
    ///
    /// # Returns
    /// A new multiset element
    #[inline]
    pub fn from_row<const M: usize>(out: &mut Multiset<N>, row: &Vector<M>) {
        Self::sort(out, row);
    }

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// use less::matrix::MatrixNormalized;
    /// let t = MatrixNormalized::<100>::from_u8(1);
    /// let c = Multiset::<100>::from_matrix(&t);
    /// ```
    ///
    /// # Parameters
    /// - `coefficients`: The input `Vec`.
    ///
    /// # Returns
    /// A new multiset element
    #[inline]
    pub fn from_matrix<const M: usize>(m: &MatrixNormalized<M>) -> [Self; M] {
        let mut tmp: [Multiset<N>; M] = core::array::from_fn(|_| Multiset::<N>::init() );
        for i in 0..N {
            Self::from_row::<M>(&mut tmp[i], &m.rows[i]);
        }
        tmp
    }

    /// row1[in]: pointer to the first row in histogram form
    /// row2[in]: pointer to the second row in histogram form
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let A: Multiset<100> = Multiset::init();
    /// let B: Multiset<100> = Multiset::init();
    /// A >= B;
    /// A <= B;
    /// A > B;
    /// A < B;
    /// ```
    ///
    /// # Returns
    ///  0 if multiset(row1) == multiset(row2)
    ///  x if multiset(row1) > multiset(row2)
    /// -x if multiset(row1) < multiset(row2)
    fn partial_cmp(a: &Self, b: &Self) -> i32 {
        let mut i: usize = 0;
        while i < N && a[i] == b[i] {
            i += 1;
        };

        (b[i] as i32) - (a[i]  as i32)
    }

    /// checks for equality
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let A: Multiset<100> = Multiset::init();
    /// let B: Multiset<100> = Multiset::init();
    /// A == B;
    /// A != B;
    /// ```
    fn eq(a: &Self, b: &Self) -> bool {
        Self::partial_cmp(a, b) == 0
    }

    /// Create an instance from a `Multiset`
    /// # Examples
    ///
    /// ```
    /// use less::multiset::Multiset;
    /// let result: Multiset<100> = Multiset::init();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    /// # Returns
    /// the dimension of a Multiset
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }

    /// sorts a row
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100> = MatrixNormalized::init();
    /// ```
    /// static internal function implementing the histogram function
    #[inline]
    fn sort<const M: usize>(out: &mut Multiset<N>, row: &Vector<M>) {
        // first clear the buffer
        for i in 0..N {
            out[i] = 0;
        }

        for i in 0..M {
            out[row[i].0 as usize] += 1;
        }
    }
}

impl<const N: usize> Index<usize> for Multiset<N> {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
impl<const N: usize> IndexMut<usize> for Multiset<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

/// TODO
#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 128;

    #[test]
    fn init() {
        let c = Multiset::<N>::init();
        for i in 0..N {
            assert_eq!(c[i], 0);
        }
    }

    #[test]
    fn from_vector() {
        let t = [1u8; N];
        let c = Multiset::<N>::from_vector(&t);
        assert_eq!(c[0], 0);
        assert_eq!(c[1], N as u8);
        for i in 2..N {
            assert_eq!(c[i], 0);
        }
    }

    #[test]
    fn from_row() {
        let t = Vector::<N>::from_u8(1);
        let mut c = Multiset::<N>::init();

        Multiset::<N>::from_row(&mut c, &t);
        assert_eq!(c[0], 0);
        assert_eq!(c[1], N as u8);
        for i in 2..N {
            assert_eq!(c[i], 0);
        }
    }
}
