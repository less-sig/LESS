#![allow(dead_code)]

use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Deref;
use std::ops::DerefMut;

use crate::fq::Fq;


#[derive(Debug)]
pub struct Vector<const N: usize>(pub [Fq; N]);

impl <const N: usize> Vector<N>{
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    /// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 100;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(0); N]);
    /// let b = Vector::<N>([Fq(1); N]);
    /// Vector::add(&mut c, &a, &b);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn add(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        if is_x86_feature_detected!("avx2") {
            assert!(N%16 == 0);
            Self::add_avx2(c, a, b);
        } else {
            for i in 0..N {
                c[i] = Fq::add(a[i], b[i]);
            }
        }
    }

    #[inline]
    pub fn add_avx2(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        let n = N / 16;
        unsafe {
            // for i in (0..N).step_by(16) {
            //     const _
            // }
        }
    }

    /// c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 100;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(1); N]);
    /// let b = Vector::<N>([Fq(0); N]);
    /// Vector::sub(&mut c, &a, &b);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn sub(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        for i in 0..N {
            c[i] = Fq::sub(a[i], b[i]);
        }
    }

    /// c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 100;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(0); N]);
    /// let b = Vector::<N>([Fq(1); N]);
    /// Vector::mul(&mut c, &a, &b);
    /// sert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn mul(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        for i in 0..N {
            c[i] = Fq::mul(a[i], b[i]);
        }
    }
}

impl<const N: usize> Index<usize> for Vector<N> {
    type Output = Fq;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const N: usize> IndexMut<usize> for Vector<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const N: usize> Deref for Vector<N> {
    type Target = [Fq];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

//impl<const N: usize> Deref for Vector<N> {
//    type Target = [Fq; N];
//
//    fn deref(&self) -> &Self::Target {
//        &self.0
//    }
//}
impl<const N: usize> DerefMut for Vector<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add() {
        const N: usize = 128;
        let mut c = Vector::<N>([Fq(0); N]);
        let a = Vector::<N>([Fq(0); N]);
        let b = Vector::<N>([Fq(1); N]);
        Vector::add(&mut c, &a, &b);
        for i in 0..N {
            assert!(c[i].0 == 1);
        }
    }

    #[test]
    fn sub() {
        const N: usize = 128;
        let mut c = Vector::<N>([Fq(0); N]);
        let a = Vector::<N>([Fq(1); N]);
        let b = Vector::<N>([Fq(0); N]);
        Vector::sub(&mut c, &a, &b);
        for i in 0..N {
            assert!(c[i].0 == 1);
        }
    }

    #[test]
    fn mul() {
        const N: usize = 128;
        let mut c = Vector::<N>([Fq(0); N]);
        let a = Vector::<N>([Fq(0); N]);
        let b = Vector::<N>([Fq(1); N]);
        Vector::mul(&mut c, &a, &b);
        for i in 0..N {
            assert!(c[i].0 == 0);
        }
    }
}
