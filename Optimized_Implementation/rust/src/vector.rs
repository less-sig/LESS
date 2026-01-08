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

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use crate::vector::Vector;
    /// use crate::fq::Fq;
    /// let t = core::array::from_fn(|_| Fq::default());
    /// let result: Vector<100> = Vector::from_vector(t);
    /// assert_eq!(result[0], 0);
    /// ```
    ///
    /// # Parameters
    /// - `coefficients`: The first number.
    ///
    /// # Returns
    /// A new vector element
    #[inline]
    pub fn from_vector(coefficients: [Fq; N]) -> Self {
        Self { 0: coefficients }
    }

    /// same as ::init()
    #[inline]
    pub fn from_u8(val: u8) -> Self {
        let mut ret = Self::new();
        for i in 0..N {
            ret.0[i] = Fq(val);
        }

        return ret;
    }

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use crate::vector::Vector;
    /// use crate::fq::Fq;
    /// let result: Vector<100> = Vector::init();
    /// assert_eq!(result[0], 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Parameters
    /// - `dimension`: Size of the vector
    ///
    /// # Returns
    /// A new vector element
    #[inline]
    pub fn init() -> Self {
        Self {
            0: core::array::from_fn(|_| Fq::default() ),
        }
    }
   
    /// same as ::init()
    #[inline]
    pub fn zero(&mut self) {
        for i in 0..N {
            self.0[i] = Fq(0);
        }
    }

    pub fn new() -> Vector<N> {
        Vector {
            0: [Fq(0); N],
        }
    }

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

    /// TODO not finished
    #[inline]
    pub fn add_avx2(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        //let n = N / 16;
        //unsafe {
        //    // for i in (0..N).step_by(16) {
        //    //     const _
        //    // }
        //}
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

    /// NOTE: non ct
    /// c[i] = a[i]^{-1} mod q, a[i] < 127 for i in 0..N
    /// # Examples
    /// ```
    /// const N: usize = 128;
    /// let mut row_out = Vector::<N>::new();
    /// let row1 = Vector::<N>::from_u8(1);
    /// gf127_row_inv(&mut row_out, &row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to accumulate
    pub fn inv(row_out: &mut Vector<N>,
               in1: &Vector<N>) {
        for col in 0..N {
            row_out[col] = Fq::inv_non_ct(in1[col]);
        }
    }
    
    /// c[i] == c[j] for all i < j < n
    /// # Examples
    /// ```
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::all_same(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to accumulate
    ///
    /// # return 
    ///     1 if all elements are the same
    ///     0 else
    pub fn all_same(row: &Vector<N>) -> bool {
        let v = row[0];
        for col in 1..N {
            if row[col] != v {
                return false;
            }
        }
        return true;
    }

    /// c[i] == 0 for all i  < n
    /// # Examples
    /// ```
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::contains_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row 
    ///
    /// # return 
    ///     1 if a zero is in the row
    ///     0 else
    pub fn contain_zero(row: &Vector<N>) -> bool {
        for col in 0..N {
            if row[col] == Fq(0) {
                return true;
            }
        }
        return false;
    }

    /// c[i] == 0 for all i  < n
    /// # Examples
    /// ```
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::count_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row 
    ///
    /// # return 
    ///     number of zeros in row
    pub fn count_zero(row: &Vector<N>) -> u32 {
        let mut c: u32 = 0;
        for col in 0..N {
            if row[col] == Fq(0) {
                c += 1;
            }
        }
        return c;
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
