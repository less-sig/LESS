use std::ops::{Add, Mul, Sub, Div};
use std::ops::Deref;

use crate::helper::compute_ct_mask;

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Fq(pub u8);

const FQ127_INV_TABLE: [u8; 128] = [0, 1, 64, 85, 32, 51, 106, 109, 16, 113, 89, 104, 53, 88, 118, 17, 8, 15, 120, 107, 108, 121, 52, 116, 90, 61, 44, 80, 59, 92, 72, 41, 4, 77, 71, 98, 60, 103, 117, 114, 54, 31, 124, 65, 26, 48, 58, 100, 45, 70, 94, 5, 22, 12, 40, 97, 93, 78, 46, 28, 36, 25, 84, 125, 2, 43, 102, 91, 99, 81, 49, 34, 30, 87, 115, 105, 122, 33, 57, 82, 27, 69, 79, 101, 62, 3, 96, 73, 13, 10, 24, 67, 29, 56, 50, 123, 86, 55, 35, 68, 47, 83, 66, 37, 11, 75, 6, 19, 20, 7, 112, 119, 110, 9, 39, 74, 23, 38, 14, 111, 18, 21, 76, 95, 42, 63, 126, 0];


impl Fq {
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    /// a mod q, a < 256
    #[must_use]
    #[inline]
    const fn conditional_sub_raw(a: u8) -> u8 {
        let sub_q = a.overflowing_sub(Fq::Q).0;
        let mask = (-((sub_q >> Fq::Q_BITS) as i8)) as u8;
        let ret = (mask & Fq::Q).overflowing_add(sub_q).0;
        return ret;
    }

    /// constant time reduction
    /// a mod q,  a < 127**2
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::red(0);
    /// assert_eq!(result.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `a`: 16bit number to reduce
    ///
    /// # Returns
    /// a mod 127
    #[must_use]
    #[inline]
    pub const fn red(a: u16) -> Fq {
        return Fq(Fq::conditional_sub_raw(((a >> Fq::Q_BITS) as u8) + ((a as u8) & Fq::Q)));
    }

    /// a+b mod q,  a,b < 127
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::add(Fq(0), Fq(0));
    /// assert_eq!(result.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `a`: first addend
    /// - `b`: second addend
    ///
    /// # Returns
    /// a mod 127
    #[must_use]
    #[inline]
    pub const fn add(a: Fq, b: Fq) -> Fq {
        assert!(a.0 < 127); assert!(b.0 < 127);
        return Fq(Fq::conditional_sub_raw(a.0 + b.0));
    }
    
    /// a-b mod q,  a,b < 127
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::sub(Fq(0), Fq(0));
    /// assert_eq!(result.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `a`: minuend
    /// - `b`: subtrahend
    ///
    /// # Returns
    /// a-b mod 127
    #[must_use]
    #[inline]
    pub const fn sub(a: Fq, b: Fq) -> Fq {
        assert!(a.0 < 127); assert!(b.0 < 127);
        return Fq(Fq::conditional_sub_raw(a.0 + Fq::Q - b.0));
    }

    /// a*b mod q,  a,b < 127
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::mul(Fq(0), Fq(17));
    /// assert_eq!(result.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `a`: multiplier
    /// - `b`: multiplicant
    ///
    /// # Returns
    /// a*b mod 127
    #[must_use]
    #[inline]
    pub const fn mul(a: Fq, b: Fq) -> Fq {
        assert!(a.0 < 127); assert!(b.0 < 127);
        return Fq::red((a.0 as u16) * (b.0 as u16));
    }

    /// a^{-1} mod q, a < 127
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::inv(Fq(1));
    /// assert_eq!(result.0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: value to invert
    ///
    /// # Returns
    /// a^{-1} mod 127
    #[must_use]
    #[inline]
    pub fn inv(a: Fq) -> Fq {
        assert!(a.0 < 127);
        let j = a.0;
        let mut r: u8 = 0;
        for i in 1..127u8{
            compute_ct_mask(i, j);
            r &= FQ127_INV_TABLE[a.0 as usize];
        }
        Fq(r)
    }

    /// a^{-1} mod q, a < 127
    /// # Examples
    ///
    /// ```
    /// use less::fq::Fq;
    /// let result = Fq::inv_non_ct(Fq(1));
    /// assert_eq!(result.0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: value to invert
    ///
    /// # Returns
    /// a^{-1} mod 127
    #[must_use]
    #[inline]
    pub const fn inv_non_ct(a: Fq) -> Fq {
        assert!(a.0 < 127);
        Fq(FQ127_INV_TABLE[a.0 as usize])
    }
}

impl Add for Fq {
    type Output = Fq;
    #[inline]
    fn add(self, r: Fq) -> Fq {
        Fq::add(self, r)
    }
}

impl Sub for Fq {
    type Output = Fq;
    #[inline]
    fn sub(self, r: Fq) -> Fq {
        Fq::sub(self, r)
    }
}

impl Mul for Fq {
    type Output = Fq;
    #[inline]
    fn mul(self, r: Fq) -> Fq {
        Fq::mul(self, r)
    }
}

impl Div for Fq {
    type Output = Fq;
    #[inline]
    fn div(self, r: Fq) -> Fq {
        Fq::mul(self, Fq::inv(r))
    }
}

impl Deref for Fq {
    type Target = u8;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conditional_sub_raw() {
        for i in 0..127 {
            let t: u8 = Fq::conditional_sub_raw(i);
            assert!(t == i);
        }
        for i in 127..254 {
            let t: u8 = Fq::conditional_sub_raw(i);
            assert!(t == i - 127);
        }
    }

    #[test]
    fn add() {
        for i in 0..127 {
            for j in 0..127 {
                let a = Fq(i);
                let b = Fq(j);
                let t: Fq = Fq::add(a, b);
                assert!(t.0 == (i+j)%127);
            }
        }
    }

    #[test]
    fn sub() {
        for i in 0..127u16 {
            for j in 0..127u16 {
                let a = Fq(i as u8);
                let b = Fq(j as u8);
                let t: Fq = Fq::sub(a, b);
                assert!(t.0 == (((i + 127) - j)%127) as u8);
            }
        }
    }

    #[test]
    fn mul() {
        for i in 0..127u16 {
            for j in 0..127u16 {
                let a = Fq(i as u8);
                let b = Fq(j as u8);
                let t: Fq = Fq::mul(a, b);
                assert!(t.0 == ((i * j)%127) as u8);
            }
        }
    }

    #[test]
    fn inv() {
        for i in 1..127u16 {
            let a = Fq(i as u8);
            let b = Fq::inv(a);
            let t: Fq = Fq::mul(a, b);
            assert!(t.0 == 1);
        }
    }
}
