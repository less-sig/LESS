#![allow(dead_code)]
#![allow(unreachable_code)]

use sha3::digest::XofReader;
use std::ops:: {
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign,
    Index, IndexMut, Deref, DerefMut
};
use std::cmp::Ordering;
use std::{ fmt::Display, fmt::Formatter, fmt::Result };
use crate::constants::{
    Q_PAD,
};
use crate::fq::Fq;
use crate::multiset::Multiset;
use crate::prng::{
    rand_range_q_state_elements,fq_star_rnd_state_elements
};

#[cfg(target_arch = "x86_64")]
use crate::opt::{
    gf127_row_add_avx2, gf127_row_add2_avx2,
    gf127_row_sub_avx2, gf127_row_sub2_avx2,
    gf127_row_mul_avx2, gf127_row_mul2_avx2,
    gf127_row_scalar_mul_avx2, gf127_row_scalar_mul2_avx2,
    gf127_row_acc_avx2, gf127_row_acc_inv_avx2,
    gf127_row_contains_zero_avx2, gf127_row_count_zero_avx2,
};

#[cfg(target_arch = "x86_64")]
use crate::opt::{
    gf127_row_add_avx512, gf127_row_add2_avx512,
    gf127_row_sub_avx512, gf127_row_sub2_avx512,
    gf127_row_mul_avx512, gf127_row_mul2_avx512,
    gf127_row_inv_avx512,
    gf127_row_scalar_mul_avx512, gf127_row_scalar_mul2_avx512,
    gf127_row_acc_avx512, gf127_row_acc_inv_avx512,
    gf127_row_contains_zero_avx512,
};

#[cfg(target_arch = "aarch64")]
use crate::opt::{
    gf127_row_add_neon, gf127_row_add2_neon,
    gf127_row_sub_neon, gf127_row_sub2_neon,
    gf127_row_mul_neon, gf127_row_mul2_neon,
    gf127_row_inv_neon,
    gf127_row_scalar_mul_neon, gf127_row_scalar_mul2_neon,
    gf127_row_acc_neon, gf127_row_acc_inv_neon,
    gf127_row_contains_zero_neon, gf127_row_count_zero_avx2,
};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(align(32))] // NOTE alignment only for avx2/avx512
pub struct Vector<const N: usize>(pub [Fq; N]);

impl <const N: usize> Vector<N>{
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// let result: Vector<100> = Vector::default();
    /// assert_eq!(result[0].0, 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Parameters
    /// - `dimension`: Size of the vector
    ///
    /// # Returns
    /// A new vector element
    #[inline]
    pub fn default() -> Self {
        Self {
            0: core::array::from_fn(|_| Fq::default() ),
        }
        // both valid
        // Vector {
        //     0: [Fq(0); N],
        // }
    }

    /// sample a random vector
    pub fn rand<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut a = Self::default();
        rand_range_q_state_elements(&mut a.0, state);
        a
    }

    /// sample a random vector
    pub fn rand_star<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut a = Self::default();
        fq_star_rnd_state_elements(&mut a.0, state);
        a
    }

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// let t = core::array::from_fn(|_| Fq::default());
    /// let result: Vector<100> = Vector::from_vector(&t);
    /// assert_eq!(result[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `coefficients`: The input `Vec`.
    ///
    /// # Returns
    /// A new vector element
    #[inline]
    pub fn from_vector(coefficients: &[Fq; N]) -> Self {
        Self { 0: *coefficients }
    }

    /// same as ::init()
    #[inline]
    pub fn from_u8(val: u8) -> Self {
        let mut ret = Self::default();
        for i in 0..N {
            ret.0[i] = Fq(val);
        }
        ret
    }

    /// partially compares two matrices
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let A: Vector<100> = Vector::default();
    /// let B: Vector<100> = Vector::default();
    /// A >= B;
    /// A <= B;
    /// A > B;
    /// A < B;
    /// ```
    ///
    /// # Returns
    ///  0 if a == b
    ///  x if a  > b
    /// -x if a  < b
    fn partial_cmp(a: &Self, b: &Self) -> Option<Ordering>{
        let mut aa = Multiset::<Q_PAD>::default();
        let mut bb = Multiset::<Q_PAD>::default();
        Multiset::from_row(&mut aa, a);
        Multiset::from_row(&mut bb, b);
        Multiset::<Q_PAD>::partial_cmp(&aa, &bb)
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let A: Vector<100> = Vector::default();
    /// let B: Vector<100> = Vector::default();
    /// A == B;
    /// A != B;
    /// ```
    ///
    /// # Returns
    ///  1 if a == b, else 0
    fn eq(a: &Self, b: &Self) -> bool {
        Self::partial_cmp(a, b).unwrap().is_eq()
    }

    /// Create an instance from a `Vector`
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let result: Vector<100> = Vector::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// The dimension of the vector
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }

    /// c[i] = a[i]+b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
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
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_add_avx512(c, a, b);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert!(N % 32 == 0);
                unsafe {
                    gf127_row_add_avx2(c, a, b);
                }
            }
            return;
        }
        #[cfg(target_feature = "neon")]
        {
            gf127_row_add_neon(c, a, b);
            return;
        }

        for i in 0..N {
            c[i] = Fq::add(a[i], b[i]);
        }
    }

    /// c[i] = c[i]+a[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(1); N]);
    /// Vector::add2(&mut c, &a);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn add2(c: &mut Vector<N>, a: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_add2_avx512(c, a);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert!(N % 32 == 0);
                unsafe {
                    gf127_row_add2_avx2(c, a);
                }
            }
            return;
        }

        #[cfg(target_feature = "neon")]
        {
            gf127_row_add2_neon(c, a);
            return;
        }

        for i in 0..N {
            c[i] = Fq::add(c[i], a[i]);
        }
    }

    /// c[i] = a[i]-b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(1); N]);
    /// let b = Vector::<N>([Fq(0); N]);
    /// Vector::sub(&mut c, &a, &b);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn sub(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_sub_avx512(c, a, b);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_sub_avx2(c, a, b);
                }
            }
            return;
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_sub_neon(c, a, b);
            }
            return;
        }

        for i in 0..N {
            c[i] = Fq::sub(a[i], b[i]);
        }
    }

    /// c[i] -= a[i] mod q,  a[i],c[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(1); N]);
    /// let a = Vector::<N>([Fq(0); N]);
    /// Vector::sub2(&mut c, &a);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn sub2(c: &mut Vector<N>, a: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_sub2_avx512(c, a);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_sub2_avx2(c, a);
                }
            }
            return;
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_sub2_neon(c, a);
            }
            return;
        }

        for i in 0..N {
            c[i] = Fq::sub(c[i], a[i]);
        }
    }

    /// c[i] = a[i]*b[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(0); N]);
    /// let b = Vector::<N>([Fq(1); N]);
    /// Vector::mul(&mut c, &a, &b);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn mul(c: &mut Vector<N>, a: &Vector<N>, b: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_mul_avx512(c, a, b);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_mul_avx2(c, a, b);
                }
            }
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_mul_neon(c, a, b);
            }
            return;
        }
        for i in 0..N {
            c[i] = Fq::mul(a[i], b[i]);
        }
    }

    /// c[i] = c[i]*a[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Vector::<N>([Fq(0); N]);
    /// Vector::mul2(&mut c, &a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn mul2(c: &mut Vector<N>, a: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert!(N % 32 == 0);
                unsafe {
                    gf127_row_mul2_avx512(c, a);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert!(N % 32 == 0);
                unsafe {
                    gf127_row_mul2_avx2(c, a);
                }
            }
            return;
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_mul2_neon(c, a);
            }
            return;
        }

        for i in 0..N {
            c[i] = Fq::mul(c[i], a[i]);
        }

    }

    /// c[i] = c[i]*a mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut out = Vector::<N>([Fq(0); N]);
    /// let c = Vector::<N>([Fq(0); N]);
    /// let a = Fq(1);
    /// Vector::scalar(&mut out, &c, a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn scalar(out: &mut Vector<N>, c: &Vector<N>, a: Fq) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_scalar_mul2_avx512(out, c, a);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_scalar_mul2_avx2(out, c, a);
                }
            }

            return;
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_scalar_mul_neon(c, a);
            }
            return;
        }

        for j in 0..N {
            out[j] = Fq::mul(c[j], a);
        }
    }

    /// c[i] = c[i]*a mod q,  a,b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let a = Fq(1);
    /// Vector::scalar2(&mut c, a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    #[inline]
    pub fn scalar2(c: &mut Vector<N>, a: Fq) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_scalar_mul_avx512(c, a);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_scalar_mul_avx2(c, a);
                }
            }

            return;
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_scalar_mul_neon(c, a);
            }
            return;
        }

        for j in 0..N {
            c[j] = Fq::mul(c[j], a);
        }
    }

    /// computes: sum_{i=0}^{N} row[i] mod q,  row[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N>([Fq(0); N]);
    /// let t = Vector::acc(&mut c);
    /// assert_eq!(t.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `row`: row to accumulate
    #[inline]
    pub fn acc(row: &Vector<N>) -> Fq {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_avx512(row));
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_avx2(row));
                }
            }
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                return Fq(gf127_row_acc_neon(row));
            }
        }

        let mut t = row[0];
        for i in 1..N {
            t = t + row[i];
        }
        t
    }


    /// NOTE: non ct
    /// c[i] = a[i]^{-1} mod q, a[i] < 127 for i in 0..N
    /// # Examples
    /// ```
    /// use less::vector::Vector;
    /// const N: usize = 128;
    /// let mut row_out = Vector::<N>::default();
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::inv_non_ct(&mut row_out, &row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to invert
    pub fn inv_non_ct(row_out: &mut Vector<N>,
                      in1: &Vector<N>) {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    gf127_row_inv_avx512(row_out, in1);
                }
            }
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                gf127_row_inv_neon(row_out, in1);
            }
            return;
        }
        for col in 0..N {
            row_out[col] = Fq::inv_non_ct(in1[col]);
        }
    }

    /// computes: sum_{i=0}^{N} a[i]^{-1} mod q, a[i] < 127 for i in 0..N
    /// # Examples
    /// ```
    /// use less::vector::Vector;
    /// const N: usize = 128;
    /// let row = Vector::<N>::new();
    /// let t = Vector::acc_inv(&row);
    /// assert_eq!(t.0, 0)
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to invert
    #[inline]
    pub fn acc_inv(row: &Vector<N>) -> Fq {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_inv_avx512(row));
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_inv_avx2(row));
                }
            }
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                return Fq(gf127_row_acc_inv_neon(row));
            }
        }

        let mut t = row[0];
        for i in 1..N {
            t = t + Fq::inv(row[i]);
        }
        t
    }
    
    /// c[i] == c[j] for all i < j < n
    /// # Examples
    /// ```
    /// use less::vector::Vector;
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::all_same(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to check if all elements are equal
    ///
    /// # Return 
    /// -   1 if all elements are the same
    /// -   0 else
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
    /// use less::vector::Vector;
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::contain_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row 
    ///
    /// # Return 
    /// -   1 if a zero is in the row
    /// -   0 else
    pub fn contain_zero(row: &Vector<N>) -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return gf127_row_contains_zero_avx512(row);
                }
            } else if is_x86_feature_detected!("avx2") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return gf127_row_contains_zero_avx2(row);
                }
            }
        }
        #[cfg(target_feature = "neon")]
        {
            unsafe {
                return gf127_row_contains_zero_neon(row);
            }
        }

        for col in 0..N {
            if row[col] == Fq(0) {
                return true;
            }
        }
        false
    }

    /// c[i] == 0 for all i  < n
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// const N: usize = 128;
    /// let row1 = Vector::<N>::from_u8(1);
    /// Vector::<N>::count_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row 
    ///
    /// # return 
    /// -   number of zeros in row
    pub fn count_zero(row: &Vector<N>) -> u32 {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512f") {
                assert_eq!(N % 32, 0);
                unsafe {
                    return gf127_row_count_zero_avx2(row);
                }
            }
        }

        #[cfg(target_feature = "neon")]
        {
            unsafe {
                return gf127_row_count_zero_neon(row);
            }
        }

        let mut c: u32 = 0;
        for col in 0..N {
            if row[col] == Fq(0) {
                c += 1;
            }
        }
        c
    }

    /// NOTE: assumes the vectors are NOT in multiset/histogram from
    /// \input row1[in]: pointer to the first row
    /// \input row2[in]: pointer to the second row
    /// \return: 0 if (row1) == (row2)
    ///          x if (row1) > (row2)
    ///         -x if (row1) < (row2)
    pub fn cmp(a: &Vector<N>, b: &Vector<N>) -> i32 {
        let mut i = 0;
        while i < (N-1) && a[i] == b[i] {
            i += 1;
        }
        (b[i].0 as i32) - (a[i].0 as i32)
    }
}


impl <const N: usize> Add for Vector<N> {
    type Output = Vector<N>;
    #[inline]
    fn add(self, r: Vector<N>) -> Vector<N> {
        let mut ret = Vector::default();
        Vector::<N>::add(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize> AddAssign for Vector<N> {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        Vector::<N>::add2(self, &other);
    }
}

impl <const N: usize> Sub for Vector<N> {
    type Output = Vector<N>;
    #[inline]
    fn sub(self, r: Vector<N>) -> Vector<N> {
        let mut ret = Vector::default();
        Vector::<N>::sub(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize> SubAssign for Vector<N> {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        Vector::<N>::sub2(self, &other);
    }
}

impl <const N: usize> Mul for Vector<N> {
    type Output = Vector<N>;
    #[inline]
    fn mul(self, r: Vector<N>) -> Vector<N> {
        let mut ret = Vector::default();
        Vector::<N>::mul(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize> MulAssign for Vector<N> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        Vector::<N>::mul2(self, &other);
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
impl<const N: usize> DerefMut for Vector<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const N: usize> Display for Vector<N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for i in 0..N {
            write!(f, "{}", self[i])?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 32;
    
    #[test]
    fn init() {
        let c = Vector::<N>::default();
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 0);
        }
    }

    #[test]
    fn from_u8() {
        let c = Vector::<N>::from_u8(1);
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 1);
        }
    }

    #[test]
    fn from_vector() {
        let t = [Fq(0u8); N];
        let c = Vector::<N>::from_vector(&t);
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 0);
        }
    }

    #[test]
    fn dimension() {
        let c = Vector::<N>::default();
        assert_eq!(c.dimension(), N);
    }

    #[test]
    fn add() {
        for i in 0..127 {
            for j in 0..127 {
                let t = (i+j)%127;
                let mut c = Vector::<N>([Fq(0); N]);
                let a = Vector::<N>([Fq(i); N]);
                let b = Vector::<N>([Fq(j); N]);
                Vector::add(&mut c, &a, &b);
                for i in 0..N {
                    assert_eq!(c[i].0, t);
                }
            }
        }
    }

    #[test]
    fn sub() {
        for i in 0..127 {
            for j in 0..127 {
                let t = ((i+127) - j)%127;
                let mut c = Vector::<N>([Fq(0); N]);
                let a = Vector::<N>([Fq(i); N]);
                let b = Vector::<N>([Fq(j); N]);
                Vector::sub(&mut c, &a, &b);
                for i in 0..N {
                    assert_eq!(c[i].0, t);
                }
            }
        }
    }

    #[test]
    fn mul() {
        for i in 0..127u16 {
            for j in 0..127u16{
                let t = ((i*j) % 127) as u8;
                let mut c = Vector::<N>([Fq(0); N]);
                let a = Vector::<N>([Fq(i as u8); N]);
                let b = Vector::<N>([Fq(j as u8); N]);
                Vector::mul(&mut c, &a, &b);
                for i in 0..N {
                    assert_eq!(c[i].0, t);
                }
            }
        }
    }

    #[test]
    fn scalar() {
        for i in 0..127u16 {
            for j in 0..127u16{
                let t = ((i*j) % 127) as u8;
                let mut c = Vector::<N>([Fq(0); N]);
                let a = Vector::<N>([Fq(i as u8); N]);
                Vector::scalar(&mut c, &a, Fq(j as u8));
                for i in 0..N {
                    assert_eq!(c[i].0, t);
                }
            }
        }
    }

    #[test]
    fn acc() {
        let c = Vector::<N>::default();
        assert_eq!(Vector::acc(&c).0, 0);

        let c = Vector::<N>::from_u8(1);
        assert_eq!(Vector::acc(&c), Fq::red(N as u16));
    }

    #[test]
    fn inv() {
        let a = Vector::<N>::default();
        let mut c = Vector::<N>::default();
        Vector::inv_non_ct(&mut c, &a);
        for i in 0..N {
            assert_eq!(c[i], Fq(0));
        }
    }

    #[test]
    fn all_same() {
        let c = Vector::<N>::default();
        assert_eq!(Vector::all_same(&c), true);
        let c = Vector::<N>::from_u8(1);
        assert_eq!(Vector::all_same(&c), true);
        let mut c = Vector::<N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::all_same(&c), false);
    }

    #[test]
    fn contain_zero() {
        let c = Vector::<N>::default();
        assert_eq!(Vector::contain_zero(&c), true);
        let c = Vector::<128>::from_u8(1);
        assert_eq!(Vector::contain_zero(&c), false);
        let mut c = Vector::<N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::contain_zero(&c), true);
    }

    #[test]
    fn count_zero() {
        let c = Vector::<N>::default();
        assert_eq!(Vector::count_zero(&c), N as u32);
        let c = Vector::<N>::from_u8(1);
        assert_eq!(Vector::count_zero(&c), 0);
        let mut c = Vector::<N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::count_zero(&c), 1);
    }
}
