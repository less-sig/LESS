#![allow(dead_code)]
#![allow(unreachable_code)]

use sha3::digest::XofReader;
use std::ops:: {
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign,
    Index, IndexMut, Deref, DerefMut
};
use std::cmp::Ordering;
use std::{ fmt::Display, fmt::Formatter, fmt::Result };
use crate::config::{
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
#[repr(align(64))] // NOTE alignment only for avx2/avx512
pub struct Vector<const N: usize, const N_PAD: usize>(pub [Fq; N_PAD]);

impl <const N: usize, const N_PAD: usize> Vector<N, N_PAD>{
    pub const Q: u8 = 127;
    pub const Q_BITS:u8 = 7;

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// let result: Vector<100, 100> = Vector::default();
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
        rand_range_q_state_elements(&mut a.0[..N], state);
        a
    }

    /// sample a random vector
    pub fn rand_star<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut a = Self::default();
        fq_star_rnd_state_elements(&mut a.0[..N], state);
        a
    }

    /// Create an instance from a `Vec`
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// let t = core::array::from_fn(|_| Fq::default());
    /// let result: Vector<100, 100> = Vector::from_vector(&t);
    /// assert_eq!(result[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `coefficients`: The input `Vec`.
    ///
    /// # Returns
    /// A new vector element
    #[inline]
    pub fn from_vector(coefficients: &[Fq; N_PAD]) -> Self {
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

    /// partially compares two vectors
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let A: Vector<100, 100> = Vector::default();
    /// let B: Vector<100, 100> = Vector::default();
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
    fn partial_cmp(a: &Self, b: &Self) -> Option<Ordering> {
        let mut aa = Multiset::<Q_PAD>::default();
        let mut bb = Multiset::<Q_PAD>::default();
        Multiset::from_row(&mut aa, a);
        Multiset::from_row(&mut bb, b);
        Some(Multiset::<Q_PAD>::partial_cmp(&aa, &bb).cmp(&0))
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let A: Vector<100, 100> = Vector::default();
    /// let B: Vector<100, 100> = Vector::default();
    /// A == B;
    /// A != B;
    /// ```
    ///
    /// # Returns
    ///  1 if a == b, else 0
    fn eq(a: &Self, b: &Self) -> bool {
        Self::partial_cmp(a, b).unwrap().is_eq()
    }

    /// returns the dimension of the vector
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// let result: Vector<100, 100> = Vector::default();
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Vector::<N, N>([Fq(0); N]);
    /// let b = Vector::<N, N>([Fq(1); N]);
    /// Vector::add(&mut c, &a, &b);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn add(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>, b: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_add_avx512(c, a, b);
                }
            } else if cfg!(target_feature = "avx2") {
                assert!(N_PAD % 32 == 0);
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

        for i in 0..N_PAD {
            c[i] = Fq::add(a[i], b[i]);
        }
    }

    /// c[i] += a[i] mod q,  a[i],b[i] < 127, for i in 0..N
    /// # Examples
    ///
    /// ```
    /// use less::vector::Vector;
    /// use less::fq::Fq;
    /// const N: usize = 64;
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Vector::<N, N>([Fq(1); N]);
    /// Vector::add2(&mut c, &a);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn add2(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_add2_avx512(c, a);
                }
            } else if cfg!(target_feature = "avx2") {
                assert!(N_PAD % 32 == 0);
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

        for i in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Vector::<N, N>([Fq(1); N]);
    /// let b = Vector::<N, N>([Fq(0); N]);
    /// Vector::sub(&mut c, &a, &b);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: second addend
    #[inline]
    pub fn sub(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>, b: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_sub_avx512(c, a, b);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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

        for i in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(1); N]);
    /// let a = Vector::<N, N>([Fq(0); N]);
    /// Vector::sub2(&mut c, &a);
    /// assert_eq!(c[0].0, 1);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: minuend
    #[inline]
    pub fn sub2(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_sub2_avx512(c, a);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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

        for i in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Vector::<N, N>([Fq(0); N]);
    /// let b = Vector::<N, N>([Fq(1); N]);
    /// Vector::mul(&mut c, &a, &b);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first multiplier
    /// - `b`: second multiplier
    #[inline]
    pub fn mul(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>, b: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_mul_avx512(c, a, b);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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
        for i in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Vector::<N, N>([Fq(0); N]);
    /// Vector::mul2(&mut c, &a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first multiplier
    /// - `b`: second multiplier
    #[inline]
    pub fn mul2(c: &mut Vector<N, N_PAD>, a: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert!(N_PAD % 32 == 0);
                unsafe {
                    gf127_row_mul2_avx512(c, a);
                }
            } else if cfg!(target_feature = "avx2") {
                assert!(N_PAD % 32 == 0);
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

        for i in 0..N_PAD {
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
    /// let mut out = Vector::<N, N>([Fq(0); N]);
    /// let c = Vector::<N, N>([Fq(0); N]);
    /// let a = Fq(1);
    /// Vector::scalar(&mut out, &c, a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: first addend
    /// - `b`: scalar
    #[inline]
    pub fn scalar(out: &mut Vector<N, N_PAD>, c: &Vector<N, N_PAD>, a: Fq) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_scalar_mul2_avx512(out, c, a);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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

        for j in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let a = Fq(1);
    /// Vector::scalar2(&mut c, a);
    /// assert_eq!(c[0].0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `c`: output value
    /// - `a`: scalar value
    #[inline]
    pub fn scalar2(c: &mut Vector<N, N_PAD>, a: Fq) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    gf127_row_scalar_mul_avx512(c, a);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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

        for j in 0..N_PAD {
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
    /// let mut c = Vector::<N, N>([Fq(0); N]);
    /// let t = Vector::acc(&mut c);
    /// assert_eq!(t.0, 0);
    /// ```
    ///
    /// # Parameters
    /// - `row`: row to accumulate
    #[inline]
    pub fn acc(row: &Vector<N, N_PAD>) -> Fq {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_avx512(row));
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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
        for i in 1..N_PAD {
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
    /// let mut row_out = Vector::<N, N>::default();
    /// let row1 = Vector::<N, N>::from_u8(1);
    /// Vector::inv_non_ct(&mut row_out, &row1);
    /// ```
    ///
    /// # Parameters
    /// - `row_out`: output row
    /// - `in1`: row to invert
    pub fn inv_non_ct(row_out: &mut Vector<N, N_PAD>,
                      in1: &Vector<N, N_PAD>) {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
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
    /// let row = Vector::<N, N>::default();
    /// let t = Vector::acc_inv(&row);
    /// assert_eq!(t.0, 0)
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to invert
    #[inline]
    pub fn acc_inv(row: &Vector<N, N_PAD>) -> Fq {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    return Fq(gf127_row_acc_inv_avx512(row));
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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
    /// let row1 = Vector::<N, N>::from_u8(1);
    /// Vector::<N, N>::all_same(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row to check if all elements are equal
    ///
    /// # Return
    /// -   1 if all elements are the same
    /// -   0 else
    pub fn all_same(row: &Vector<N, N_PAD>) -> bool {
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
    /// let row1 = Vector::<N, N>::from_u8(1);
    /// Vector::<N, N>::contain_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row
    ///
    /// # Return
    /// -   1 if a zero is in the row
    /// -   0 else
    pub fn contain_zero(row: &Vector<N, N_PAD>) -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
                unsafe {
                    return gf127_row_contains_zero_avx512(row);
                }
            } else if cfg!(target_feature = "avx2") {
                assert_eq!(N_PAD % 32, 0);
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
    /// let row1 = Vector::<N, N>::from_u8(1);
    /// Vector::<N, N>::count_zero(&row1);
    /// ```
    ///
    /// # Parameters
    /// - `a`: row
    ///
    /// # return
    /// -   number of zeros in row
    pub fn count_zero(row: &Vector<N, N_PAD>) -> u32 {
        #[cfg(target_arch = "x86_64")]
        {
            if cfg!(all(target_feature = "avx512f", target_feature = "avx512bw")) {
                assert_eq!(N_PAD % 32, 0);
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
    pub fn cmp(a: &Vector<N, N_PAD>, b: &Vector<N, N_PAD>) -> i32 {
        let mut i = 0;
        while i < (N-1) && a[i] == b[i] {
            i += 1;
        }
        (b[i].0 as i32) - (a[i].0 as i32)
    }
}


impl <const N: usize, const N_PAD: usize> Add for Vector<N, N_PAD> {
    type Output = Vector<N, N_PAD>;
    #[inline]
    fn add(self, r: Vector<N, N_PAD>) -> Vector<N, N_PAD> {
        let mut ret = Vector::default();
        Vector::<N, N_PAD>::add(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize, const N_PAD: usize> AddAssign for Vector<N, N_PAD> {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        Vector::<N, N_PAD>::add2(self, &other);
    }
}

impl <const N: usize, const N_PAD: usize> Sub for Vector<N, N_PAD> {
    type Output = Vector<N, N_PAD>;
    #[inline]
    fn sub(self, r: Vector<N, N_PAD>) -> Vector<N, N_PAD> {
        let mut ret = Vector::default();
        Vector::<N, N_PAD>::sub(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize, const N_PAD: usize> SubAssign for Vector<N, N_PAD> {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        Vector::<N, N_PAD>::sub2(self, &other);
    }
}

impl <const N: usize, const N_PAD: usize> Mul for Vector<N, N_PAD> {
    type Output = Vector<N, N_PAD>;
    #[inline]
    fn mul(self, r: Vector<N, N_PAD>) -> Vector<N, N_PAD> {
        let mut ret = Vector::default();
        Vector::<N, N_PAD>::mul(&mut ret, &self, &r);
        ret
    }
}
impl <const N: usize, const N_PAD: usize> MulAssign for Vector<N, N_PAD> {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        Vector::<N, N_PAD>::mul2(self, &other);
    }
}


impl<const N: usize, const N_PAD: usize> Index<usize> for Vector<N, N_PAD> {
    type Output = Fq;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
impl<const N: usize, const N_PAD: usize> IndexMut<usize> for Vector<N, N_PAD> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const N: usize, const N_PAD: usize> Deref for Vector<N, N_PAD> {
    type Target = [Fq];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<const N: usize, const N_PAD: usize> DerefMut for Vector<N, N_PAD> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const N: usize, const N_PAD: usize> Display for Vector<N, N_PAD> {
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
    use crate::config::Xof;
    use crate::prng::csprng_initialize;
    const N: usize = 32;
    
    #[test]
    fn init() {
        let c = Vector::<N, N>::default();
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 0);
        }
    }

    #[test]
    fn from_u8() {
        let c = Vector::<N, N>::from_u8(1);
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 1);
        }
    }

    #[test]
    fn from_vector() {
        let t = [Fq(0u8); N];
        let c = Vector::<N, N>::from_vector(&t);
        for i in 0..c.dimension() {
            assert_eq!(c[i].0, 0);
        }
    }

    #[test]
    fn dimension() {
        let c = Vector::<N, N>::default();
        assert_eq!(c.dimension(), N);
    }

    #[test]
    fn add() {
        for i in 0..127 {
            for j in 0..127 {
                let t = (i+j)%127;
                let mut c = Vector::<N, N>([Fq(0); N]);
                let a = Vector::<N, N>([Fq(i); N]);
                let b = Vector::<N, N>([Fq(j); N]);
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
                let mut c = Vector::<N, N>([Fq(0); N]);
                let a = Vector::<N, N>([Fq(i); N]);
                let b = Vector::<N, N>([Fq(j); N]);
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
                let mut c = Vector::<N, N>([Fq(0); N]);
                let a = Vector::<N, N>([Fq(i as u8); N]);
                let b = Vector::<N, N>([Fq(j as u8); N]);
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
                let mut c = Vector::<N, N>([Fq(0); N]);
                let a = Vector::<N, N>([Fq(i as u8); N]);
                Vector::scalar(&mut c, &a, Fq(j as u8));
                for i in 0..N {
                    assert_eq!(c[i].0, t);
                }
            }
        }
    }

    #[test]
    fn acc() {
        let c = Vector::<N, N>::default();
        assert_eq!(Vector::acc(&c).0, 0);

        let c = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::acc(&c), Fq::red(N as u16));
    }

    #[test]
    fn inv() {
        let a = Vector::<N, N>::default();
        let mut c = Vector::<N, N>::default();
        Vector::inv_non_ct(&mut c, &a);
        for i in 0..N {
            assert_eq!(c[i], Fq(0));
        }
    }

    #[test]
    fn all_same() {
        let c = Vector::<N, N>::default();
        assert_eq!(Vector::all_same(&c), true);
        let c = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::all_same(&c), true);
        let mut c = Vector::<N, N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::all_same(&c), false);
    }

    #[test]
    fn contain_zero() {
        let c = Vector::<N, N>::default();
        assert_eq!(Vector::contain_zero(&c), true);
        let c = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::contain_zero(&c), false);
        let mut c = Vector::<N, N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::contain_zero(&c), true);
    }

    #[test]
    fn count_zero() {
        let c = Vector::<N, N>::default();
        assert_eq!(Vector::count_zero(&c), N as u32);
        let c = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::count_zero(&c), 0);
        let mut c = Vector::<N, N>::from_u8(1);
        c[N-1] = Fq(0);
        assert_eq!(Vector::count_zero(&c), 1);

        // a length which is not a multiple of 32, the scalar fallback
        let mut c = Vector::<126, 126>::from_u8(5);
        c[0] = Fq(0);
        c[125] = Fq(0);
        assert_eq!(Vector::count_zero(&c), 2);
    }

    #[test]
    fn rand() {
        let mut s = csprng_initialize::<Xof>(&[1u8; 32]);
        let a = Vector::<N, N>::rand(&mut s);
        for i in 0..N {
            assert!(a[i].0 < 127);
        }

        let mut s = csprng_initialize::<Xof>(&[1u8; 32]);
        let b = Vector::<N, N>::rand(&mut s);
        assert_eq!(a, b);
    }

    #[test]
    fn rand_star() {
        for seed in 0..32u8 {
            let mut s = csprng_initialize::<Xof>(&[seed; 32]);
            let a = Vector::<N, N>::rand_star(&mut s);
            for i in 0..N {
                assert!(a[i].0 >= 1);
                assert!(a[i].0 < 127);
            }
        }
    }

    #[test]
    fn add2() {
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[3u8; 32]));
        let b = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[4u8; 32]));
        let mut c = Vector::<N, N>::from_vector(&a.0);
        Vector::add2(&mut c, &b);
        for i in 0..N {
            assert_eq!(c[i], Fq::add(a[i], b[i]));
        }
    }

    #[test]
    fn sub2() {
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[5u8; 32]));
        let b = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[6u8; 32]));
        let mut c = Vector::<N, N>::from_vector(&a.0);
        Vector::sub2(&mut c, &b);
        for i in 0..N {
            assert_eq!(c[i], Fq::sub(a[i], b[i]));
        }
    }

    #[test]
    fn mul2() {
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[7u8; 32]));
        let b = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[8u8; 32]));
        let mut c = Vector::<N, N>::from_vector(&a.0);
        Vector::mul2(&mut c, &b);
        for i in 0..N {
            assert_eq!(c[i], Fq::mul(a[i], b[i]));
        }
    }

    #[test]
    fn scalar2() {
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[9u8; 32]));
        for s in 0..127u8 {
            let mut c = Vector::<N, N>::from_vector(&a.0);
            Vector::scalar2(&mut c, Fq(s));
            for i in 0..N {
                assert_eq!(c[i], Fq::mul(a[i], Fq(s)));
            }
        }
    }

    #[test]
    fn acc_inv() {
        let a = Vector::<N, N>::default();
        assert_eq!(Vector::acc_inv(&a).0, 0);

        let a = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::acc_inv(&a), Fq::red(N as u16));

        // every element is inverted, the accumulator starts at zero
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[11u8; 32]));
        let mut s = Fq(0);
        for i in 0..N {
            s = Fq::add(s, Fq::inv_non_ct(a[i]));
        }
        assert_eq!(Vector::acc_inv(&a), s);
    }

    #[test]
    fn cmp() {
        let a = Vector::<N, N>::from_u8(1);
        let b = Vector::<N, N>::from_u8(1);
        assert_eq!(Vector::cmp(&a, &b), 0);

        let mut c = Vector::<N, N>::from_u8(1);
        c[0] = Fq(2);
        assert!(Vector::cmp(&a, &c) > 0);
        assert!(Vector::cmp(&c, &a) < 0);
    }

    #[test]
    fn partial_cmp_and_eq() {
        let a = Vector::<N, N>::from_u8(1);
        let mut b = Vector::<N, N>::from_u8(1);
        assert!(Vector::eq(&a, &b));

        b[0] = Fq(2);
        b[1] = Fq(2);
        assert!(!Vector::eq(&a, &b));
        assert_eq!(Vector::partial_cmp(&b, &a), Some(Ordering::Greater));
        assert_eq!(Vector::partial_cmp(&a, &b), Some(Ordering::Less));

        // a permutation does not change the multiset
        let mut c = Vector::<N, N>::default();
        for i in 0..N { c[i] = Fq(i as u8 % 3); }
        let mut d = Vector::<N, N>::default();
        for i in 0..N { d[N - 1 - i] = c[i]; }
        assert!(Vector::eq(&c, &d));
    }

    #[test]
    fn operators() {
        let a = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[0u8; 32]));
        let b = Vector::<N, N>::rand(&mut csprng_initialize::<Xof>(&[1u8; 32]));

        let c = Vector::<N, N>::from_vector(&a.0) + Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::add(a[i], b[i])); }

        let c = Vector::<N, N>::from_vector(&a.0) - Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::sub(a[i], b[i])); }

        let c = Vector::<N, N>::from_vector(&a.0) * Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::mul(a[i], b[i])); }

        let mut c = Vector::<N, N>::from_vector(&a.0);
        c += Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::add(a[i], b[i])); }

        let mut c = Vector::<N, N>::from_vector(&a.0);
        c -= Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::sub(a[i], b[i])); }

        let mut c = Vector::<N, N>::from_vector(&a.0);
        c *= Vector::<N, N>::from_vector(&b.0);
        for i in 0..N { assert_eq!(c[i], Fq::mul(a[i], b[i])); }
    }

    #[test]
    fn index_and_deref() {
        let mut a = Vector::<N, N>::default();
        a[3] = Fq(17);
        assert_eq!(a[3].0, 17);
        assert_eq!(a.len(), N);
        assert_eq!(a.iter().filter(|x| x.0 == 17).count(), 1);
        a.iter_mut().for_each(|x| *x = Fq(1));
        assert_eq!(a[3].0, 1);
    }

    #[test]
    fn display() {
        let a = Vector::<3, 3>::from_u8(7);
        assert_eq!(format!("{}", a), "777");
    }

    /// the derived Ord is the plain lexicographic order, the column sort
    /// relies on it
    #[test]
    fn derived_ord() {
        let a = Vector::<N, N>::from_u8(1);
        let mut b = Vector::<N, N>::from_u8(1);
        b[5] = Fq(2);
        assert!(a < b);
        assert!(b > a);
    }
}
