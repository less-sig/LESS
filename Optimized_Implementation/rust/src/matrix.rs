#![allow(dead_code)]

use std::{
    cmp,
    fmt,
    ops::{Add, Index, IndexMut},
};

use crate::fq::Fq;
use crate::monomial::Monomial;
use crate::vector::Vector;


const Q_PAD: usize = 128;
const Q: u8 = 127;

pub struct Matrix<const N: usize, const M: usize> {
    rows: [Vector<M>; N],
}

impl<const N: usize, const M: usize> Matrix<N, M> {

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::init();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A new Matrix
    #[inline]
    pub fn init() -> Self {
        Self {
            rows: core::array::from_fn(|_| Vector::<M>::init() ),
        }
    }

    /// same as init
    #[inline]
    pub fn new() -> Matrix<N, M> {
        Matrix {
            rows: core::array::from_fn(|_| Vector::<M>::init() ),
        }
    }

    /// Identity Matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::id();
    /// assert_eq!(result[0][0].0, 1);
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A new identity matrix
    pub fn id() -> Self {
        let mut a = Self::init();
        for i in 0..cmp::min(N, M) {
            a[i][i] = Fq(1);
        }
        a
    }

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::init();
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A size/dimension of the vector
    #[inline]
    pub fn dimension(&self) -> (usize, usize) {
        (N, M)
    }

    /// Return the number of columns
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::init();
    /// assert_eq!(result.ncols(), 100);
    /// ```
    ///
    /// # Returns
    /// number of columns
    #[inline]
    pub fn ncols(&self) -> usize {
        M
    }

    /// Return the number of rows
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::init();
    /// assert_eq!(result.nrows(), 100);
    /// ```
    ///
    /// # Returns
    /// n_rows
    #[inline]
    pub fn nrows(&self) -> usize {
        N
    }
    
    /// constant time addition
    ///     r += l
    ///
    /// # Examples
    ///
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `l`: NxM matrix
    #[inline]
    pub fn add(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<M>::add2(&mut r[i],&l[i]);
        }
    }

    /// constant time subtraction
    /// r -= l
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let mut C: Matrix<100, 100> = Matrix::init();
    /// let A: Matrix<100, 100> = Matrix::init();
    /// Matrix<100, 100>::sub(&mut C, &A)
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `l`: NxM matrix
    #[inline]
    pub fn sub(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<M>::sub2(&mut r[i],&l[i]);
        }
    }

    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    #[inline]
    pub fn swap_rows(r: &mut Self, r1: usize, r2: usize) {
        for j in 0..M {
            let tmp = r.rows[r1][j];
            r.rows[r1][j] = r.rows[r2][j];
            r.rows[r2][j] = tmp;
        } 
    }

    /// TODO not ct
    pub fn monomial_mul(res: &mut Self, g: &Self, monom: &Monomial<N>) {
        for src_col_idx in 0..M {
            for row_idx in 0..N {
                let pos = monom.perms[src_col_idx];
                let a = g.rows[row_idx][src_col_idx];
                let b = monom.coeffs[src_col_idx];
                res.rows[row_idx][pos as usize] = Fq::mul(a, b);
            }
        }
    }

    pub fn rref_ct(g: &mut Self,
                   is_pivot_column: &mut [u8; M],
                   was_pivot_column: &[u8; M],
                   pvt_reuse_limit: u32) {
        let mut pvt_reuse_limit: u32 = 0;

        if pvt_reuse_limit != 0 {
            // TODO: loop boundaries propably not correct
            for preproc_col in N - 1 .. 0 {
                if was_pivot_column[preproc_col as usize] == 1 {
                    let mut pivot_el_row = 0;
                    for row in 0..N {
                        if g.rows[row][preproc_col] != Fq(0) {
                            pivot_el_row = row;
                        }
                    }

                    Self::swap_rows(g, preproc_col, pivot_el_row);
                }
            }
        }

        for row_to_reduce in 0..N {
            // TODO
        }
    }

}


impl<const N: usize, const M: usize> Add for Matrix<N, M> {
    type Output = Matrix<N, M>;
    fn add(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::add(&mut t, &r);
        t
    }
}

impl<'a, const N: usize, const M: usize> Add for &'a mut Matrix<N, M> {
    type Output = &'a mut Matrix<N, M>;
    fn add(self, r: Self) -> Self::Output {
        Matrix::<N, M>::add(self, r);
        self
    }
}


/// Direct access to a row
impl<const N: usize, const M: usize> Index<usize> for Matrix<N, M> {
    type Output = Vector<M>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}
/// Direct access to a row (mutable)
impl<const N: usize, const M: usize> IndexMut<usize> for Matrix<N, M> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.rows[index]
    }
}
impl<const N: usize, const M: usize> fmt::Debug for Matrix<N, M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:?}", self.rows)
    }
}


#[derive(Debug)]
pub struct MatrixNonIS<const N: usize> {
    rows: [[Fq; N]; N],
}

impl<const N: usize> MatrixNonIS<N> {
    /// generates a new matriz init with 0
    pub fn new() -> MatrixNonIS<N> {
        MatrixNonIS {
            rows: [[Fq(0); N]; N],
        }
    }
   
    /// generates a new matrix init with q-1
    pub fn new_large() -> MatrixNonIS<N> {
        MatrixNonIS {
            rows: [[Fq(Q-1); N]; N],
        }
    }

    /// static internal function implementing the histogram function
    fn sort(out: &mut [u8; Q_PAD], row: &[Fq; N]) {
        for i in 0..N {
            out[row[i].0 as usize] += 1;
        }
    }

    /// static internal function comparing two rows (histogram form)
    fn compare_rows(a: &[u8; Q_PAD], b: &[u8; Q_PAD]) -> i32 {
        let mut i: usize = 0;
        while ((i as u8) < (Q-1)) && (a[i] == b[i]) {
            i += 1;
        };

        return (b[i] as i32) - (a[i]  as i32);
    }
   
    // not working as i cannot borrow each row seperatily as mut.
    // fn row_mul(row: &mut [Fq; N], s: Fq) {
    //     for i in 0..N {
    //         row[i] = Fq::mul(row[i], s);
    //     }
    // }

    /// rows[i] *= s 
    fn row_mul(&mut self, i: usize, s: Fq) {
        for j in 0..N {
            self.rows[i][j] = Fq::mul(self.rows[i][j], s);
        }
    }

    /// \return sum(row)
    fn row_acc(row: &[Fq; N]) -> Fq {
        let mut t = row[0];
        for i in 1..N {
            t = t + row[i];
        }
        t 
    }

    /// \return sum(row**{-1})
    fn row_acc_inv(row: &[Fq; N]) -> Fq {
        let mut t = Fq::inv(row[0]);
        for i in 1..N {
            t = t + Fq::inv(row[i]);
        }
        t 
    }

    /// insert sort 
    fn sort_rows(&mut self) -> i32 {

        let mut tmp = [[0u8; Q_PAD]; N];

        for i in 0..N {
            Self::sort(&mut tmp[i], &mut self.rows[i]);
        }

        for i in 1..N {
            let mut j = i;
            while (j > 0) && (Self::compare_rows(&tmp[j-1], &tmp[j]) < 0) {
                if let Ok([a, b]) = &mut tmp.get_disjoint_mut([j-1, j]) {
                    std::mem::swap(a, b);
                }
                j -= 1;
             }
        }
        0
    }

    fn sort_cols(&mut self) -> i32 {
        0
    }

    fn sort_cf(&mut self) -> i32 {
        if self.sort_rows() != 0 {
            return 0;
        }
        self.sort_cols();
        1
    }

    ///
    fn scale_cf(&mut self) -> i32 {
        for i in 0..N  {
            let mut s = Self::row_acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::row_acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_mul(i as usize, s);
        };

        return self.sort_cf();
    }

    /// by floyd
    fn scale_cf_preprocess_v1(&mut self, z: usize, m: & [u8; Q_PAD]) -> i32 {
        let mut tmp = [0u8; Q_PAD];
        for i in 0..z  {
            let mut s = Self::row_acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::row_acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_mul(i as usize, s);
            Self::sort(&mut tmp, &self.rows[i]);
            if Self::compare_rows(&tmp, m) < 0 {
                return 0;
            }
        };

        return 0;
    }

    /// by luke
    fn scale_cf_preprocess_v2(&mut self, z: usize, m: &mut [u8; Q_PAD]) -> i32 {
        let mut tmp = [0u8; Q_PAD];
        let mut ret = 0;
        for i in 0..z  {
            let mut s = Self::row_acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::row_acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_mul(i as usize, s);
            Self::sort(&mut tmp, &self.rows[i]);
            if Self::compare_rows(&tmp, m) < 0 {
                ret = 1;
                for i in 0..Q_PAD {
                    m[i] = tmp[i];
                }
            }
        };

        ret
    }

    pub fn cf_opt1(&mut self) {
        let mut M = Self::new_large();
        let mut J = [0; N];
        let mut z = 0;
        let mut num_zeros = 0;
    }
}
