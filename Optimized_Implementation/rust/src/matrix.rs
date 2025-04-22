#![allow(dead_code)]

use std::ops::Add;
use crate::fq::Fq;


const Q_PAD: usize = 128;
const Q: u8 = 127;

#[derive(Debug)]
pub struct Matrix<const N: usize, const M: usize> {
    rows: [[Fq; M]; N],
}

impl<const N: usize, const M: usize> Matrix<N, M> {
    pub fn new() -> Matrix<N, M> {
        Matrix {
            rows: [[Fq(0); M]; N],
        }
    }
    
    // unneded
    // fn sub(r: &mut Self, l: &Self) {
    //     for i in 0..N {
    //         for j in 0..M {
    //             r.rows[i][j] = Fq::sub(r.rows[i][j], l.rows[i][j]);
    //         } 
    //     }
    // }

    fn add(r: &mut Self, l: &Self) {
        for i in 0..N {
            for j in 0..M {
                r.rows[i][j] = Fq::add(r.rows[i][j], l.rows[i][j]);
            } 
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
