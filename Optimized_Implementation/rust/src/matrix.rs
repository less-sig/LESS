// TODO remove
#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use std::{
    cmp,
    fmt,
    ops::{
        Add, Sub, Mul, Index, IndexMut
    },
};
use std::{ fmt::Display, fmt::Formatter, fmt::Result };
use std::ptr::copy_nonoverlapping;
use sha3::{
    Digest, Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::config::{
    K, N8,
};
use crate::constants::{
    Q, Q_PAD,
};
use crate::fq::Fq;
use crate::monomial::{Monomial, Permutation};
use crate::vector::Vector;
use crate::multiset::Multiset;
use crate::prng::rand_range_q_state_elements;

/// row major form
/// N: number of rows
/// M: number of cols
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Matrix<const N: usize, const M: usize, const M_PAD: usize> {
    rows: [Vector<M, M_PAD>; N],
}


impl<const N: usize, const M: usize, const M_PAD: usize> Matrix<N, M, M_PAD> {

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::default();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A new Matrix
    #[inline]
    pub fn default() -> Self {
        Self {
            rows: core::array::from_fn(|_| Vector::<M, M_PAD>::default() ),
        }
    }

    /// same as init, but inits all values with `v`
    #[inline]
    pub fn from_u8(v: u8) -> Matrix<N, M, M_PAD> {
        Matrix {
            rows: core::array::from_fn(|_| Vector::<M, M_PAD>::from_u8(v) ),
        }
    }

    /// translation of `generator_rref_expand`
    /// Expands a compressed RREF generator matrix into a full one
    /// \param full[out]: output generator matrix (K \times N)
    /// \param compact[out]: input compressed generator matrix (K \times N-K)
    pub fn from_compressed_rref<const M_PRIME: usize, const M_PRIME_PAD: usize>(column_pos: &[u16; M_PRIME_PAD],compact: &Matrix<N, M_PRIME, M_PRIME_PAD>) -> Self {
        // static_assert!(M <= M_PRIME);
        let mut placed_dense_cols = 0usize;
        let mut full = Self::default();
        // only the first `crate::config::N` columns belong to the code, the
        // remaining ones of the padded row stay zero
        for col_idx in 0..crate::config::N {
            if (placed_dense_cols < crate::config::N - crate::config::K) && (col_idx as u16 == column_pos[placed_dense_cols]) {
                // non-pivot column, restore one full column
                for row_idx in 0..N {
                    full.rows[row_idx][col_idx] = compact.rows[row_idx][placed_dense_cols];
                }
                placed_dense_cols += 1;
            } else {
                // regenerate the appropriate pivot column
                for row_idx in 0..N {
                    full.rows[row_idx][col_idx] = Fq((row_idx == col_idx - placed_dense_cols) as u8);
                }
            }
        }

        full
    }

    /// Identity Matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::id();
    /// assert_eq!(result[0][0].0, 1);
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A new identity matrix
    pub fn id() -> Self {
        let mut a = Self::default();
        for i in 0..cmp::min(N, M) {
            a[i][i] = Fq(1);
        }
        a
    }

    /// sample a random matrix
    /// # Parameters
    /// - `seed`: a already initializzed XofReader
    pub fn rand<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut a = Self::default();
        for i in 0..N {
            rand_range_q_state_elements(&mut a.rows[i].0[..M], state);
        }
        a
    }

    /// sample a random matrix
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// ```
    ///
    /// # Parameters
    /// - `seed`: NxM matrix
    pub fn rand_from_seed<S>(seed: &[u8]) -> Self
    where 
        S: ExtendableOutput + Default + Clone
    {
        // static_assert(seed.len() % 32 == 0);
        let mut hasher = S::default();
        hasher.update(seed);
        let mut reader = hasher.finalize_xof();
        Self::rand(&mut reader)
    }

    /// compares two matrices
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let A: Matrix<100, 100, 100> = Matrix::default();
    /// let B: Matrix<100, 100, 100> = Matrix::default();
    /// A == B;
    /// A <= B;
    /// A >= B;
    /// A < B;
    /// A > B;
    /// ```
    ///
    /// # Returns
    ///  0 if a == b
    ///  x if a  > b
    /// -x if a  < b
    fn partial_cmp(a: &Self, b: &Self) -> i32 {
        for row in 0..a.nrows() {
            let mut i: usize = 0;
            while i < a.ncols() && a[row][i] == b[row][i] {
                i += 1;
            }

            if i == a.ncols() {
                continue;
            }

            return (b[row][i].0 as i32) - (a[row][i].0 as i32);
        }
        0
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let A: Matrix<100, 100, 100> = Matrix::default();
    /// let B: Matrix<100, 100, 100> = Matrix::default();
    /// A == B;
    /// ```
    ///
    /// # Returns
    ///  1 if a == b, else 0
    fn eq(a: &Self, b: &Self) -> bool {
        Self::partial_cmp(a, b) == 0
    }

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::default();
    /// assert_eq!(result.is_zero(), true);
    /// ```
    ///
    /// # Returns
    /// A size/dimension of the vector
    #[inline]
    #[must_use]
    pub fn is_zero(&self) -> bool {
        for row in 0..self.nrows() {
            for col in 0..self.ncols() {
                if self[row][col].0 > 0 {
                    return false
                }
            }
        }

        true
    }

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::default();
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A size/dimension of the vector
    #[inline]
    #[must_use]
    pub fn dimension(&self) -> (usize, usize) {
        (N, M)
    }

    /// Return the number of columns
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::default();
    /// assert_eq!(result.ncols(), 100);
    /// ```
    ///
    /// # Returns
    /// number of columns
    #[inline]
    #[must_use]
    pub fn ncols(&self) -> usize {
        M
    }

    /// Return the number of rows
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100, 100> = Matrix::default();
    /// assert_eq!(result.nrows(), 100);
    /// ```
    ///
    /// # Returns
    /// n_rows
    #[inline]
    #[must_use]
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
            Vector::<M, M_PAD>::add2(&mut r[i],&l[i]);
        }
    }

    /// constant time subtraction
    /// r -= l
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// const N: usize = 128;
    /// let mut C = Matrix::<N, N, N>::default();
    /// let A = Matrix::<N, N, N>::default();
    /// Matrix::sub(&mut C, &A)
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `l`: NxM matrix
    #[inline]
    pub fn sub(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<M, M_PAD>::sub2(&mut r[i],&l[i]);
        }
    }

    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let mut a = Matrix::<32, 32, 32>::default();
    /// a.swap_rows(10, 11);
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    #[inline]
    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        for j in 0..M {
            let tmp = self.rows[r1][j];
            self.rows[r1][j] = self.rows[r2][j];
            self.rows[r2][j] = tmp;
        } 
    }

    /// TODO not ct
    /// row swap
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// use less::monomial::Monomial;
    /// const N: usize = 32;
    /// let g = Matrix::<N, N, N>::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// let monom = Monomial::<N, N>::rand_from_seed::<Shake128>(&[1u8; 32]);
    /// let mut res = Matrix::<N, N, N>::default();
    /// Matrix::monomial_mul(&mut res, &g, &monom);
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn monomial_mul<const M_prime: usize, const M_prime_PAD: usize>(res: &mut Self, g: &Self, monom: &Monomial<M_prime, M_prime_PAD>) {
        // static_assert!(M <= M_prime);
        // `M_prime` is the length of the code, the padding is not permuted
        for src_col_idx in 0..M_prime {
            for row_idx in 0..N {
                let pos = monom.perms[src_col_idx];
                let a = g.rows[row_idx][src_col_idx];
                let b = monom.coeffs[src_col_idx];
                res.rows[row_idx][pos as usize] = Fq::mul(a, b);
            }
        }
    }

    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// const N: usize = 32;
    /// let mut a = Matrix::<N, N, N>::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// let mut is_pivot_column = [0u8; N];
    /// let mut was_pivot_column = [0u8; N];
    /// assert!(a.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, 0));
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    #[must_use]
    pub fn rref<const CT: bool>(&mut self,
                                is_pivot_column: &mut [u8; M],
                                was_pivot_column: &mut [u8; M],
                                pvt_reuse_limit: usize) -> bool {
        let mut pvt_reuse_cnt: usize = 0;

        // row swap pre-process - swap previous pivot elements to corresponding
        // row to reduce likelihood of corruption
        if pvt_reuse_limit != 0 {
            // `N` is the number of rows, the scan goes from the last row
            // down to the first one.
            for preproc_col in (0..N).rev() {
                if was_pivot_column[preproc_col as usize] == 1 {
                    let mut pivot_el_row = 0;
                    for row in 0..N {
                        if self.rows[row][preproc_col] != Fq(0) {
                            pivot_el_row = row;
                        }
                    }

                    Self::swap_rows(self, preproc_col, pivot_el_row);
                }
            }
        }

     for row_to_reduce in 0..N {
            let mut pivot_row = row_to_reduce;
            let mut pivot_column = row_to_reduce;

            /* Search pivot */
            while pivot_column < M && self.rows[pivot_row][pivot_column] == Fq::ZERO {
                while pivot_row < N && self.rows[pivot_row][pivot_column] == Fq::ZERO {
                    pivot_row += 1;
                }

                if pivot_row >= N {
                    pivot_column += 1;
                    pivot_row = row_to_reduce;
                }
            }

            if pivot_column >= M {
                return false;
            }

            is_pivot_column[pivot_column] = 1;

            /* Swap rows if needed */
            if row_to_reduce != pivot_row {
                // pivot no longer reusable
                was_pivot_column[pivot_row] = 0;
                Self::swap_rows(self, row_to_reduce, pivot_row);
            }

            pivot_row = row_to_reduce;

            /* Pivot reuse shortcut */
            if was_pivot_column[pivot_column] == 1
                && pvt_reuse_cnt < pvt_reuse_limit
                && pivot_column < N
            {
                pvt_reuse_cnt += 1;
                continue;
            }

            /* Rescale pivot row */
            let scaling_factor = match CT {
                true => Fq::inv(self[pivot_row][pivot_column]),
                false => Fq::inv_non_ct(self[pivot_row][pivot_column]),
            };

            // rescaling and elimination run over all `M` columns
            for col in pivot_column..M {
                self[pivot_row][col] = scaling_factor * self[pivot_row][col];
            }

            /* Eliminate pivot column from other rows */
            for row_idx in 0..N {
                if row_idx != pivot_row {
                    let multiplier = self[row_idx][pivot_column];
                    for col_idx in 0..M {
                        let tmp = multiplier * self[pivot_row][col_idx];
                        self[row_idx][col_idx] -= tmp;
                    }
                }
            }
        }

        true
    }

    /// translation of `UnpackCosetRep`
    /// NOTE: not constant time
    /// \param res[out]: G*c: a generator matrix: K \times N-K
    /// \param g[in]: current generator matrix: K \times N-K
    /// \param c[in]: compressed cf action
    /// \param initial_g_col_pivot[in]: input IS
    /// \param permuted_g_col_pivot[out]: output IS, to keep track of the pivot cols
    pub fn unpack_coset_rep(res: &mut Self,
                            g: &Self,
                            c: &[u8; N8],
                            initial_g_col_pivot: &[u8; M],
                            permuted_g_col_pivot: &mut [u8; M]) {
        let mut l = 0usize;
        let mut r = 0usize;

        'outer: for i in 0..N8 {
            for j in 0..8 {
                if (i*8 + j) >= M { break 'outer; }

                let bit = (c[i] >> j) & 1u8;
                let pos = if bit != 0 {
                    let p = l;
                    l += 1;
                    p
                } else {
                    let p = K + r;
                    r += 1;
                    p
                };

                permuted_g_col_pivot[pos] = initial_g_col_pivot[i*8 + j];

                // copy the column
                for k in 0..N {
                    res[k][pos] = g[k][i*8 + j];
                }
            }
        }
    }

    /// Compresses a generator matrix in RREF into a byte array
    ///
    /// * `compressed` – output buffer (RREF_MAT_PACKEDBYTES bytes)
    /// * `is_pivot_column` – length-N array, exactly K entries set to 1
    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// use less::config::{K, N, N_pad, RREF_MAT_PACKEDBYTES};
    /// let mut a = Matrix::<K, N, N_pad>::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// let mut is_pivot_column = [0u8; N];
    /// let mut was_pivot_column = [0u8; N];
    /// assert!(a.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, K));
    ///
    /// // an RREF generator has exactly K pivot columns, so the packed form
    /// // fills RREF_MAT_PACKEDBYTES exactly
    /// let mut compressed = [0u8; RREF_MAT_PACKEDBYTES];
    /// a.compress(&mut compressed, &is_pivot_column);
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn compress(
        &self,
        compressed: &mut [u8],
        is_pivot_column: &[u8; M],
    ) {
        /* Compress pivot flags */
        // the code has `crate::config::N` columns, the rest of the row is padding
        for col_byte in 0..(crate::config::N / 8) {
            compressed[col_byte] =
                (is_pivot_column[8 * col_byte + 0] << 0) |
                (is_pivot_column[8 * col_byte + 1] << 1) |
                (is_pivot_column[8 * col_byte + 2] << 2) |
                (is_pivot_column[8 * col_byte + 3] << 3) |
                (is_pivot_column[8 * col_byte + 4] << 4) |
                (is_pivot_column[8 * col_byte + 5] << 5) |
                (is_pivot_column[8 * col_byte + 6] << 6) |
                (is_pivot_column[8 * col_byte + 7] << 7);
        }

        /* Handle category-dependent tail */
        let mut compress_idx: usize;
        if crate::config::N == 252 || crate::config::N == 548 {
            compressed[crate::config::N / 8] =
                (is_pivot_column[crate::config::N - 4] << 0) |
                (is_pivot_column[crate::config::N - 3] << 1) |
                (is_pivot_column[crate::config::N - 2] << 2) |
                (is_pivot_column[crate::config::N - 1] << 3);

            compress_idx = crate::config::N / 8 + 1;
        } else {
            compress_idx = crate::config::N / 8;
        }

        /* Compress non-pivot columns row-by-row */
        let mut encode_state: u8 = 0;

        for row_idx in 0..N {
            for col_idx in 0..crate::config::N {
                if is_pivot_column[col_idx] == 0 {
                    let val: u8 = self[row_idx][col_idx].0;

                    match encode_state {
                        0 => {
                            compressed[compress_idx] = val;
                        }
                        1 => {
                            compressed[compress_idx] |= val << 7;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 1;
                        }
                        2 => {
                            compressed[compress_idx] |= val << 6;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 2;
                        }
                        3 => {
                            compressed[compress_idx] |= val << 5;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 3;
                        }
                        4 => {
                            compressed[compress_idx] |= val << 4;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 4;
                        }
                        5 => {
                            compressed[compress_idx] |= val << 3;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 5;
                        }
                        6 => {
                            compressed[compress_idx] |= val << 2;
                            compress_idx += 1;
                            compressed[compress_idx] = val >> 6;
                        }
                        7 => {
                            compressed[compress_idx] |= val << 1;
                            compress_idx += 1;
                        }
                        _ => unreachable!(),
                    }

                    if encode_state != 7 {
                        encode_state += 1;
                    } else {
                        encode_state = 0;
                    }
                }
            }
        }
    }

    /// Expands a compressed RREF generator matrix into a full one
    ///
    /// * `compressed` – input byte stream
    /// * `is_pivot_column` – length-N array, will be initialized here
    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::Matrix;
    /// use less::config::{K, N, N_pad, RREF_MAT_PACKEDBYTES};
    /// let mut a = Matrix::<K, N, N_pad>::rand_from_seed::<Shake128>(&[0u8; 32]);
    /// let mut is_pivot_column = [0u8; N];
    /// let mut was_pivot_column = [0u8; N];
    /// assert!(a.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, K));
    ///
    /// let mut compressed = [0u8; RREF_MAT_PACKEDBYTES];
    /// a.compress(&mut compressed, &is_pivot_column);
    ///
    /// // expand is the exact inverse of compress on an RREF matrix
    /// let mut b = Matrix::<K, N, N_pad>::default();
    /// let mut flags = [0u8; N];
    /// b.expand(&compressed, &mut flags);
    /// assert_eq!(flags, is_pivot_column);
    /// for row in 0..K {
    ///     for col in 0..N {
    ///         assert_eq!(b[row][col], a[row][col]);
    ///     }
    /// }
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn expand(
        &mut self,
        compressed: &[u8],
        is_pivot_column: &mut [u8; M]) {
        /* Initialize pivot flags */
        for i in 0..crate::config::N {
            is_pivot_column[i] = 0;
        }

        /* Decompress pivot flags */
        for col_byte in 0..(crate::config::N / 8) {
            let b = compressed[col_byte];
            is_pivot_column[col_byte * 8 + 0] = (b >> 0) & 0x1;
            is_pivot_column[col_byte * 8 + 1] = (b >> 1) & 0x1;
            is_pivot_column[col_byte * 8 + 2] = (b >> 2) & 0x1;
            is_pivot_column[col_byte * 8 + 3] = (b >> 3) & 0x1;
            is_pivot_column[col_byte * 8 + 4] = (b >> 4) & 0x1;
            is_pivot_column[col_byte * 8 + 5] = (b >> 5) & 0x1;
            is_pivot_column[col_byte * 8 + 6] = (b >> 6) & 0x1;
            is_pivot_column[col_byte * 8 + 7] = (b >> 7) & 0x1;
        }

        /* Category-dependent tail */
        let mut compress_idx: usize;
        if crate::config::N == 252 || crate::config::N == 548 {
            let b = compressed[crate::config::N / 8];
            is_pivot_column[crate::config::N - 4] = (b >> 0) & 0x1;
            is_pivot_column[crate::config::N - 3] = (b >> 1) & 0x1;
            is_pivot_column[crate::config::N - 2] = (b >> 2) & 0x1;
            is_pivot_column[crate::config::N - 1] = (b >> 3) & 0x1;

            compress_idx = crate::config::N / 8 + 1;
        } else {
            compress_idx = crate::config::N / 8;
        }

        /* Decompress matrix row-by-row */
        let mut decode_state: u8 = 0;

        for row_idx in 0..N {
            let mut pivot_idx: usize = 0;

            for col_idx in 0..crate::config::N {
                if is_pivot_column[col_idx] == 0 {
                    let val: u8 = match decode_state {
                        0 => {
                            compressed[compress_idx] & Fq::Q_MASK
                        }
                        1 => {
                            let v = ((compressed[compress_idx] >> 7)
                                | (compressed[compress_idx + 1] << 1))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        2 => {
                            let v = ((compressed[compress_idx] >> 6)
                                | (compressed[compress_idx + 1] << 2))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        3 => {
                            let v = ((compressed[compress_idx] >> 5)
                                | (compressed[compress_idx + 1] << 3))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        4 => {
                            let v = ((compressed[compress_idx] >> 4)
                                | (compressed[compress_idx + 1] << 4))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        5 => {
                            let v = ((compressed[compress_idx] >> 3)
                                | (compressed[compress_idx + 1] << 5))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        6 => {
                            let v = ((compressed[compress_idx] >> 2)
                                | (compressed[compress_idx + 1] << 6))
                                & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        7 => {
                            let v = (compressed[compress_idx] >> 1) & Fq::Q_MASK;
                            compress_idx += 1;
                            v
                        }
                        _ => unreachable!(),
                    };

                    self[row_idx][col_idx] = Fq(val);

                    if decode_state != 7 {
                        decode_state += 1;
                    } else {
                        decode_state = 0;
                    }
                } else {
                    /* Pivot column */
                    let v: u8 = if row_idx == pivot_idx { 1 } else { 0 };
                    self[row_idx][col_idx] = Fq(v);
                    pivot_idx += 1;
                }
            }
        }
    }
}


impl<const N: usize, const M: usize, const M_PAD: usize> Add for Matrix<N, M, M_PAD> {
    type Output = Matrix<N, M, M_PAD>;
    fn add(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::add(&mut t, &r);
        t
    }
}
impl<'a, const N: usize, const M: usize, const M_PAD: usize> Add for &'a mut Matrix<N, M, M_PAD> {
    type Output = &'a mut Matrix<N, M, M_PAD>;
    fn add(self, r: Self) -> Self::Output {
        Matrix::<N, M, M_PAD>::add(self, r);
        self
    }
}
impl<const N: usize, const M: usize, const M_PAD: usize> Sub for Matrix<N, M, M_PAD> {
    type Output = Matrix<N, M, M_PAD>;
    fn sub(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::sub(&mut t, &r);
        t
    }
}
impl<'a, const N: usize, const M: usize, const M_PAD: usize> Sub for &'a mut Matrix<N, M, M_PAD> {
    type Output = &'a mut Matrix<N, M, M_PAD>;
    fn sub(self, r: Self) -> Self::Output {
        Matrix::<N, M, M_PAD>::sub(self, r);
        self
    }
}
impl<const N: usize, const M: usize, const M_PAD: usize> Mul<Monomial<M, M_PAD>> for Matrix<N, M, M_PAD> {
    type Output = Matrix<N, M, M_PAD>;
    fn mul(self, r: Monomial<M, M_PAD>) -> Self::Output {
        let mut t = Matrix::default();
        Self::monomial_mul(&mut t, &self, &r);
        t
    }
}


/// Direct access to a row
impl<const N: usize, const M: usize, const M_PAD: usize> Index<usize> for Matrix<N, M, M_PAD> {
    type Output = Vector<M, M_PAD>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}
/// Direct access to a row (mutable)
impl<const N: usize, const M: usize, const M_PAD: usize> IndexMut<usize> for Matrix<N, M, M_PAD> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.rows[index]
    }
}
impl<const N: usize, const M: usize, const M_PAD: usize> fmt::Debug for Matrix<N, M, M_PAD> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:?}", self.rows)
    }
}

impl<const N: usize, const M: usize, const M_PAD: usize> Display for Matrix<N, M, M_PAD> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for i in 0..N {
            write!(f, "{}", self[i])?;
        }

        Ok(())
    }
}


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct MatrixNormalized<const N: usize, const N_PAD: usize> {
    pub rows: [Vector<N, N_PAD>; N_PAD],
}

impl<const N: usize, const N_PAD: usize> MatrixNormalized<N, N_PAD> {
    /// generates a new matrix init with 0
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::default();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    pub fn default() -> MatrixNormalized<N, N_PAD> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N, N_PAD>::default() ),
        }
    }
   
    /// generates a new matrix init with q-1
    pub fn default_large() -> MatrixNormalized<N, N_PAD> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N, N_PAD>::from_u8(Q-1) ),
        }
    }

    /// same as init, but inits all values with `v`
    #[inline]
    pub fn from_u8(v: u8) -> MatrixNormalized<N, N_PAD> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N, N_PAD>::from_u8(v) ),
        }
    }

    /// translation of 'void normalized_copy_from_generator_non_information_set'
    /// \param A[out]: pointer to allocated normalized struct, which get filled with the
    ///     non-IS of the generator matrix G
    /// \param G[in]: generator matrix to extract the non-IS from.
    /// \param is_pivot_column[in]: array identifying a pivot column via a 1
    #[inline]
    pub fn from_information_set<const N_PRIME: usize, const M_PRIME: usize, const M_PRIME_PAD: usize>(g: &Matrix<N_PRIME, M_PRIME, M_PRIME_PAD>, is_pivot_column: &[u8; M_PRIME]) -> MatrixNormalized<N, N_PAD> {
        let mut ret = MatrixNormalized::default();
        // we simply copy the last N-K columns even if they are not the information set.
        for i in 0..K {
            for j in 0..N {
                ret[i][j] = g[i][K + j];
            }
        }
        // now we scan if we need to fix the non information set
        let mut ctr = 0;
        while ctr < K && is_pivot_column[ctr] == 1 {
            ctr += 1;
        }

        // easy part: the last N-K columns are the non IS
        if ctr == K { return ret; }

        // "hard" part: copy all remaining columns < K into the non information set part
        let mut j = 0;
        while j < N && is_pivot_column[ctr] == 0 {
            // copy column
            for k in 0..K {
                ret[k][j] = g[k][ctr];
            }
            ctr += 1;
            j += 1;
        }

        ret
    }

    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let A: MatrixNormalized<100, 100> = MatrixNormalized::id();
    /// let B: MatrixNormalized<100, 100> = MatrixNormalized::id();
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
    fn partial_cmp(a: &Self, b: &Self) -> i32 {
        for row in 0..a.nrows() {
            let mut i: usize = 0;
            while i < a.ncols() && a[row][i] == b[row][i] {
                i += 1;
            }

            if i == a.ncols() {
                continue;
            }

            return (b[row][i].0 as i32) - (a[row][i].0 as i32);
        }
        0
    }

    /// checks for equality
    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let A: MatrixNormalized<100, 100> = MatrixNormalized::id();
    /// let B: MatrixNormalized<100, 100> = MatrixNormalized::id();
    /// A == B;
    /// ```
    ///
    /// # Returns
    ///  1 if a == b, else 0
    fn eq(a: &Self, b: &Self) -> bool {
        Self::partial_cmp(a, b) == 0
    }

    /// Identity Matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::id();
    /// assert_eq!(result[0][0].0, 1);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// A new identity matrix
    pub fn id() -> Self {
        let mut a = Self::default();
        for i in 0..N {
            a[i][i] = Fq(1);
        }
        a
    }

    /// sample a random matrix
    ///
    /// # Parameters
    /// - `state`: already init Xof
    pub fn rand<S>(state: &mut S) -> Self
    where
        S: XofReader
    {
        let mut a = Self::default();
        for row in 0..a.nrows() {
            a[row] = Vector::rand(state);
        }
        a
    }

    /// sample a random matrix
    /// # Examples
    ///
    /// ```
    /// use sha3::Shake128;
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::rand_from_seed::<Shake128>(&[0u8; 32]);
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

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::default();
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    /// # Returns
    /// A size/dimension of the vector
    #[inline]
    pub fn dimension(&self) -> usize {
        N
    }

    /// Return the number of columns
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::default();
    /// assert_eq!(result.ncols(), 100);
    /// ```
    ///
    /// # Returns
    /// number of columns
    #[inline]
    pub fn ncols(&self) -> usize {
        N
    }

    /// Return the number of rows
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100, 100> = MatrixNormalized::default();
    /// assert_eq!(result.nrows(), 100);
    /// ```
    ///
    /// # Returns
    /// n_rows
    #[inline]
    pub fn nrows(&self) -> usize {
        N
    }

    /// matrix addition
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// const N: usize = 128;
    /// let mut C: MatrixNormalized<N, N> = MatrixNormalized::default();
    /// let A : MatrixNormalized<N, N> = MatrixNormalized::default();
    /// MatrixNormalized::add(&mut C, &A);
    /// ```
    ///
    pub fn add(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<N, N_PAD>::add2(&mut r[i], &l[i]);
        }
    }

    /// matrix subtraction
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// const N: usize = 128;
    /// let mut C: MatrixNormalized<N, N> = MatrixNormalized::default();
    /// let A : MatrixNormalized<N, N> = MatrixNormalized::default();
    /// MatrixNormalized::add(&mut C, &A);
    /// ```
    ///
    pub fn sub(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<N, N_PAD>::sub2(&mut r[i], &l[i]);
        }
    }

    /// matrix scalar multiplication
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// use less::fq::Fq;
    /// const N: usize = 128;
    /// let mut C: MatrixNormalized<N, N> = MatrixNormalized::default();
    /// MatrixNormalized::scalar(&mut C, Fq(1));
    /// ```
    ///
    pub fn scalar(r: &mut Self, l: Fq) {
        for i in 0..N {
            Vector::<N, N_PAD>::scalar2(&mut r[i], l);
        }
    }

    /// row addition
    /// rows[i] += rows[j]
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut result: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// MatrixNormalized::row_add(&mut result, 0, 1);
    /// ```
    ///
    pub fn row_add(&mut self, i: usize, j: usize) {
        // what we want to implement
        // Vector::<N, N_PAD>::add2(&mut self.rows[i], &mut self.rows[j]);
        if i == j {
            return;
        }

        let [a, b] = self.rows.get_disjoint_mut([i, j]).unwrap();
        Vector::<N, N_PAD>::add2(a, b);
    }

    /// row subtraction
    /// rows[i] -= rows[j]
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut result: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// MatrixNormalized::row_sub(&mut result, 0, 1);
    /// ```
    ///
    pub fn row_sub(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        let [a, b] = self.rows.get_disjoint_mut([i, j]).unwrap();
        Vector::<N, N_PAD>::sub2(a, b);
    }

    /// row scalar multiplication
    /// rows[i] *= s 
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// use less::fq::Fq;
    /// let mut result: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// MatrixNormalized::row_scalar(&mut result, 0, Fq(1));
    /// ```
    ///
    pub fn row_scalar(&mut self, i: usize, s: Fq) {
        Vector::<N, N_PAD>::scalar2(&mut self.rows[i], s);
    }

    /// transpose: a = b^T
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut A: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// MatrixNormalized::transpose(&mut B, &A);
    /// ```
    pub fn transpose(a: &mut Self, b: &Self) {
        // TODO optimize
        for row in 0..a.ncols() {
            for col in 0..b.ncols() {
                a[col][row] = b[row][col];
            }
        }
    }

    /// NOTE: absolutely not constant time.
    /// NOTE: only operates on ptrs
    /// NOTE: not constant time
    /// \param G[in/out]: generator matrix to sort
    /// \param n[in] number of elements to sort
    /// \param L[in]: pointer to the currently shortest row
    /// \return 1 on success
    ///			0 if two rows generate the same multiset
    fn SortRows(&mut self, z: usize, L: &mut Multiset<Q_PAD>) -> i32 {
        // first sort each row into a tmp buffer
        let tmp: [Multiset<Q_PAD>; N_PAD] = Multiset::<Q_PAD>::from_matrix::<N, N_PAD>(self);

        // the sorting is done on a list of pointers into `tmp` alongside the
        // permutation `P`, both are swapped together, so `ptr[i] == tmp[P[i]]`
        // holds at any time and `P` alone is enough here.
        let mut p: [usize; N_PAD] = core::array::from_fn(|i| i);

        let mut max_zeros = 0u8;
        for i in 0..z {
            if tmp[i][0] > max_zeros { max_zeros = tmp[i][0]; }
        }

        // early exit if the shortest row of this matrix is already longer
        // than the currently shortest one
        if max_zeros < L[0] { return 0; }

        // quick sort with an explicit stack
        let mut stack = [0i64; 128];
        let mut s: i64 = -1;
        s += 1; stack[s as usize] = 0;
        s += 1; stack[s as usize] = (z as i64) - 1;

        while s >= 0 {
            let h = stack[s as usize]; s -= 1;
            let l = stack[s as usize]; s -= 1;

            // hoare partition
            let pivot = p[h as usize];
            let mut i: i64 = l - 1;
            for j in l..h {
                if Multiset::<Q_PAD>::partial_cmp(&tmp[p[j as usize]], &tmp[pivot]) < 0 {
                    i += 1;
                    if i == j { continue; }
                    p.swap(i as usize, j as usize);
                }
            }
            if i + 1 != h {
                p.swap((i + 1) as usize, h as usize);
            }
            let q = i + 1;

            if q - 1 > l {
                s += 1; stack[s as usize] = l;
                s += 1; stack[s as usize] = q - 1;
            }
            if q + 1 < h {
                s += 1; stack[s as usize] = q + 1;
                s += 1; stack[s as usize] = h;
            }
        }

        // apply the permutation
        for t in 0..z {
            let mut ind = p[t];
            while ind < t { ind = p[ind]; }

            if t != ind {
                self.rows.swap(t, ind);
            }
        }

        1
    }

    /// TODO example and test
    fn SortCols(&mut self, z: usize) -> i32 {
        let _ = z;
        let mut VT = MatrixNormalized::default();
        MatrixNormalized::transpose(&mut VT, &self);

        // same pointer/permutation trick as in `SortRows`
        let mut p: [usize; N_PAD] = core::array::from_fn(|i| i);

        // quick sort with an explicit stack
        let mut stack = [0i64; 64];
        let mut s: i64 = -1;
        s += 1; stack[s as usize] = 0;
        s += 1; stack[s as usize] = (N as i64) - 1;

        while s >= 0 {
            let h = stack[s as usize]; s -= 1;
            let l = stack[s as usize]; s -= 1;

            // hoare partition
            let pivot = p[h as usize];
            let mut i: i64 = l - 1;
            for j in l..h {
                if Vector::cmp(&VT[p[j as usize]], &VT[pivot]) > 0 {
                    i += 1;
                    if i == j { continue; }
                    p.swap(i as usize, j as usize);
                }
            }
            if i + 1 != h {
                p.swap((i + 1) as usize, h as usize);
            }
            let q = i + 1;

            if q - 1 > l {
                s += 1; stack[s as usize] = l;
                s += 1; stack[s as usize] = q - 1;
            }
            if q + 1 < h {
                s += 1; stack[s as usize] = q + 1;
                s += 1; stack[s as usize] = h;
            }
        }

        // apply the permutation
        for t in 0..N {
            let mut ind = p[t];
            while ind < t { ind = p[ind]; }

            if t != ind {
                VT.rows.swap(t, ind);
            }
        }

        MatrixNormalized::transpose(self, &VT);
        1
    }

    /// by floyd
    /// TODO example and test
    fn scale_cf_preprocess_v1(&mut self, z: usize, m: &Multiset<Q_PAD>) -> i32 {
        let mut tmp = Multiset::<Q_PAD>::default();
        for i in 0..z  {
            let mut s = Vector::<N, N_PAD>::acc(&self.rows[i]);
            if s != Fq(0) {
                s = Fq::inv_non_ct(s);
            } else {
                s = Vector::<N, N_PAD>::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    return 0;
                }
            }

            self.row_scalar(i as usize, s);
            // compute the histogram
            Multiset::<Q_PAD>::from_row(&mut tmp, &self.rows[i]);
            if Multiset::<Q_PAD>::partial_cmp(&tmp, m) < 0 {
                return 1;
            }
        };

        0
    }

    /// by luke
    /// TODO example and test
    fn ScaleCFSubPreprocess(&mut self, z: usize, m: &mut Multiset<Q_PAD>) -> i32 {
        let mut tmp = Multiset::<Q_PAD>::default();
        let mut ret = 0;
        for i in 0..z  {
            let mut s = Vector::<N, N_PAD>::acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv_non_ct(s);
            } else {
                s = Vector::<N, N_PAD>::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_scalar(i as usize, s);
            // compute the histogram
            Multiset::<Q_PAD>::from_row(&mut tmp, &self.rows[i]);
            if Multiset::<Q_PAD>::partial_cmp(&tmp, m) < 0  {
                ret = 1;
                // copy new smallest row
                for i in 0..Q_PAD {
                    m[i] = tmp[i];
                }
            }
        };

        ret
    }

    /// Implementation of the `SortCF` algorithm from the specification.
    /// NOTE: non-constant time
    /// NOTE: computes the result inplace
    /// NOTE: early exist from the `SortRows` function if `L` is already shorter.
    /// NOTE: first sort the rows, then the columns
    /// \param G[in/out]: pointer to the non IS-part of a generator matrix.
    /// The rows and then the columns are sorted inplace.
    /// \param L[in]: pointer the currently shortest row.
    /// \return 0 on failure (identical rows, which create the same multiset)
    /// 		1 on success
    fn SortCF(&mut self, L: &mut Multiset<Q_PAD>) -> i32 {
        if self.SortRows(N, L) == 0 {
            return 0;
        }
        self.SortCols(N_PAD); // NOTE: must be N_PAD
        1
    }

    /// Implementation of the `ScaleCF` algorithm from the specification.
    /// NOTE: non-constant time
    /// NOTE: computes the result inplace
    /// NOTE: does not early exit. Only `SortCF` can do that.
    /// \param G[in/out]: pointer to the non IS-part of a generator matrix.
    /// \param L[in]: pointer the currently shortest row, in histogram form.
    /// \return 0 on failure:
    /// 			- compute_power_column fails.
    /// 			- identical rows, which create the same multiset
    /// 		1 on success
    fn ScaleCF(&mut self, L: &mut Multiset<Q_PAD>) -> i32 {
        for row in 0..N {
            if Vector::all_same(&self[row]) {
                continue;
            }

            let mut s = Vector::acc(&self[row]);
            if s.0 != 0 {
                s = Fq::inv_non_ct(s);
            } else {
                s = Vector::<N, N_PAD>::acc_inv(&self.rows[row]);
                if s.0 == 0 {
                    return 0;
                }
            }
            Vector::scalar2(&mut self[row], s);
        }
        self.SortCF(L)
    }

    //pub fn cf_opt1(&mut self) {
    //    let mut M = Self::new_large();
    //    let mut J = [0; N];
    //    let mut z = 0;
    //    let mut num_zeros = 0;
    //}

    /// res = g * monom
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut A: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// // TODO
    /// ```
    pub fn monomial_mul(res: &mut Self, g: &Self, monom: &Monomial<N, N_PAD>) {
        for src_col_idx in 0..N {
            for row_idx in 0..N {
                let pos = monom.perms[src_col_idx];
                let a = g.rows[row_idx][src_col_idx];
                let b = monom.coeffs[src_col_idx];
                res.rows[row_idx][pos as usize] = Fq::mul(a, b);
            }
        }
    }

    /// res = g * monom
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut A: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64, 64> = MatrixNormalized::default();
    /// // TODO
    /// ```
    pub fn permutation_mul(res: &mut Self, g: &Self, perm: &Permutation<N>) {
        for src_col_idx in 0..N {
            for row_idx in 0..N {
                let pos = perm.perms[src_col_idx];
                res.rows[row_idx][pos as usize] = g.rows[row_idx][src_col_idx];
            }
        }
    }

    /// TODO doc
    pub fn blind<S>(&mut self, state: &mut S)
    where
        S: XofReader
    {
        // the randomness for `right` is drawn before the one for `left`
        let right = Monomial::<N, N_PAD>::rand(state);
        let left = Monomial::<N, N_PAD>::rand(state);

        // right multiplication
        let mut tmp = MatrixNormalized::default();
        Self::monomial_mul(&mut tmp, self, &right);

        // left multiplication
        for i in 0..N {
            let pos = left.perms[i] as usize;
            let val = left.coeffs[i];
            Vector::scalar(&mut self[i], &tmp[pos], val);
        }
    }

    /// NOTE: non-constant time
    /// NOTE: computes the result inplace
    /// This is the second-fastest implementation of the canonical form function.
    /// Fallback implementation if `cf` fails
    /// \param G[in/out] non IS part of a generator matrix
    /// \return 0 on failure
    /// 		1 on success
    pub fn ImprovedCFBase(&mut self) -> i32 {
        let mut M = MatrixNormalized::<N, N_PAD>::default_large();
        let mut touched = 0i32;
        let mut J = [0u32; N];
        let mut z = 0usize;

        // count zeros in each row;
        let mut max_zeros = 0u32;
        let mut Z = [false; N];

        for row in 0..N {
            let num_zeros = Vector::count_zero(&self[row]);
            Z[row] = num_zeros > 0;
            if num_zeros < max_zeros {
                continue
            }
            if num_zeros > max_zeros {
                z = 1;                      // reset z
                J[0] = row as u32;
                max_zeros = num_zeros;
                continue;
            }
            // num_zeros == max_zeros
            J[z] = row as u32;
            z += 1;
        }

        // NOTE: fallback solution if everything falls apart.
        // NOTE: reaching this needs every row to contain the same number of
        // zeros, `CFOriginal` is not part of this port.
        if z == N {
            return 0;
        }

        let mut B = MatrixNormalized::<N, N_PAD>::default();
        let mut row_inv_data = Vector::<N, N_PAD>::default();
        // NOTE: this is already "sorted"
        let mut L = Multiset::<Q_PAD>::default();

        // Check smallest rows of all matrices to find the smallest candidate
        for row in 0..N {
            if Z[row] { continue; }

            Vector::inv_non_ct(&mut row_inv_data, &self[row]);
            for row2 in 0..z {
                Vector::mul(&mut B[row2], &self[J[row2] as usize], &row_inv_data);
            }

            if B.scale_cf_preprocess_v1(z, &L) == 1 {
                for row2 in 0..N {
                    Vector::mul(&mut B[row2], &self[row2], &row_inv_data);
                }

                let ret = B.ScaleCF(&mut L);
                // NOTE: `compare_matrices` has the opposite sign convention
                // of `partial_cmp`, so this is `compare_matrices(B, M) < 0`
                if ret == 1 && MatrixNormalized::partial_cmp(&B, &M) > 0 {
                    Multiset::<Q_PAD>::from_row(&mut L, &B[0]);
                    touched = 1;
                    for row2 in 0..N {
                        M[row2].0 = B[row2].0;
                    }
                }
            }
        }

        for row in 0..N {
            self[row].0 = M[row].0;
        }

        touched
    }

    /// translation of `ImprovedCF`
    /// NOTE: non-constant time
    /// NOTE: computes the result inplace
    /// This is the fastest implementation of the canonical form function.
    /// This function scans the matrix for the rows which probably lead
    /// to the shortest canonical form. If this process fails it falls
    /// back to either `CFOriginal` or `ImprovedCF`. The very slow `CFOriginal`
    /// is only called if every row in the matrix contains zeros.
    /// \param G[in/out] non IS part of a generator matrix
    /// \return 0 on failure
    /// 		1 on success/
    pub fn cf(&mut self) -> i32 {
        let mut J = [0u32; N];
        let mut z = 0usize;
        let mut smallest_scaling_row = 0usize;

        // count zeros in each row;
        let mut max_zeros = 0u32;
        let mut Z = [false; N];

        for row in 0..N {
            let num_zeros = Vector::count_zero(&self[row]);
            Z[row] = num_zeros > 0;
            if num_zeros < max_zeros {
                continue
            }
            if num_zeros > max_zeros {
                z = 1;                      // reset z
                J[0] = row as u32;
                max_zeros = num_zeros;
                continue;
            }
            // num_zeros == max_zeros
            J[z] = row as u32;
            z += 1;
        }

        // NOTE: fallback solution if everything falls apart, see `ImprovedCFBase`
        if z == N {
            return 0;
        }

        let mut B = MatrixNormalized::<N, N_PAD>::default();
        let mut row_inv_data = Vector::<N, N_PAD>::default();
        // NOTE: this is already "sorted"
        let mut L = Multiset::<Q_PAD>::default();

        // Check smallest rows of all matricies to find the smallest candidate
        for row in 0..N {
            if Z[row] { continue; }

            Vector::inv_non_ct(&mut row_inv_data, &self[row]);
            for row2 in 0..z {
                Vector::mul(&mut B[row2], &self[J[row2] as usize], &row_inv_data);
            }

            if B.ScaleCFSubPreprocess(z, &mut L) == 1 {
                smallest_scaling_row = row;
            }
        }

        // Calculate CF for best candidate
        Vector::inv_non_ct(&mut row_inv_data, &self[smallest_scaling_row]);
        for row2 in 0..N {
            Vector::mul(&mut B[row2], &self[row2], &row_inv_data);
        }
        let ret = B.ScaleCF(&mut L);
        // If candidate was not valid, fall back to regular approach
        if ret != 1 {
            // Best scaled row was not valid;
            return self.ImprovedCFBase();
        }

        for row in 0..N {
            self[row].0 = B[row].0;
        }
        ret
    }
}

impl<const N: usize, const N_PAD: usize> Add for MatrixNormalized<N, N_PAD> {
    type Output = MatrixNormalized<N, N_PAD>;
    fn add(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::add(&mut t, &r);
        t
    }
}
impl<'a, const N: usize, const N_PAD: usize> Add for &'a mut MatrixNormalized<N, N_PAD> {
    type Output = &'a mut MatrixNormalized<N, N_PAD>;
    fn add(self, r: Self) -> Self::Output {
        MatrixNormalized::<N, N_PAD>::add(self, r);
        self
    }
}
impl<const N: usize, const N_PAD: usize> Sub for MatrixNormalized<N, N_PAD> {
    type Output = MatrixNormalized<N, N_PAD>;
    fn sub(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::sub(&mut t, &r);
        t
    }
}
impl<'a, const N: usize, const N_PAD: usize> Sub for &'a mut MatrixNormalized<N, N_PAD> {
    type Output = &'a mut MatrixNormalized<N, N_PAD>;
    fn sub(self, r: Self) -> Self::Output {
        MatrixNormalized::<N, N_PAD>::sub(self, r);
        self
    }
}

impl<const N: usize, const N_PAD: usize> Mul<Permutation<N>> for MatrixNormalized<N, N_PAD> {
    type Output = MatrixNormalized<N, N_PAD>;
    fn mul(self, r: Permutation<N>) -> Self::Output {
        let mut t = MatrixNormalized::default();
        Self::permutation_mul(&mut t, &self, &r);
        t
    }
}

/// Direct access to a row
impl<const N: usize, const N_PAD: usize> Index<usize> for MatrixNormalized<N, N_PAD> {
    type Output = Vector<N, N_PAD>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}
/// Direct access to a row (mutable)
impl<const N: usize, const N_PAD: usize> IndexMut<usize> for MatrixNormalized<N, N_PAD> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.rows[index]
    }
}


/// use less::matrix::Matrix;
#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 32;
    const M: usize = 32;
    
    #[test]
    fn init() {
        let _ = Matrix::<N, M, M>::default();
    }
    #[test]
    fn from_u8() {
        let a = Matrix::<N, M, M>::from_u8(1);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert_eq!(a[row][col].0, 1);
            }
        }
    }

    #[test]
    fn rand() {
        let a = Matrix::<N, M, M>::rand_from_seed::<Shake128>(&[0u8; 32]);
        println!("{}", a);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert!(a[row][col].0 <= Fq::Q);
            }
        }

        assert_eq!(a.is_zero(), false);
    }
    #[test]
    fn rand_from_seed() {
        let seed = [0u8; 32];
        let a = Matrix::<N, M, M>::rand_from_seed::<Shake128>(&seed);
        println!("{}", a);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert!(a[row][col].0 <= Fq::Q);
            }
        }

        assert_eq!(a.is_zero(), false);
    }

    #[test]
    fn id() {
        let a = Matrix::<N, M, M>::id();
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                let b = (row == col) as u8;
                assert_eq!(a[row][col].0, b);
            }
        }
    }
    #[test]
    fn ncols() {
        let a = Matrix::<N, M, M>::default();
        assert_eq!(a.ncols(), M);
    }
    #[test]
    fn nrows() {
        let a = Matrix::<N, M, M>::default();
        assert_eq!(a.nrows(), N);
    }

    #[test]
    fn add() {
        let a = Matrix::<N, M, M>::default();
        let b = Matrix::<N, M, M>::from_u8(1);
        let c = a.clone().add(b.clone());
        for i in 0..N {
            assert_eq!(c[i][i].0, 1);
        }
        let c = a + b;
        for i in 0..N {
            assert_eq!(c[i][i].0, 1);
        }
    }

    #[test]
    fn sub() {
        let a = Matrix::<N, M, M>::default();
        let b = Matrix::<N, M, M>::from_u8(1);
        let c = b - a;
        for i in 0..N {
            assert_eq!(c[i][i].0, 1);
        }
    }

    /// res = g*monom permutes the columns of `g` and scales them
    #[test]
    fn monomial_mul() {
        use crate::prng::csprng_initialize;
        use crate::config::Xof;

        let g = Matrix::<N, M, M>::rand_from_seed::<Xof>(&[0u8; 32]);
        let monom = Monomial::<M, M>::rand_from_seed::<Xof>(&[1u8; 32]);

        let mut res = Matrix::<N, M, M>::default();
        Matrix::monomial_mul(&mut res, &g, &monom);

        // column `col` of `g` lands on column `monom.perms[col]`, scaled by
        // `monom.coeffs[col]`
        for col in 0..M {
            let dst = monom.perms[col] as usize;
            for row in 0..N {
                assert_eq!(res[row][dst], Fq::mul(g[row][col], monom.coeffs[col]),
                           "row {} col {} -> {}", row, col, dst);
            }
        }

        // multiplying with the identity monomial only scales
        let mut id = Monomial::<M, M>::default();
        for i in 0..M { id.coeffs[i] = Fq(1); }
        let mut res = Matrix::<N, M, M>::default();
        Matrix::monomial_mul(&mut res, &g, &id);
        for row in 0..N {
            for col in 0..M {
                assert_eq!(res[row][col], g[row][col]);
            }
        }

        // and the result is a permutation of the columns, nothing is lost
        let mut s = csprng_initialize::<Xof>(&[2u8; 32]);
        let monom = Monomial::<M, M>::rand(&mut s);
        let mut res = Matrix::<N, M, M>::default();
        Matrix::monomial_mul(&mut res, &g, &monom);
        let mut seen = [false; M];
        for col in 0..M {
            let dst = monom.perms[col] as usize;
            assert!(!seen[dst]);
            seen[dst] = true;
        }
    }

    #[test]
    fn swap_rows() {
        let mut a = Matrix::<N, M, M>::default();
        let i1 = 0; let i2 = N/2;
        for i in 0..a.ncols() {
            a[i1][i] = Fq(1);
        }
    
        a.swap_rows(i1, i2);
        for i in 0..N {
            assert_eq!(a[i1][i].0, 0);
            assert_eq!(a[i2][i].0, 1);
        }
    }

    #[test]
    fn rref() {
        let mut a = Matrix::<N, M, M>::rand_from_seed::<Shake128>(&[0u8; 32]);
        let mut is_pivot_column = [0u8; M];
        let mut was_pivot_column = [0u8; M];
        let b = a.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, M);
        assert!(b);
    }
}
#[cfg(test)]
mod normalized_tests {
    use super::*;
    const N: usize = 32;

    #[test]
    fn init() {
        let _ = MatrixNormalized::<N, N>::default();
    }

    #[test]
    fn from_u8() {
        let a = MatrixNormalized::<N, N>::from_u8(1);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert_eq!(a[row][col], Fq(1));
            }
        }
    }

    #[test]
    fn transpose() {
        let a = MatrixNormalized::<N, N>::rand_from_seed::<Shake128>(&[0u8; 32]);
        let mut b = MatrixNormalized::<N, N>::default();
        MatrixNormalized::transpose(&mut b, &a);

        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert_eq!(a[col][row], b[row][col]);
            }
        }
    }

    #[test]
    fn blind() {
        let mut hasher = Shake128::default();
        hasher.update(b"123");
        let mut reader = hasher.finalize_xof();

        let mut a = MatrixNormalized::<N, N>::default();
        MatrixNormalized::blind(&mut a, &mut reader);

    }

    /// a matrix whose rows all have a distinct histogram
    fn distinct_rows() -> MatrixNormalized<N, N> {
        let mut a = MatrixNormalized::<N, N>::default();
        for row in 0..N {
            for col in 0..N {
                a[row][col] = Fq(if col <= row { 1 } else { 2 });
            }
        }
        a
    }

    #[test]
    fn SortRows() {
        let mut a = distinct_rows();
        let mut l = Multiset::<Q_PAD>::default();
        assert_eq!(a.SortRows(N, &mut l), 1);

        // the rows end up in descending histogram order
        let h = Multiset::<Q_PAD>::from_matrix::<N, N>(&a);
        for row in 1..N {
            assert!(h[row - 1] >= h[row], "row {} is out of order", row);
        }
    }

    /// sorting is a permutation of the rows, nothing is lost or duplicated
    #[test]
    fn SortRows_is_a_permutation() {
        let a = distinct_rows();
        let mut b = distinct_rows();
        let mut l = Multiset::<Q_PAD>::default();
        assert_eq!(b.SortRows(N, &mut l), 1);

        for row in 0..N {
            let mut found = false;
            for other in 0..N {
                if a[other] == b[row] { found = true; break; }
            }
            assert!(found, "row {} is not a row of the input", row);
        }
    }

    /// `L` is the early exit bound: if the shortest row of this matrix is
    /// already longer than the best one so far, give up
    #[test]
    fn SortRows_early_exit() {
        let mut a = distinct_rows();

        let mut l = Multiset::<Q_PAD>::default();
        l[0] = (N + 1) as u8;
        assert_eq!(a.SortRows(N, &mut l), 0);

        // a bound that cannot be beaten lets the sort run
        let mut l = Multiset::<Q_PAD>::default();
        assert_eq!(a.SortRows(N, &mut l), 1);
    }

    #[test]
    fn SortCols() {
        let mut a = distinct_rows();
        assert_eq!(a.SortCols(N), 1);

        // the columns end up in ascending order
        let mut t = MatrixNormalized::<N, N>::default();
        MatrixNormalized::transpose(&mut t, &a);
        for row in 1..N {
            assert!(t[row - 1] <= t[row], "column {} is out of order", row);
        }
    }

    /// the canonical form must be invariant under a monomial action: for any
    /// monomial Q, CF(G) == CF(G*Q)
    #[test]
    fn cf_is_invariant_under_a_monomial_action() {
        use crate::config::Xof;
        use crate::prng::csprng_initialize;

        let mut hits = 0;
        for seed in 0..8u8 {
            let mut s = csprng_initialize::<Xof>(&[seed; 32]);
            let g = MatrixNormalized::<N, N>::rand(&mut s);

            let mut a = MatrixNormalized::<N, N>::default();
            for row in 0..N { a[row].0 = g[row].0; }

            let monom = Monomial::<N, N>::rand(&mut s);
            let mut b = MatrixNormalized::<N, N>::default();
            MatrixNormalized::monomial_mul(&mut b, &g, &monom);

            if a.cf() == 1 && b.cf() == 1 {
                hits += 1;
                for row in 0..N {
                    assert_eq!(a[row], b[row], "seed {}: row {} differs", seed, row);
                }
            }
        }
        assert!(hits > 0, "cf never succeeded");
    }

    /// sorting twice must not move anything the second time
    #[test]
    fn SortCols_is_idempotent() {
        let mut a = distinct_rows();
        assert_eq!(a.SortCols(N), 1);

        let mut b = MatrixNormalized::<N, N>::default();
        for row in 0..N { b[row].0 = a[row].0; }

        assert_eq!(a.SortCols(N), 1);
        for row in 0..N {
            assert_eq!(a[row], b[row], "row {} moved on the second sort", row);
        }
    }
}