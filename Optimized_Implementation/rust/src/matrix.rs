// TODO remove
#![allow(dead_code)]
#![allow(unused_imports)]

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
pub struct Matrix<const N: usize, const M: usize> {
    rows: [Vector<M>; N],
}


impl<const N: usize, const M: usize> Matrix<N, M> {

    /// Zero Initialised
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::default();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), (100, 100));
    /// ```
    ///
    /// # Returns
    /// A new Matrix
    #[inline]
    pub fn default() -> Self {
        Self {
            rows: core::array::from_fn(|_| Vector::<M>::default() ),
        }
    }

    /// same as init, but inits all values with `v`
    #[inline]
    pub fn from_u8(v: u8) -> Matrix<N, M> {
        Matrix {
            rows: core::array::from_fn(|_| Vector::<M>::from_u8(v) ),
        }
    }

    /// translation of `generator_rref_expand`
    /// Expands a compressed RREF generator matrix into a full one
    /// \param full[out]: output generator matrix (K \times N)
    /// \param compact[out]: input compressed generator matrix (K \times N-K)
    pub fn from_compressed_rref<const M_PRIME: usize>(column_pos: &[u16; M_PRIME],compact: &Matrix<N, M_PRIME>) -> Self {
        // static_assert!(M <= M_PRIME);
        let mut placed_dense_cols = 0usize;
        let mut full = Self::default();
        for col_idx in 0..M {
            if (placed_dense_cols < M - N) && (col_idx as u16 == column_pos[placed_dense_cols]) {
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
    /// let result: Matrix<100, 100> = Matrix::id();
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
            rand_range_q_state_elements(&mut a.rows[i].0, state);
        }
        a
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
    /// let A: Matrix<100, 100> = Matrix::default();
    /// let B: Matrix<100, 100> = Matrix::default();
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
    /// let A: Matrix<100, 100> = Matrix::default();
    /// let B: Matrix<100, 100> = Matrix::default();
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
    /// let result: Matrix<100, 100> = Matrix::default();
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
    /// let result: Matrix<100, 100> = Matrix::default();
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
    /// let result: Matrix<100, 100> = Matrix::default();
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
    /// let result: Matrix<100, 100> = Matrix::default();
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
            Vector::<M>::add2(&mut r[i],&l[i]);
        }
    }

    /// constant time subtraction
    /// r -= l
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// const N: usize = 128;
    /// let mut C = Matrix::<N, N>::default();
    /// let A = Matrix::<N, N>::default();
    /// Matrix::sub(&mut C, &A)
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
    /// let mut a = Matrix::<32, 32>::default();
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
    /// use less::matrix::Matrix;
    /// TODO
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn monomial_mul<const M_prime: usize>(res: &mut Self, g: &Self, monom: &Monomial<M_prime>) {
        // static_assert!(M <= M_prime);
        for src_col_idx in 0..M {
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
    /// use less::matrix::Matrix;
    /// TODO
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

        if pvt_reuse_limit != 0 {
            // TODO: loop boundaries probably not correct
            for preproc_col in N - 1 .. 0 {
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

            for col in pivot_column..N {
                self[pivot_row][col] = scaling_factor * self[pivot_row][col];
            }

            /* Eliminate pivot column from other rows */
            for row_idx in 0..N {
                if row_idx != pivot_row {
                    let multiplier = self[row_idx][pivot_column];
                    for col_idx in 0..N {
                        let tmp = multiplier * self[pivot_row][col_idx];
                        self[row_idx][col_idx] -= tmp;
                    }
                }
            }
        }

        true
    }

    /// Compresses a generator matrix in RREF into a byte array
    ///
    /// * `compressed` – output buffer (RREF_MAT_PACKEDBYTES bytes)
    /// * `is_pivot_column` – length-N array, exactly K entries set to 1
    /// constant time row swap
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// TODO
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
        for col_byte in 0..(M / 8) {
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
        if M == 252 || M == 548 {
            compressed[N / 8] =
                (is_pivot_column[N - 4] << 0) |
                (is_pivot_column[N - 3] << 1) |
                (is_pivot_column[N - 2] << 2) |
                (is_pivot_column[N - 1] << 3);

            compress_idx = N / 8 + 1;
        } else {
            compress_idx = N / 8;
        }

        /* Compress non-pivot columns row-by-row */
        let mut encode_state: u8 = 0;

        for row_idx in 0..N {
            for col_idx in 0..M {
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
    /// use less::matrix::Matrix;
    /// TODO
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn expand(
        &mut self,
        compressed: &[u8],
        is_pivot_column: &mut [u8; N]) {
        /* Initialize pivot flags */
        for i in 0..N {
            is_pivot_column[i] = 0;
        }

        /* Decompress pivot flags */
        for col_byte in 0..(N / 8) {
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
        if M == 252 || M == 548 {
            let b = compressed[N / 8];
            is_pivot_column[N - 4] = (b >> 0) & 0x1;
            is_pivot_column[N - 3] = (b >> 1) & 0x1;
            is_pivot_column[N - 2] = (b >> 2) & 0x1;
            is_pivot_column[N - 1] = (b >> 3) & 0x1;

            compress_idx = N / 8 + 1;
        } else {
            compress_idx = N / 8;
        }

        /* Decompress matrix row-by-row */
        let mut decode_state: u8 = 0;

        for row_idx in 0..N {
            let mut pivot_idx: usize = 0;

            for col_idx in 0..M {
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
impl<const N: usize, const M: usize> Sub for Matrix<N, M> {
    type Output = Matrix<N, M>;
    fn sub(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::sub(&mut t, &r);
        t
    }
}
impl<'a, const N: usize, const M: usize> Sub for &'a mut Matrix<N, M> {
    type Output = &'a mut Matrix<N, M>;
    fn sub(self, r: Self) -> Self::Output {
        Matrix::<N, M>::sub(self, r);
        self
    }
}
impl<const N: usize, const M: usize> Mul<Monomial<M>> for Matrix<N, M> {
    type Output = Matrix<N, M>;
    fn mul(self, r: Monomial<M>) -> Self::Output {
        let mut t = Matrix::default();
        Self::monomial_mul(&mut t, &self, &r);
        t
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

impl<const N: usize, const M: usize> Display for Matrix<N, M> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for i in 0..N {
            write!(f, "{}", self[i])?;
        }

        Ok(())
    }
}


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct MatrixNormalized<const N: usize> {
    pub rows: [Vector<N>; N],
}

impl<const N: usize> MatrixNormalized<N> {
    /// generates a new matrix init with 0
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100> = MatrixNormalized::default();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    pub fn default() -> MatrixNormalized<N> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N>::default() ),
        }
    }
   
    /// generates a new matrix init with q-1
    pub fn default_large() -> MatrixNormalized<N> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N>::from_u8(Q-1) ),
        }
    }

    /// same as init, but inits all values with `v`
    #[inline]
    pub fn from_u8(v: u8) -> MatrixNormalized<N> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N>::from_u8(v) ),
        }
    }

    /// translation of 'void normalized_copy_from_generator_non_information_set'
    /// \param A[out]: pointer to allocated normalized struct, which get filled with the
    ///     non-IS of the generator matrix G
    /// \param G[in]: generator matrix to extract the non-IS from.
    /// \param is_pivot_column[in]: array identifying a pivot column via a 1
    #[inline]
    pub fn from_information_set<const N_PRIME: usize, const M_PRIME: usize>(g: &Matrix<N_PRIME, M_PRIME>, is_pivot_column: &[u8; M_PRIME]) -> MatrixNormalized<N> {
        let mut ret = MatrixNormalized::default();
        // we simply copy the last N-K columns even if they are not the information set.
        for i in 0..N {
            for j in 0..N {
                ret[i][j] = g[i][j];
            }
        }
        // now we scan if we need to fix the non information set
        let mut ctr = 0;
        for _ in 0.. N {
            if is_pivot_column[ctr] == 1 {
                ctr += 1;
            }
        }

        // easy part: the last N-K columns are the non IS
        if ctr == N { return ret; }

        // "hard" part: copy all remaining columns < K into the non information set part
        for i in 0..N {
            if is_pivot_column[ctr] == 0 { break; }

            // copy column
            for j in 0..N {
                ret[j][i] = g[j][ctr];
            }
            ctr += 1
        }

        ret
    }

    /// NOTE: not constant time
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let A: MatrixNormalized<100> = MatrixNormalized::id();
    /// let B: MatrixNormalized<100> = MatrixNormalized::id();
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
    /// let A: MatrixNormalized<100> = MatrixNormalized::id();
    /// let B: MatrixNormalized<100> = MatrixNormalized::id();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::id();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::rand_from_seed::<Shake128>(&[0u8; 32]);
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::default();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::default();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::default();
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
    /// let mut C: MatrixNormalized<N> = MatrixNormalized::default();
    /// let A : MatrixNormalized<N> = MatrixNormalized::default();
    /// MatrixNormalized::add(&mut C, &A);
    /// ```
    ///
    pub fn add(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<N>::add2(&mut r[i], &l[i]);
        }
    }

    /// matrix subtraction
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// const N: usize = 128;
    /// let mut C: MatrixNormalized<N> = MatrixNormalized::default();
    /// let A : MatrixNormalized<N> = MatrixNormalized::default();
    /// MatrixNormalized::add(&mut C, &A);
    /// ```
    ///
    pub fn sub(r: &mut Self, l: &Self) {
        for i in 0..N {
            Vector::<N>::sub2(&mut r[i], &l[i]);
        }
    }

    /// matrix scalar multiplication
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// use less::fq::Fq;
    /// const N: usize = 128;
    /// let mut C: MatrixNormalized<N> = MatrixNormalized::default();
    /// MatrixNormalized::scalar(&mut C, Fq(1));
    /// ```
    ///
    pub fn scalar(r: &mut Self, l: Fq) {
        for i in 0..N {
            Vector::<N>::scalar2(&mut r[i], l);
        }
    }

    /// row addition
    /// rows[i] += rows[j]
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut result: MatrixNormalized<64> = MatrixNormalized::default();
    /// MatrixNormalized::row_add(&mut result, 0, 1);
    /// ```
    ///
    pub fn row_add(&mut self, i: usize, j: usize) {
        // what we want to implement
        // Vector::<N>::add2(&mut self.rows[i], &mut self.rows[j]);
        if i == j {
            return;
        }

        let [a, b] = self.rows.get_disjoint_mut([i, j]).unwrap();
        Vector::<N>::add2(a, b);
    }

    /// row subtraction
    /// rows[i] -= rows[j]
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut result: MatrixNormalized<64> = MatrixNormalized::default();
    /// MatrixNormalized::row_sub(&mut result, 0, 1);
    /// ```
    ///
    pub fn row_sub(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        let [a, b] = self.rows.get_disjoint_mut([i, j]).unwrap();
        Vector::<N>::sub2(a, b);
    }

    /// row scalar multiplication
    /// rows[i] *= s 
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// use less::fq::Fq;
    /// let mut result: MatrixNormalized<64> = MatrixNormalized::default();
    /// MatrixNormalized::row_scalar(&mut result, 0, Fq(1));
    /// ```
    ///
    pub fn row_scalar(&mut self, i: usize, s: Fq) {
        Vector::<N>::scalar2(&mut self.rows[i], s);
    }

    /// transpose: a = b^T
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let mut A: MatrixNormalized<64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64> = MatrixNormalized::default();
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
        // TODO L is ignored,
        let mut tmp: [Multiset<Q_PAD>; N] = Multiset::<Q_PAD>::from_matrix::<N>(self);
        for i in 1..z {
            let mut j = i;
            while (j > 0) && (tmp[j-1] < tmp[j]) {
                if let Ok([a, b]) = &mut tmp.get_disjoint_mut([j-1, j]) {
                    std::mem::swap(a, b);
                    if let Ok([a, b]) = &mut self.rows.get_disjoint_mut([j-1, j]) {
                        std::mem::swap(a, b);
                    }
                }
                j -= 1;
             }
        }
        0
    }

    /// TODO example and test
    fn SortCols(&mut self, z: usize) -> i32 {
        let mut VT = MatrixNormalized::default();
        MatrixNormalized::transpose(&mut VT, &self);

        for i in 1..z {
            let mut j = i;
            while (j > 0) && (VT[j-1] < VT[j]) {
                if let Ok([a, b]) = &mut VT.rows.get_disjoint_mut([j-1, j]) {
                    std::mem::swap(a, b);
                }
                j -= 1;
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
            let mut s = Vector::<N>::acc(&self.rows[i]);
            if s != Fq(0) {
                s = Fq::inv(s);
            } else {
                s = Vector::<N>::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_scalar(i as usize, s);
            Multiset::<Q_PAD>::from_row(&mut tmp, &self.rows[i]);
            if tmp < *m {
                return 0;
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
            let mut s = Vector::<N>::acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Vector::<N>::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.row_scalar(i as usize, s);
            Multiset::<Q_PAD>::from_row(&mut tmp, &self.rows[i]);
            if tmp < *m  {
                ret = 1;
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
        self.SortCols(N); // NOTE: must be N_PAD
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
                s = Vector::<N>::acc_inv(&self.rows[row]);
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
    /// let mut A: MatrixNormalized<64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64> = MatrixNormalized::default();
    /// // TODO
    /// ```
    pub fn monomial_mul(res: &mut Self, g: &Self, monom: &Monomial<N>) {
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
    /// let mut A: MatrixNormalized<64> = MatrixNormalized::default();
    /// let mut B: MatrixNormalized<64> = MatrixNormalized::default();
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
        let left = Monomial::<N>::rand(state);
        let right = Monomial::<N>::rand(state);

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
        let mut M = MatrixNormalized::<N>::default_large();
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

        if z == N {
            assert!(false);
            // TODO return self.CFOriginal();
        }

        let mut row_inv_data = Vector::<N>::default();
        // NOTE: this is already "sorted"
        let mut L = Multiset::<Q_PAD>::default();

        // Check smallest rows of all matrices to find the smallest candidate
        for row in 0..N {
            if Z[row] { continue; }

            let mut B = MatrixNormalized::<N>::default();
            Vector::inv_non_ct(&mut row_inv_data, &self[row]);
            for row2 in 0..z {
                Vector::mul(&mut B[row2], &self[J[row2] as usize], &row_inv_data);
            }

            let ret = B.ScaleCF(&mut L);
            if ret == 1 && B < M {
                // TODO sort(B[0])
                touched = 1;
                M = B;
            }
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

        if z == N {
            // TODO
            assert!(false);
            return 0;
        }

        let mut B = MatrixNormalized::<N>::default();
        let mut row_inv_data = Vector::<N>::default();
        // NOTE: this is already "sorted"
        let mut L = Multiset::<Q_PAD>::default();

        // Check smallest rows of all matricies to find the smallest candidate
        for row in 0..N {
            if Z[row] { continue; }

            Vector::inv_non_ct(&mut row_inv_data, &self[row]);
            for row2 in 0..z {
                Vector::mul(&mut B[row2], &self[J[row2] as usize], &row_inv_data);
            }

            if B.ScaleCFSubPreprocess(z, &mut L) == 0 {
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

        // TODO self = B;
        ret
    }
}

impl<const N: usize> Add for MatrixNormalized<N> {
    type Output = MatrixNormalized<N>;
    fn add(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::add(&mut t, &r);
        t
    }
}
impl<'a, const N: usize> Add for &'a mut MatrixNormalized<N> {
    type Output = &'a mut MatrixNormalized<N>;
    fn add(self, r: Self) -> Self::Output {
        MatrixNormalized::<N>::add(self, r);
        self
    }
}
impl<const N: usize> Sub for MatrixNormalized<N> {
    type Output = MatrixNormalized<N>;
    fn sub(self, r: Self) -> Self::Output {
        let mut t = self;
        Self::sub(&mut t, &r);
        t
    }
}
impl<'a, const N: usize> Sub for &'a mut MatrixNormalized<N> {
    type Output = &'a mut MatrixNormalized<N>;
    fn sub(self, r: Self) -> Self::Output {
        MatrixNormalized::<N>::sub(self, r);
        self
    }
}

impl<const N: usize> Mul<Permutation<N>> for MatrixNormalized<N> {
    type Output = MatrixNormalized<N>;
    fn mul(self, r: Permutation<N>) -> Self::Output {
        let mut t = MatrixNormalized::default();
        Self::permutation_mul(&mut t, &self, &r);
        t
    }
}

/// Direct access to a row
impl<const N: usize> Index<usize> for MatrixNormalized<N> {
    type Output = Vector<N>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}
/// Direct access to a row (mutable)
impl<const N: usize> IndexMut<usize> for MatrixNormalized<N> {
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
        let _ = Matrix::<N, M>::default();
    }
    #[test]
    fn from_u8() {
        let a = Matrix::<N, M>::from_u8(1);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert_eq!(a[row][col].0, 1);
            }
        }
    }

    #[test]
    fn rand() {
        let a = Matrix::<N, M>::rand_from_seed::<Shake128>(&[0u8; 32]);
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
        let a = Matrix::<N, M>::rand_from_seed::<Shake128>(&seed);
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
        let a = Matrix::<N, M>::id();
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                let b = (row == col) as u8;
                assert_eq!(a[row][col].0, b);
            }
        }
    }
    #[test]
    fn ncols() {
        let a = Matrix::<N, M>::default();
        assert_eq!(a.ncols(), M);
    }
    #[test]
    fn nrows() {
        let a = Matrix::<N, M>::default();
        assert_eq!(a.nrows(), N);
    }

    #[test]
    fn add() {
        let a = Matrix::<N, M>::default();
        let b = Matrix::<N, M>::from_u8(1);
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
        let a = Matrix::<N, M>::default();
        let b = Matrix::<N, M>::from_u8(1);
        let c = b - a;
        for i in 0..N {
            assert_eq!(c[i][i].0, 1);
        }
    }

    #[test]
    fn monomial_mul() {
        let mut a = Matrix::<N, M>::default();
        let t = Monomial::<N>::default();
        // TODO
    }

    #[test]
    fn swap_rows() {
        let mut a = Matrix::<N, M>::default();
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
        let mut a = Matrix::<N, M>::rand_from_seed::<Shake128>(&[0u8; 32]);
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
        let _ = MatrixNormalized::<N>::default();
    }

    #[test]
    fn from_u8() {
        let a = MatrixNormalized::<N>::from_u8(1);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert_eq!(a[row][col], Fq(1));
            }
        }
    }

    #[test]
    fn transpose() {
        let a = MatrixNormalized::<N>::rand_from_seed::<Shake128>(&[0u8; 32]);
        let mut b = MatrixNormalized::<N>::default();
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

        let mut a = MatrixNormalized::<N>::default();
        MatrixNormalized::blind(&mut a, &mut reader);

    }
}