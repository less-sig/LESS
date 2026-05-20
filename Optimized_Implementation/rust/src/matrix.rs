// TODO remove
#![allow(dead_code)]

use std::{
    cmp,
    fmt,
    ops::{
        Add, Sub, Index, IndexMut
    },
};
use std::{ fmt::Display, fmt::Formatter, fmt::Result };
use sha3::{
    Digest, Sha3_256, Sha3_512, Shake128, Shake256, Shake128ReaderCore,
    digest::{ExtendableOutput, Update, XofReader},
    digest::core_api::XofReaderCoreWrapper,
};

use crate::fq::Fq;
use crate::monomial::Monomial;
use crate::vector::Vector;
use crate::prng::rand_range_q_state_elements;

const Q_PAD: usize = 128;
const Q: u8 = 127;

/// row major form
/// N: number of rows
/// M: number of cols
#[derive(Clone)]
pub struct Matrix<const N: usize, const M: usize> {
    rows: [Vector<M>; N],
}

/// N: number of rows
pub struct MatrixRREF<const N: usize> {
    rows: [Vector<N>; N],
    column_pos: [u16; N],
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

    /// same as init, zero initalized
    #[inline]
    pub fn new() -> Matrix<N, M> {
        Matrix {
            rows: core::array::from_fn(|_| Vector::<M>::init() ),
        }
    }

    /// same as init, but inits all values with `v`
    #[inline]
    pub fn from_u8(v: u8) -> Matrix<N, M> {
        Matrix {
            rows: core::array::from_fn(|_| Vector::<M>::from_u8(v) ),
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

    /// sample a random matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// TODO
    /// ```
    ///
    /// # Parameters
    /// - `seed`: NxM matrix
    pub fn rand() -> Self {
        let mut a = Self::init();
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                a[row][col] = Fq::rand();
            }
        }
        a
    }

    /// sample a random matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// TODO
    /// ```
    ///
    /// # Parameters
    /// - `seed`: NxM matrix
    pub fn rand_from_seed<S>(&mut self, seed: [u8; 32]) 
    where 
        S: ExtendableOutput + Default + Clone
    {
        let mut hasher = Shake128::default();
        hasher.update(&seed);
        let mut reader = hasher.finalize_xof();
        // TODO: is it possible to remove the `clone`?
        for i in 0..N {
            rand_range_q_state_elements(&mut reader, &mut self.rows[i].0);
        }
    }

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::Matrix;
    /// let result: Matrix<100, 100> = Matrix::init();
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
    /// let result: Matrix<100, 100> = Matrix::init();
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
    /// let result: Matrix<100, 100> = Matrix::init();
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
    /// let result: Matrix<100, 100> = Matrix::init();
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
    /// let mut C = Matrix::<N, N>::init();
    /// let A = Matrix::<N, N>::init();
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
    /// let mut a = Matrix::<32, 32>::new();
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
    #[must_use]
    pub fn rref<const CT: bool>(&mut self,
                                is_pivot_column: &mut [u8; M],
                                was_pivot_column: &mut [u8; M],
                                pvt_reuse_limit: u32) -> bool {
        let mut pvt_reuse_cnt: u32 = 0;

        if pvt_reuse_limit != 0 {
            // TODO: loop boundaries propably not correct
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
    /// ```
    ///
    /// # Parameters
    /// - `r`: NxM matrix
    /// - `r1`: index of the first row
    /// - `r2`: index of the second row
    pub fn compress(
        &self,
        compressed: &mut [u8],
        is_pivot_column: &[u8; N],
    ) {
        /* Compress pivot flags */
        for col_byte in 0..(N / 8) {
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

    /// Expands a compressed RREF generator matrix into a full one
    /// \param full[out]: output generator matrix (K \times N) 
    /// \param compact[out]: input compressed generator matrix (K \times N-K) 
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
    pub fn expand_to_rref(&mut self,
                          compact: MatrixRREF<N>) {
        let mut placed_dense_cols: usize = 0;
        for col_idx in 0..M {
            if placed_dense_cols < N-M && col_idx as u16 == compact.column_pos[placed_dense_cols] {
                for row_idx in 0..N {
                    self[row_idx][col_idx] = compact.rows[row_idx][placed_dense_cols];
                }
                placed_dense_cols += 1;
            } else {
                for row_idx in 0..N {
                    self[row_idx][col_idx] = Fq((row_idx == (col_idx - placed_dense_cols)) as u8)
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


#[derive(Debug)]
pub struct MatrixNormalized<const N: usize> {
    rows: [Vector<N>; N],
}

impl<const N: usize> MatrixNormalized<N> {
    /// generates a new matrix init with 0
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100> = MatrixNormalized::init();
    /// assert_eq!(result[0][0].0, 0);
    /// assert_eq!(result.dimension(), 100);
    /// ```
    ///
    #[inline]
    pub fn init() -> Self {
        Self {
            rows: core::array::from_fn(|_| Vector::<N>::init() ),
        }
    }

    /// same as init
    pub fn new() -> MatrixNormalized<N> {
        MatrixNormalized {
            rows: core::array::from_fn(|_| Vector::<N>::init() ),
        }
    }
   
    /// generates a new matrix init with q-1
    pub fn new_large() -> MatrixNormalized<N> {
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
        let mut a = Self::init();
        for i in 0..N {
            a[i][i] = Fq(1);
        }
        a
    }

    /// The size/dimension of the matrix
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<100> = MatrixNormalized::init();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::init();
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
    /// let result: MatrixNormalized<100> = MatrixNormalized::init();
    /// assert_eq!(result.nrows(), 100);
    /// ```
    ///
    /// # Returns
    /// n_rows
    #[inline]
    pub fn nrows(&self) -> usize {
        N
    }

    /// TODO doc, opt
    /// static internal function implementing the histogram function
    #[inline]
    fn sort(out: &mut [u8; Q_PAD], row: &Vector<N>) {
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

    /// row addition
    /// rows[i] += rows[j]
    /// # Examples
    ///
    /// ```
    /// use less::matrix::MatrixNormalized;
    /// let result: MatrixNormalized<64> = MatrixNormalized::init();
    /// ```
    ///
    fn add(&mut self, i: usize, j: usize) {
        // what we want tonimplement
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
    /// let result: MatrixNormalized<64> = MatrixNormalized::init();
    /// ```
    ///
    fn sub(&mut self, i: usize, j: usize) {
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
    /// let result: MatrixNormalized<64> = MatrixNormalized::init();
    /// ```
    ///
    fn scalar(&mut self, i: usize, s: Fq) {
        Vector::<N>::scalar2(&mut self.rows[i], s);
    }

    /// \return sum(row)
    fn acc(row: &Vector<N>) -> Fq {
        Vector::<N>::acc(row)
    }

    /// \return sum(row**{-1})
    fn acc_inv(row: &Vector<N>) -> Fq {
        Vector::<N>::acc_inv(row)
    }

    /// insert sort 
    fn sort_rows(&mut self) -> i32 {
        let mut tmp = [[0u8; Q_PAD]; N];

        for i in 0..N {
            Self::sort(&mut tmp[i], &self.rows[i]);
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
            let mut s = Self::acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.scalar(i as usize, s);
        };

        return self.sort_cf();
    }

    /// by floyd
    fn scale_cf_preprocess_v1(&mut self, z: usize, m: & [u8; Q_PAD]) -> i32 {
        let mut tmp = [0u8; Q_PAD];
        for i in 0..z  {
            let mut s = Self::acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.scalar(i as usize, s);
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
            let mut s = Self::acc(&self.rows[i]);
            if s.0 != 0 {
                s = Fq::inv(s);
            } else {
                s = Self::acc_inv(&self.rows[i]);
                if s.0 == 0 {
                    continue;
                }
            }

            self.scalar(i as usize, s);
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

    //pub fn cf_opt1(&mut self) {
    //    let mut M = Self::new_large();
    //    let mut J = [0; N];
    //    let mut z = 0;
    //    let mut num_zeros = 0;
    //}

    ///
    pub fn blind<S>(&mut self, state: S) {
        // TODO
    }

    ///
    pub fn CanonicalForm(&mut self) {
        // TODO
    }
}

//impl<const N: usize> Add for MatrixNormalized<N> {
//    type Output = MatrixNormalized<N>;
//    fn add(self, r: Self) -> Self::Output {
//        let mut t = self;
//        Self::add(&mut t, &r);
//        t
//    }
//}
//impl<'a, const N: usize> Add for &'a mut MatrixNormalized<N> {
//    type Output = &'a mut MatrixNormalized<N>;
//    fn add(self, r: Self) -> Self::Output {
//        Matrix::<N, M>::add(self, r);
//        self
//    }
//}
//impl<const N: usize> Sub for MatrixNormalized<N> {
//    type Output = MatrixNormalized<N>;
//    fn sub(self, r: Self) -> Self::Output {
//        let mut t = self;
//        Self::sub(&mut t, &r);
//        t
//    }
//}
//impl<'a, const N: usize> Sub for &'a mut MatrixNormalized<N> {
//    type Output = &'a mut MatrixNormalized<N>;
//    fn sub(self, r: Self) -> Self::Output {
//        MatrixNormalized::<N>::sub(self, r);
//        self
//    }
//}
//

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
//impl<const N: usize> fmt::Debug for MatrixNormalized<N> {
//    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//        writeln!(f, "{:?}", self.rows)
//    }
//}


/// use less::matrix::Matrix;
#[cfg(test)]
mod tests {
    use super::*;
    const N: usize = 32;
    const M: usize = 32;
    
    #[test]
    fn init() {
        let _ = Matrix::<N, M>::init();
    }
    #[test]
    fn new() {
        let _ = Matrix::<N, M>::new();
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
        let a = Matrix::<N, M>::rand();
        println!("{}", a);
        for row in 0..a.nrows() {
            for col in 0..a.ncols() {
                assert!(a[row][col].0 <= Fq::Q);
            }
        }

        assert!(a.is_zero() == false);
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
        let a = Matrix::<N, M>::new();
        assert_eq!(a.ncols(), M);
    }
    #[test]
    fn nrows() {
        let a = Matrix::<N, M>::new();
        assert_eq!(a.nrows(), N);
    }

    #[test]
    fn add() {
        let a = Matrix::<N, M>::new();
        let b = Matrix::<N, M>::from_u8(1);
        let c = a.clone().add(b.clone());
        for i in 0..N {
            assert!(c[i][i].0 == 1);
        }
        let c = a + b;
        for i in 0..N {
            assert!(c[i][i].0 == 1);
        }
    }

    #[test]
    fn sub() {
        let a = Matrix::<N, M>::new();
        let b = Matrix::<N, M>::from_u8(1);
        let c = b - a;
        for i in 0..N {
            assert!(c[i][i].0 == 1);
        }
    }

    #[test]
    fn monomial_mul() {
        let mut a = Matrix::<N, M>::new();
        let t = Monomial::<N>::new();
        // TODO
    }

    #[test]
    fn swap_rows() {
        let mut a = Matrix::<N, M>::new();
        let i1 = 0; let i2 = N/2;
        for i in 0..a.ncols() {
            a[i1][i] = Fq(1);
        }
    
        a.swap_rows(i1, i2);
        for i in 0..N {
            assert!(a[i1][i].0 == 0);
            assert!(a[i2][i].0 == 1);
        }
    }

    #[test]
    fn rref() {
        let mut a = Matrix::<N, M>::rand();
        let mut is_pivot_column = [0u8; M];
        let mut was_pivot_column = [0u8; M];
        let pvt_reuse_limit = M as u32;
        let b = a.rref::<false>(&mut is_pivot_column, &mut was_pivot_column, pvt_reuse_limit);
        assert!(b);
    }
}
