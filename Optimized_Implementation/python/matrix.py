#!/usr/bin/env python3
""" super simple matrix implementation. The only goal is to have zero dependencies """
import random
from functools import cmp_to_key

from config import K, N, N8, Q, RREF_MAT_PACKEDBYTES
from fq import Fq
from monomial import Monomial
from multiset import Multiset
from prng import csprng_initialize
from vector import Vector


def monomial_mul(rows: list[Vector], monom: Monomial) -> list[Vector]:
    """ res = g * monom """
    ret = [Vector.default(len(row)) for row in rows]
    for src_col_idx, (pos, coeff) in enumerate(zip(monom.perms, monom.coeffs)):
        for out, row in zip(ret, rows):
            out[pos] = row[src_col_idx] * coeff
    return ret


class Matrix:
    """ simple matrix class
    NOTE: row major
    NOTE: always zero initialized
    """
    def __init__(self, nrows: int, ncols: int, q: int = Q) -> None:
        """ zero initialized """
        self.nrows = nrows
        self.ncols = ncols
        self.q = q
        self.rows = [Vector.default(ncols) for _ in range(nrows)]

    @classmethod
    def from_rows(cls, rows: list[Vector]):
        ret = cls(len(rows), len(rows[0]))
        ret.rows = rows
        return ret

    def __getitem__(self, tup):
        """ nice access function
        NOTE: access it via:
            A[i, j] or  (returns field element)
            A[i]        (returns row)
        """
        if isinstance(tup, tuple):
            return self.rows[tup[0]][tup[1]]
        return self.rows[tup]

    def __setitem__(self, tup, data):
        """ nice access function
        NOTE: access it via:
            A[i, j] or  (returns field element)
            A[i]        (returns row)
        """
        if isinstance(tup, tuple):
            self.rows[tup[0]][tup[1]] = Fq(data)
        else:
            self.rows[tup] = data

    def __mul__(self, B):
        if isinstance(B, Matrix):
            return self.mul(B)
        return self.from_rows(monomial_mul(self.rows, B))

    def __eq__(self, B) -> bool:
        return self.rows == B.rows

    def __str__(self) -> str:
        return "".join(" ".join(str(v).rjust(3) for v in row) + "\n" for row in self.rows)

    def print(self):
        print(self)

    def id(self) -> 'Matrix':
        """ identity matrix """
        self.rows = [Vector(Fq(i == j) for j in range(self.ncols)) for i in range(self.nrows)]
        return self

    def zero(self) -> 'Matrix':
        """ zeros all elements"""
        self.rows = [Vector.default(self.ncols) for _ in range(self.nrows)]
        return self

    def set(self, d: list[list[int]]) -> 'Matrix':
        """ sets every element within in the matrix to `d`"""
        self.rows = [Vector(Fq(v) for v in row) for row in d]
        return self

    def random(self) -> 'Matrix':
        """ generates a random matrix """
        self.rows = [Vector(Fq(random.randrange(self.q)) for _ in range(self.ncols)) for _ in range(self.nrows)]
        return self

    def random_row_with_weight(self, row: int, w: int) -> 'Matrix':
        """ generates a random weight w row """
        assert 0 < w < self.ncols
        self.zero()
        self.rows[row][:w] = [Fq(random.randrange(1, self.q)) for _ in range(w)]
        random.shuffle(self.rows[row])
        return self

    def transpose(self) -> 'Matrix':
        """ simple transpose, not inplace """
        return Matrix.from_rows([Vector(c) for c in zip(*self.rows)])

    def mul(self, B: 'Matrix') -> 'Matrix':
        """ simple multiplication """
        assert self.ncols == B.nrows
        cols = [Vector(c) for c in zip(*B.rows)]
        return Matrix.from_rows([Vector(row.mul(col).acc() for col in cols) for row in self.rows])

    def scale_row(self, row: int, v: Fq) -> 'Matrix':
        """ scales the row `row` by v"""
        self.rows[row] = self.rows[row].scalar_mul(v)
        return self

    @staticmethod
    def rand(state, nrows: int, ncols: int) -> 'Matrix':
        """ sample a random matrix """
        return Matrix.from_rows([Vector.rand(state, ncols) for _ in range(nrows)])

    @staticmethod
    def rand_from_seed(seed: bytes, nrows: int, ncols: int) -> 'Matrix':
        return Matrix.rand(csprng_initialize(seed), nrows, ncols)

    @staticmethod
    def from_compressed_rref(column_pos: list[int], compact: 'Matrix') -> 'Matrix':
        """
        translation of `generator_rref_expand`
        Expands a compressed RREF generator matrix into a full one
        :param compact[in]: input compressed generator matrix (K \times N-K)
        :return full: output generator matrix (K \times N)
        """
        full = Matrix(compact.nrows, N)
        placed_dense_cols = 0
        for col_idx in range(N):
            if placed_dense_cols < N - K and col_idx == column_pos[placed_dense_cols]:
                # non-pivot column, restore one full column
                for row_idx in range(full.nrows):
                    full.rows[row_idx][col_idx] = compact.rows[row_idx][placed_dense_cols]
                placed_dense_cols += 1
            else:
                # regenerate the appropriate pivot column
                for row_idx in range(full.nrows):
                    full.rows[row_idx][col_idx] = Fq(row_idx == col_idx - placed_dense_cols)
        return full

    def rref(self,
             is_pivot_column: list[int],
             was_pivot_column: list[int],
             pvt_reuse_limit: int) -> bool:
        """ TODO doc string """
        rows, n, m = self.rows, self.nrows, self.ncols
        pvt_reuse_cnt = 0

        # row swap pre-process - swap previous pivot elements to corresponding
        # row to reduce likelihood of corruption
        if pvt_reuse_limit != 0:
            for preproc_col in reversed(range(n)):
                if was_pivot_column[preproc_col] == 1:
                    pivot_el_row = 0
                    for row in range(n):
                        if rows[row][preproc_col] != 0:
                            pivot_el_row = row
                    rows[preproc_col], rows[pivot_el_row] = rows[pivot_el_row], rows[preproc_col]

        for row_to_reduce in range(n):
            pivot_row = pivot_column = row_to_reduce

            # Search pivot
            while pivot_column < m and rows[pivot_row][pivot_column] == 0:
                while pivot_row < n and rows[pivot_row][pivot_column] == 0:
                    pivot_row += 1

                if pivot_row >= n:
                    pivot_column += 1
                    pivot_row = row_to_reduce

            if pivot_column >= m:
                return False

            is_pivot_column[pivot_column] = 1

            # Swap rows if needed
            if row_to_reduce != pivot_row:
                # pivot no longer reusable
                was_pivot_column[pivot_row] = 0
                rows[row_to_reduce], rows[pivot_row] = rows[pivot_row], rows[row_to_reduce]

            pivot_row = row_to_reduce

            # Pivot reuse shortcut
            if was_pivot_column[pivot_column] == 1 and pvt_reuse_cnt < pvt_reuse_limit and pivot_column < n:
                pvt_reuse_cnt += 1
                continue

            # Rescale pivot row
            scaling_factor = rows[pivot_row][pivot_column].inv()
            p = rows[pivot_row]
            p[pivot_column:] = [scaling_factor * v for v in p[pivot_column:]]

            # Eliminate pivot column from other rows
            for row_idx in range(n):
                multiplier = rows[row_idx][pivot_column]
                if row_idx != pivot_row and multiplier != 0:
                    rows[row_idx] = Vector(a - multiplier * b for a, b in zip(rows[row_idx], p))

        return True

    @staticmethod
    def unpack_coset_rep(g: 'Matrix',
                         c: bytes,
                         initial_g_col_pivot: list[int],
                         permuted_g_col_pivot: list[int]) -> 'Matrix':
        """
        translation of `UnpackCosetRep`
        :param g[in]: current generator matrix: K \times N-K
        :param c[in]: compressed cf action
        :param initial_g_col_pivot[in]: input IS
        :param permuted_g_col_pivot[out]: output IS, to keep track of the pivot cols
        :return res: G*c: a generator matrix: K \times N-K
        """
        res = Matrix(g.nrows, g.ncols)
        l, r = 0, K
        for col in range(N):
            if (c[col // 8] >> (col % 8)) & 1:
                pos, l = l, l + 1
            else:
                pos, r = r, r + 1

            permuted_g_col_pivot[pos] = initial_g_col_pivot[col]

            # copy the column
            for k in range(g.nrows):
                res.rows[k][pos] = g.rows[k][col]
        return res

    def compress(self, is_pivot_column: list[int]) -> bytes:
        """ Compresses a generator matrix in RREF into a byte array """
        # Compress pivot flags
        flags = sum(b << col for col, b in enumerate(is_pivot_column[:N]))

        # Compress non-pivot columns row-by-row
        values = [v.get() for row in self.rows for col, v in enumerate(row) if is_pivot_column[col] == 0]
        packed = sum(v << (7 * i) for i, v in enumerate(values))
        return flags.to_bytes(N8, "little") + packed.to_bytes(RREF_MAT_PACKEDBYTES - N8, "little")

    @staticmethod
    def expand(compressed: bytes, is_pivot_column: list[int]) -> 'Matrix':
        """
        Expands a compressed RREF generator matrix into a full one
        :param is_pivot_column: length-N array, will be initialized here
        """
        # Decompress pivot flags
        flags = int.from_bytes(compressed[:N8], "little")
        for col in range(N):
            is_pivot_column[col] = (flags >> col) & 1

        # Decompress matrix row-by-row
        packed = int.from_bytes(compressed[N8:RREF_MAT_PACKEDBYTES], "little")
        ret = Matrix(K, N)
        i = 0
        for row_idx, row in enumerate(ret.rows):
            pivot_idx = 0
            for col in range(N):
                if is_pivot_column[col] == 0:
                    row[col] = Fq((packed >> (7 * i)) & Q)
                    i += 1
                else:
                    # Pivot column
                    row[col] = Fq(row_idx == pivot_idx)
                    pivot_idx += 1
        return ret


class MatrixNormalized(Matrix):
    @classmethod
    def from_information_set(cls, g: Matrix, is_pivot_column: list[int]) -> 'MatrixNormalized':
        """
        translation of 'void normalized_copy_from_generator_non_information_set'
        :param g[in]: generator matrix to extract the non-IS from.
        :param is_pivot_column[in]: array identifying a pivot column via a 1
        """
        # we simply copy the last N-K columns even if they are not the information set.
        ret = cls.from_rows([Vector(row[K:N]) for row in g.rows[:K]])
        # now we scan if we need to fix the non information set
        ctr = 0
        while ctr < K and is_pivot_column[ctr] == 1:
            ctr += 1

        # easy part: the last N-K columns are the non IS
        if ctr == K:
            return ret

        # "hard" part: copy all remaining columns < K into the non information set part
        j = 0
        while j < N - K and is_pivot_column[ctr] == 0:
            # copy column
            for k in range(K):
                ret.rows[k][j] = g.rows[k][ctr]
            ctr += 1
            j += 1
        return ret

    def blind(self, state) -> None:
        # the randomness for `right` is drawn before the one for `left`
        right = Monomial.rand(state, self.nrows)
        left = Monomial.rand(state, self.nrows)

        # right multiplication
        tmp = monomial_mul(self.rows, right)

        # left multiplication
        self.rows = [tmp[pos].scalar_mul(val) for pos, val in zip(left.perms, left.coeffs)]

    def SortRows(self, z: int, L: Multiset) -> int:
        """
        NOTE: absolutely not constant time.
        :param z[in] number of elements to sort
        :param L[in]: pointer to the currently shortest row
        :return 1 on success
                0 if two rows generate the same multiset
        """
        # first sort each row into a tmp buffer
        tmp = [Multiset.from_row(row) for row in self.rows[:z]]

        # early exit if the shortest row of this matrix is already longer
        # than the currently shortest one
        if max(t[0] for t in tmp) < L[0]:
            return 0

        p = sorted(range(z), key=cmp_to_key(lambda a, b: tmp[a].partial_cmp(tmp[b])))
        self.rows[:z] = [self.rows[i] for i in p]
        return 1

    def SortCols(self) -> int:
        self.rows = [Vector(row) for row in zip(*sorted(zip(*self.rows)))]
        return 1

    def scale_cf_preprocess_v1(self, z: int, m: Multiset) -> int:
        """ by floyd """
        for i in range(z):
            s = self.rows[i].acc()
            if s != 0:
                s = s.inv()
            else:
                s = self.rows[i].acc_inv()
                if s == 0:
                    return 0

            self.rows[i] = self.rows[i].scalar_mul(s)
            # compute the histogram
            if Multiset.from_row(self.rows[i]).partial_cmp(m) < 0:
                return 1
        return 0

    def ScaleCFSubPreprocess(self, z: int, m: Multiset) -> int:
        """ by luke """
        ret = 0
        for i in range(z):
            s = self.rows[i].acc()
            if s != 0:
                s = s.inv()
            else:
                s = self.rows[i].acc_inv()
                if s == 0:
                    continue

            self.rows[i] = self.rows[i].scalar_mul(s)
            # compute the histogram
            tmp = Multiset.from_row(self.rows[i])
            if tmp.partial_cmp(m) < 0:
                ret = 1
                # copy new smallest row
                m[:] = tmp
        return ret

    def SortCF(self, L: Multiset) -> int:
        """
        Implementation of the `SortCF` algorithm from the specification.
        NOTE: computes the result inplace
        NOTE: early exist from the `SortRows` function if `L` is already shorter.
        NOTE: first sort the rows, then the columns
        :param L[in]: pointer the currently shortest row.
        :return 0 on failure (identical rows, which create the same multiset)
                1 on success
        """
        if self.SortRows(self.nrows, L) == 0:
            return 0
        return self.SortCols()

    def ScaleCF(self, L: Multiset) -> int:
        """
        Implementation of the `ScaleCF` algorithm from the specification.
        NOTE: computes the result inplace
        NOTE: does not early exit. Only `SortCF` can do that.
        :param L[in]: pointer the currently shortest row, in histogram form.
        :return 0 on failure:
                     - compute_power_column fails.
                     - identical rows, which create the same multiset
                 1 on success
        """
        for row in range(self.nrows):
            if self.rows[row].all_same():
                continue

            s = self.rows[row].acc()
            if s != 0:
                s = s.inv()
            else:
                s = self.rows[row].acc_inv()
                if s == 0:
                    return 0
            self.rows[row] = self.rows[row].scalar_mul(s)
        return self.SortCF(L)

    def max_zero_rows(self) -> tuple[list[int], list[int]]:
        # count zeros in each row
        zeros = [row.count_zero() for row in self.rows]
        max_zeros = max(zeros)
        return zeros, [row for row, z in enumerate(zeros) if z == max_zeros]

    def scaled_rows(self, row_inv_data: Vector, J) -> 'MatrixNormalized':
        return MatrixNormalized.from_rows([self.rows[j].mul(row_inv_data) for j in J])

    def ImprovedCFBase(self) -> int:
        """
        NOTE: non-constant time
        NOTE: computes the result inplace
        This is the second-fastest implementation of the canonical form function.
        Fallback implementation if `cf` fails
        \\return 0 on failure
                 1 on success
        """
        M = [Vector([Fq(Q - 1)] * self.ncols) for _ in range(self.nrows)]
        touched = 0
        zeros, J = self.max_zero_rows()

        # NOTE: fallback solution if everything falls apart.
        if len(J) == self.nrows:
            return 0

        # NOTE: this is already "sorted"
        L = Multiset.default()

        # Check smallest rows of all matrices to find the smallest candidate
        for row in range(self.nrows):
            if zeros[row] > 0:
                continue

            row_inv_data = self.rows[row].inv_non_ct()
            B = self.scaled_rows(row_inv_data, J)
            if B.scale_cf_preprocess_v1(len(J), L) == 1:
                B = self.scaled_rows(row_inv_data, range(self.nrows))
                if B.ScaleCF(L) == 1 and B.rows < M:
                    L = Multiset.from_row(B.rows[0])
                    touched = 1
                    M = B.rows

        self.rows = M
        return touched

    def cf(self) -> int:
        """
        translation of `ImprovedCF`
        NOTE: non-constant time
        NOTE: computes the result inplace

        This is the fastest implementation of the canonical form function.
        This function scans the matrix for the rows which probably lead
        to the shortest canonical form. If this process fails it falls
        back to either `CFOriginal` or `ImprovedCF`. The very slow `CFOriginal`
        is only called if every row in the matrix contains zeros.
        :return 0 on failure
                1 on success
        """
        smallest_scaling_row = 0
        zeros, J = self.max_zero_rows()

        if len(J) == self.nrows:
            return 0

        # NOTE: this is already "sorted"
        L = Multiset.default()

        # Check smallest rows of all matricies to find the smallest candidate
        for row in range(self.nrows):
            if zeros[row] > 0:
                continue

            B = self.scaled_rows(self.rows[row].inv_non_ct(), J)
            if B.ScaleCFSubPreprocess(len(J), L) == 1:
                smallest_scaling_row = row

        # Calculate CF for best candidate
        B = self.scaled_rows(self.rows[smallest_scaling_row].inv_non_ct(), range(self.nrows))
        ret = B.ScaleCF(L)
        # If candidate was not valid, fall back to regular approach
        if ret != 1:
            # Best scaled row was not valid;
            return self.ImprovedCFBase()

        self.rows = B.rows
        return ret


def test():
    def scale_row(A, i: int) -> Matrix:
        ret = A[i][0]
        for j in range(1, len(A[i])):
            ret += A[i][j]
        t = ret.inv()
        A.scale_row(i, t)
        return A

    def scale_rows(A) -> Matrix:
        for i in range(A.nrows):
            scale_row(A, i)
        return A

    nc, nr, q = 4, 3, 7
    A = Matrix(nr, nc, q)
    # A[0][0], A[0][1], A[0][2], A[0][3] =
    # A[1][0], A[1][1], A[1][2], A[1][3] =
    # A[2][0], A[2][1], A[2][2], A[2][3] =
    A.set([[1, 2, 4, 5], [0, 2, 3, 1], [2, 4, 2, 0]])
    #A.print()

    P = Matrix(nr, nr, q)
    P.set([[0, 1, 0], [4, 0, 0], [0, 0, 3]])

    print(A)
    C = P.mul(A);
    print(C)
    C = scale_rows(C)
    A = scale_rows(A)
    print(A)
    print(C)


def test_less():
    # expand is the exact inverse of compress on an RREF matrix
    a = Matrix.rand_from_seed(bytes(32), K, N)
    is_pivot_column, was_pivot_column = [0] * N, [0] * N
    assert a.rref(is_pivot_column, was_pivot_column, K)
    assert sum(is_pivot_column) == K
    assert all(isinstance(v, Fq) for row in a.rows for v in row)
    flags = [0] * N
    b = Matrix.expand(a.compress(is_pivot_column), flags)
    assert flags == is_pivot_column and b == a

    # the canonical form is invariant under blinding
    state = csprng_initialize(bytes(16))
    A = MatrixNormalized.from_information_set(a, is_pivot_column)
    B = MatrixNormalized.from_rows([Vector(row) for row in A.rows])
    B.blind(state)
    assert B != A
    assert A.cf() == 1 and B.cf() == 1
    assert A == B


if __name__ == "__main__":
    test()
    test_less()
    exit(0)
    nc, nr, q, w = 10, 5, 2, 2
    A = Matrix(nr, nc, q)
    A.print()
    print()

    A.random()
    A.print()
    print()

    rank = A.rref([0] * nc, [0] * nc, 0)
    A.print()
    print("rank", rank)

    e = Matrix(1, nc, q)
    e.random_row_with_weight(0, w)
    e.print()
    print()

    eT = e.transpose()
    eT.print()
    print()

    C = A.mul(eT)
    C.print()
