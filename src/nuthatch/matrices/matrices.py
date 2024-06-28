"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.

AUTHORS:
- Benjamin Antieau (2024): initial version.
- Luca Zerega (2024): major revisions.
"""

import flint
import numpy as np
import ray
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.complexes import CC
from nuthatch.rings.rationals import QQ


class _MatrixGenericData:
    """
    This class must be called with all arguments. It is assumed, but not
    checked, that the shape of `entries` conforms with this. In the case that
    `nrows==0` or `ncols==0`, one should have `entries=[]`.
    """

    def __init__(self, *, base_ring, nrows, ncols, entries):
        self.base_ring = base_ring
        self.nrows = nrows
        self.ncols = ncols
        self.entries = entries

    def det(self):
        """Alias for `determinant` method."""
        return self.determinant()

    def determinant(self):
        if self.nrows != self.nrows:
            raise ValueError("Matrix must be square.")
        if self.nrows == 0:
            return self.base_ring(0)
        if self.nrows == 1:
            return self.entries[0][0]
        entries = self.entries
        def minor(array,i,j):
            c = array
            c = c[:i] + c[i+1:]
            for k in range(0,len(c)):
                c[k] = c[k][:j]+c[k][j+1:]
            return c
        def det(array,n):
            if n == 1 :
                return array[0][0]
            if n == 2 :
                return array[0][0]*array[1][1] - array[0][1]*array[1][0]
            det = self.base_ring(base_ring=RR, ngens=3, prefix="x", packed=False)
            det = det.zero
            for i in range(0,n):
                m = minor(array,0,i)
                det = det + RR((-1)**i) * array[0][i] * det(m,n-1)
            return det
        return det(entries, self.nrows)

    def transpose(self):
        """Returns copy of the transposed matrix data."""
        entries = self.entries
        new_entries = []
        for i in range(self.ncols):
            new_entries.append([])
            for j in range(self.nrows):
                new_entries[-1].append(entries[j][i])

        return _MatrixGenericData(
            base_ring=self.base_ring,
            nrows=self.ncols,
            ncols=self.nrows,
            entries=new_entries
        )
    def size(self):
        """Returns size of a matrix as a tuple (rows, cols)."""
        return (self.nrows, self.ncols)

    def __add__(self, other):
        if self.nrows != other.nrows or self.ncols != other.ncols:
            raise ValueError(
                f"Cannot add matrix of size {self.nrows}x{self.ncols} to matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=other.nrows,
                ncols=other.ncols,
                entries=[],
            )
        new_entries = [self.entries[i].copy() for i in range(self.nrows)]
        for i in range(self.nrows):
            for j in range(self.ncols):
                new_entries[i][j] += other.entries[i][j]

        return other.__class__(
            base_ring=other.base_ring,
            nrows=other.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __mul__(self, other):
        if (str(type(other)) == "<class 'nuthatch.rings.reals.RealNumber'>"
            or str(type(other)) == "<class 'nuthatch.rings.complexes.ComplexNumber'>"
            or str(type(other)) == "<class 'nuthatch.rings.integers.Integer'>"
            or str(type(other)) == "<class 'nuthatch.rings.rationals.Rational'>"):
            new_entries = []
            for i in self.entries:
                new_entries.append([])
                for j in i:
                    new_entries[-1].append(other * j)

            return self.__class__(
                base_ring=self.base_ring,
                nrows=self.nrows, 
                ncols=self.ncols,
                entries=new_entries
            )
        if self.ncols != other.nrows:
            raise ValueError(
                f"Cannot multiply matrix of size {self.nrows}x{self.ncols} with matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:

            return other.__class__(
                base_ring=other.base_ring,
                nrows=self.nrows,
                ncols=other.ncols,
                entries=[],
            )

        # STRASSEN
        elif (self.ncols % 2 == 0
              and self.nrows % 2 == 0
              and other.ncols % 2 == 0
              and other.nrows % 2 == 0
              and self.size() == other.size()
              and self.size()[0] > 250):
            @ray.remote
            def strassen(A, B):
                n = len(A.entries)
                if n <= 2:  # Base case
                    return B.__class__(
                        base_ring=B.base_ring,
                        nrows=n,
                        ncols=n,
                        entries=np.dot(A.entries, B.entries),
                    )

                # Partitions
                mid = n // 2
                A11 = A[:mid, :mid]
                A12 = A[:mid, mid:]
                A21 = A[mid:, :mid]
                A22 = A[mid:, mid:]
                B11 = B[:mid, :mid]
                B12 = B[:mid, mid:]
                B21 = B[mid:, :mid]
                B22 = B[mid:, mid:]
                # Recursions
                futures = [
                    strassen(A11, B12 - B22),
                    strassen(A11 + A12, B22),
                    strassen(A21 + A22, B11),
                    strassen(A22, B21 - B11),
                    strassen(A11 + A22, B11 + B22),
                    strassen(A12 - A22, B21 + B22),
                    strassen(A11 - A21, B11 + B12)
                    ]
                
                P1 = futures[0]
                P2 = futures[1]
                P3 = futures[2]
                P4 = futures[3]
                P5 = futures[4]
                P6 = futures[5]
                P7 = futures[6]

                # Combine results to form C
                C11 = (P5 + P4 - P2 + P6).entries
                C12 = (P1 + P2).entries
                C21 = (P3 + P4).entries
                C22 = (P5 + P1 - P3 - P7).entries

                C1 = np.hstack((np.array(C11), np.array(C12)))
                C2 = np.hstack((np.array(C21), np.array(C22)))
                C = np.vstack((C1, C2)).tolist()

                # Combine quadrants to form C
                return _MatrixGenericData(
                    base_ring=A.base_ring,
                    nrows=len(C),
                    ncols=len(C[0]),
                    entries=C
                )
            
            # Set up paralleization
            mid = self.nrows // 2
            A11 = self[:mid, :mid]
            A12 = self[:mid, mid:]
            A21 = self[mid:, :mid]
            A22 = self[mid:, mid:]
            B11 = other[:mid, :mid]
            B12 = other[:mid, mid:]
            B21 = other[mid:, :mid]
            B22 = other[mid:, mid:]
            # Recursions (parallel)

            P1 = strassen(A11, B12 - B22)
            P2 = strassen(A11 + A12, B22)
            P3 = strassen(A21 + A22, B11)
            P4 = strassen(A22, B21 - B11)
            P5 = strassen(A11 + A22, B11 + B22)
            P6 = strassen(A12 - A22, B21 + B22)
            P7 = strassen(A11 - A21, B11 + B12)

            # Combine results to form C
            C11 = (P5 + P4 - P2 + P6).entries
            C12 = (P1 + P2).entries
            C21 = (P3 + P4).entries
            C22 = (P5 + P1 - P3 - P7).entries

            C1 = np.hstack((np.array(C11), np.array(C12)))
            C2 = np.hstack((np.array(C21), np.array(C22)))
            C = np.vstack((C1, C2)).tolist()

            # Combine quadrants to form C
            return _MatrixGenericData(
                base_ring=other.base_ring,
                nrows=len(C),
                ncols=len(C[0]),
                entries=C
            )

        # standard matrix multiplication
        else:
            new_entries: list = []
            for i in range(self.nrows):
                new_entries.append([])
                for j in range(other.ncols):
                    new_entry = self.entries[i][0] * other.entries[0][j]
                    for k in range(1, self.ncols):
                        new_entry += self.entries[i][k] * other.entries[k][j]
                    new_entries[i].append(new_entry)


        return other.__class__(
            base_ring=other.base_ring,
            nrows=self.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __sub__(self, other):
        if self.nrows != other.nrows or self.ncols != other.ncols:
            raise ValueError(
                f"Cannot add matrix of size {self.nrows}x{self.ncols} to matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=other.nrows,
                ncols=other.ncols,
                entries=[],
            )

        new_entries = [self.entries[i].copy() for i in range(self.nrows)]
        for i in range(self.nrows):
            for j in range(self.ncols):
                new_entries[i][j] -= other.entries[i][j]

        return other.__class__(
            base_ring=other.base_ring,
            nrows=other.nrows,
            ncols=other.ncols,
            entries=new_entries,
        )

    def __str__(self):
        return str(self.entries)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self.nrows != other.nrows or self.ncols != other.ncols:
            return False
        if self.nrows == 0 or self.ncols == 0:
            return self.base_ring == other.base_ring
        return self.entries == other.entries

    def __getitem__(self, args):
        r, c = args
        entries = self.entries
        new_entries = entries[r]
        new_data = []

        if not isinstance(new_entries[0], list):
            new_data = [new_entries[c]]
            if not isinstance(new_data[0], list):
                new_data = [new_data]
        else:
            for i in new_entries:
                if not isinstance(i[c], list):
                    new_data.append([i[c]])
                else:
                    new_data.append(i[c])

        ncols = len(new_data[0])
        nrows = len(new_data)
        if nrows == 1 and ncols == 1:
            return entries[0][0]
        
        return self.__class__(
            base_ring=self.base_ring,
            nrows=nrows,
            ncols=ncols,
            entries=new_data,
        )



class _MatrixPythonData:
    def __init__(self, *, base_ring, shape, entries):
        self.base_ring = base_ring
        self.shape = shape
        self.entries = entries



class Matrix:
    """
    Base class for matrices of numbers or polynomials.
    Numerical matrices are built off of flint mat objects
    while generic (polynomial) matrices use the 
    _GenericMatrixData format.
    """
    def __init__(
        self,
        *,
        base_ring,
        nrows=None,
        ncols=None,
        entries=None,
        data=None,
    ):
        
        self.base_ring = base_ring
        if self.base_ring in {ZZ, QQ, RR, CC, RR_py, ZZ_py}:
            self._is_generic = False
            self._is_python = False

            if self.base_ring in {RR_py, ZZ_py}:
                self._is_python = True
        else:
            self._is_generic = True
            self._is_python = False

        # If data is provided, there is a fast constructor.
        if data is not None:
            self.data = data
            if self._is_generic:
                self.nrows = data.nrows
                self.ncols = data.ncols
            else:
                if self._is_python:
                    self.nrows = data.shape[0]
                    self.ncols = data.shape[1]
                else:

                    self.nrows = data.nrows()
                    self.ncols = data.ncols()
        # Else, if the rows or columns are zero, construct an empty matrix of
        # the appropriate size.
        elif nrows == 0 or ncols == 0:
            if nrows is None or ncols is None:
                raise ValueError("Not enough information provided to construct Matrix.")
            self.nrows = nrows
            self.ncols = ncols

            if self.base_ring == ZZ:
                self.data = flint.fmpz_mat(self.nrows, self.ncols)
            elif self.base_ring == QQ:
                self.data = flint.fmpq_mat(self.nrows, self.ncols)
            elif self.base_ring == RR:
                self.data = flint.arb_mat(self.nrows, self.ncols)
            elif self.base_ring == CC:
                self.data = flint.acb_mat(self.nrows, self.ncols)
            elif self.base_ring == ZZ_py:
                self.data = np.ndarray([self.nrows, self.ncols])
            elif self.base_ring == RR_py:
                self.data = np.ndarray([self.nrows, self.ncols])


            else:
                self.data = _MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    entries=[],
                )
        else:
            if (nrows is None or ncols is None) and entries is None:
                raise ValueError("Not enough information provided to construct Matrix.")

            if entries is None:
                self.nrows = nrows
                self.ncols = ncols
                new_entries = self._zero_entries()
            else:
                if isinstance(entries, dict):
                    if nrows is None or ncols is None:
                        raise ValueError(
                            "Not enough information to determine the size of the desired matrix."
                        )
                    self.nrows = nrows
                    self.ncols = ncols
                    new_entries = self._zero_entries()
                    # COERCE the given entries into the base ring, or it's data
                    # class in the non-generic case.
                    for t, e in entries.items():
                        if self._is_generic:
                            new_entries[t[0]][t[1]] = base_ring(e)
                        else:
                            new_entries[t[0]][t[1]] = base_ring(e).data

                # Else, assume that it is a list of lists of base_ring elements.
                else:
                    self.nrows = len(entries)
                    if self.nrows == 0:
                        self.ncols = 0
                    else:
                        self.ncols = len(entries[0])
                    # Now, COERCE the given entries into the base_ring.
                    if self._is_generic:
                        new_entries = []
                        for i in range(self.nrows):
                            new_entries.append([])
                            for j in range(self.ncols):
                                new_entries[i].append(self.base_ring(entries[i][j]))
                    # Unless, the base_ring is special, in which case we use
                    # the underlying data to later construct a FLINT matrix.
                    else:
                        new_entries = []
                        for i in range(self.nrows):
                            new_entries.append([])
                            for j in range(self.ncols):
                                # TODO: fix the coercion here. If the input
                                # entries happen to already be NUTHATCH
                                # classes wrapping FLINT classes, then the
                                # following code does not work!
                                new_entries[i].append(
                                    self.base_ring(entries[i][j]).data
                                )

            if self.base_ring == ZZ:
                self.data = flint.fmpz_mat(new_entries)
            elif self.base_ring == QQ:
                self.data = flint.fmpq_mat(new_entries)
            elif self.base_ring == RR:
                self.data = flint.arb_mat(new_entries)
            elif self.base_ring == CC:
                self.data = flint.acb_mat(new_entries)
            elif self.base_ring == ZZ_py:
                self.data = np.array(new_entries)
            elif self.base_ring == RR_py:
                self.data = np.array(new_entries)

            else:
                self.data = _MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    entries=new_entries,
                )

    # Constructor helper function.
    def _zero_entries(self):
        entries = []
        for i in range(self.nrows):
            i_row = []
            for j in range(self.ncols):
                if self._is_generic:
                    i_row.append(self.base_ring(0))
                else:
                    i_row.append(self.base_ring(0).data)
            entries.append(i_row)
        return entries

    def det(self):
        """Alias for `determinenant` method."""
        return self.determinant()

    def determinant(self):
        """Determinant method for matrices."""
        if self.nrows != self.ncols:
            raise ValueError("Matrix must be square.")
        if not self._is_generic:
            if self._is_python:
                return self.base_ring(np.linalg.det(self.data))
            return self.base_ring(self.data.det())
        return self.data.det()

    def T(self):
        """Alias for the transpose() method."""
        return self.transpose()

    def transpose(self):
        """Returns a transposed copy of self."""
        return Matrix(base_ring=self.base_ring,
                      nrows=self.ncols,
                      ncols=self.nrows,
                      data=self.data.transpose())

    def size(self):
        """Returns size of a matrix as a tuple (rows, cols)."""
        return (self.nrows, self.ncols)

    def __add__(self, other):
        """Returns self + other with base ring that of other."""
        return other.__class__(
            base_ring=other.base_ring,
            data=self.data + other.data,
        )

    def __sub__(self, other):
        """Returns self - other with base ring that of other."""
        return other.__class__(
            base_ring=other.base_ring,
            data=self.data - other.data,
        )

    def __mul__(self, other):
        if (str(type(other)) == "<class 'nuthatch.rings.reals.RealNumber'>"
            or str(type(other)) == "<class 'nuthatch.rings.complexes.ComplexNumber'>"
            or str(type(other)) == "<class 'nuthatch.rings.integers.Integer'>"
            or str(type(other)) == "<class 'nuthatch.rings.rationals.Rational'>"):
            if str(self.base_ring) == "<class 'nuthatch.rings.polynomials.PolynomialRing'>":
                return self.__class__(
                    base_ring=self.base_ring,
                    nrows=self.nrows,
                    ncols=self.ncols,
                    data=self.data * other
                )

            return self.__class__(
                base_ring=self.base_ring,
                nrows=self.nrows,
                ncols=self.ncols,
                data=self.data * other.data
            )

        if self.ncols != other.nrows:
            raise ValueError(
                f"Cannot multiply matrix of size {self.nrows}x{self.ncols} with matrix of size {other.nrows}x{other.ncols}."
            )

        if self.nrows == 0 or other.nrows == 0 or other.ncols == 0:
            return other.__class__(
                base_ring=other.base_ring,
                nrows=self.nrows,
                ncols=other.ncols,
            )

        return other.__class__(
            base_ring=other.base_ring,
            data=self.data * other.data,
        )

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __eq__(self, other):
        if self._is_python and other._is_python:
            return self.data.all() == other.data.all()
        return self.data == other.data
    

    def __getitem__(self, args):
        if not isinstance(args, tuple):
            return ValueError (
                'The matrix slice method takes 2 args [rows, columns], but 1 were given.'
            )

        elif len(args) == 2:
            r, c = args
            if not self._is_generic:
                r, c = args
                entries = self.data.tolist()
            else:
                entries = self.data.entries
            new_entries = entries[r]
            new_data = []

            if not isinstance(new_entries[0], list):
                new_data = [new_entries[c]]
                if not isinstance(new_data[0], list):
                    new_data = [new_data]
            else:
                for i in new_entries:
                    if not isinstance(i[c], list):
                        new_data.append([i[c]])
                    else:
                        new_data.append(i[c])

            ncols = len(new_data[0])
            nrows = len(new_data)
            if not self._is_generic:
                r, c = args
                entries = self.data.tolist()
                return self.__class__(
                base_ring=self.base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=new_data,
                )

            return self.__class__(
                base_ring=self.base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=new_data,
                data=_MatrixGenericData(
                    base_ring=self.base_ring,
                    nrows=nrows,
                    ncols=ncols,
                    entries=new_data,
                    )
                )

        return ValueError (
            f'The matrix slice method takes 2 args [rows, columns], but {len(args)} were given.'
        )


# Functions for matricies

def generate(value, nrows, ncols):
    """
    Generates a repeating matrix of size (rows, cols). Can be generic or not.
    """
    entries = []
    for i in range(nrows):
        entries.append([])
        for j in range(ncols):
            entries[-1].append(value)

    return Matrix(
                base_ring=value.ring,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                data=_MatrixGenericData(
                    base_ring=value.ring,
                    nrows=nrows,
                    ncols=ncols,
                    entries=entries,
                    )
                )

def random(base_ring, max_val, nrows, ncols, poly=0, poly_ring=0, ngens=1):
    """
    Generates a random matrix of size nrows x ncols using np.random.
    Values go from [0, max_val]. 
    Polynomials are of form n0*x0**r0 + n1*x1**r1 ....
    """
    entries = []

    if poly:
        p = poly_ring(base_ring=base_ring, ngens=ngens, prefix="x", packed=True)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                val = p(0)
                for k in range(ngens):
                    if base_ring == RR:
                        val = val + base_ring(np.random.random() * max_val) * p.gens[k] ** ZZ(int(np.random.random() * max_val))
                    elif base_ring == ZZ:
                        val = val + base_ring(int(np.random.random() * max_val)) * p.gens[k] ** ZZ(int(np.random.random() * max_val))

                entries[-1].append(val)


        return Matrix(
                base_ring=poly_ring,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                data=_MatrixGenericData(
                    base_ring=poly_ring,
                    nrows=nrows,
                    ncols=ncols,
                    entries=entries,
                    )
                )
    if base_ring == RR:
        vals = np.random.random(nrows * ncols) * max_val
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(base_ring(vals[(i + 1) * (j)]))

    elif base_ring == ZZ:
        vals = np.random.randint(max_val, size=nrows * ncols)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(base_ring(int(vals[(i + 1) * (j)])))

    elif base_ring == RR_py:
        vals = np.random.random(nrows * ncols) * max_val
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(vals[(i + 1) * (j)])

        return Matrix(
                    base_ring=base_ring,
                    entries=entries,
                    )
    elif base_ring == ZZ_py:
        vals = np.random.randint(max_val, size=nrows * ncols)
        for i in range(nrows):
            entries.append([])
            for j in range(ncols):
                entries[-1].append(int(vals[(i + 1) * (j)]))

        return Matrix(
                    base_ring=base_ring,
                    entries=entries,
                    )

    return Matrix(
                base_ring=base_ring,
                nrows=nrows,
                ncols=ncols,
                entries=entries,
                )