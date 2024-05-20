"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.
"""

import flint
from nuthatch.rings.integers import ZZ
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
            raise ValueError("matrix must be square")
        if self.nrows == 0:
            return self.base_ring(0)
        if self.nrows == 1:
            return self.entries[0][0]
        # TODO: implement cofactor definition of determinant for starters.
        return NotImplemented

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

        new_entries = []
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
        # TODO: when should two empty matrices be equal?
        if self.nrows != other.nrows or self.ncols != other.ncols:
            return False
        if self.nrows == 0 or self.ncols == 0:
            return self.base_ring == other.base_ring
        return self.entries == other.entries


class Matrix:
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
        if self.base_ring in {ZZ, QQ}:
            self._is_generic = False
        else:
            self._is_generic = True

        # If data is provided, there is a fast constructor.
        if data is not None:
            self.data = data
            if self._is_generic:
                self.nrows = data.nrows
                self.ncols = data.ncols
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
        if self.nrows != self.ncols:
            raise ValueError("matrix must be square")
        return self.base_ring(self.data.det())

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
        return self.data == other.data
