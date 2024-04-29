"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `_data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.
"""

import flint
from nuthatch.rings.integers import ZZ
from nuthatch.rings.rationals import QQ

class _MatrixGenericData:
    def __init__(self,*,base_ring,nrows,ncols,entries):
        self.base_ring = base_ring
        self.nrows = nrows
        self.ncols = ncols
        self.entries = entries

    def __mul__(self, other):
        if self.ncols != other.nrows:
            raise ValueError(f"Cannot multiply matrix of size {self.nrows}x{self.ncols} with matrix of size {other.nrows}x{other.ncols}.")
        new_entries = [[None] * other.ncols] * self.nrows
        for i in range(self.nrows):
            for j in range(other.ncols):
                new_entry = self.base_ring(0)
                for k in range(self.ncols):
                    new_entry += self.entries[i][k] * other.entries[k][j]
                new_entries[i][j] = new_entry
        return other.__class__(base_ring=other.base_ring,nrows=self.nrows,ncols=other.ncols,entries=new_entries)

    def __str__(self):
        return str(self.entries)

    def __repr__(self):
        return self.__str__()


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
        if not data is None:
            self._data = data
            if self._is_generic:
                self.nrows = data.nrows
                self.ncols = data.ncols
            else:
                self.nrows = data.nrows()
                self.ncols = data.ncols()
        else:
            if (nrows is None or ncols is None) and entries is None:
                raise ValueError("Not enough information provided to construct Matrix.")

            if entries is None:
                self.nrows = nrows
                self.ncols = ncols
                new_entries = self._zero_entries()
            else:
                if type(entries) == dict:
                    if nrows is None or ncols is None:
                        raise ValueError("Not enough information to determine the size of the desired matrix.")
                    self.nrows = nrows
                    self.ncols = ncols
                    new_entries = self._zero_entries()
                    # COERCE the given entries into the base ring, or it's data
                    # class in the non-generic case.
                    for t,e in entries.items():
                        if self._is_generic:
                            new_entries[t[0]][t[1]] = base_ring(e)
                        else:
                            new_entries[t[0]][t[1]] = base_ring(e).data

                # Else, assume that it's a list of lists of base_ring elements.
                else:
                    self.nrows = len(entries)
                    if self.nrows == 0:
                        self.ncols = 0
                    else:
                        self.ncols = len(entries[0])
                    # Now, COERCE the given entries into the base_ring.
                    if self._is_generic:
                        new_entries = [[None] * self.ncols] * self.nrows
                        for i in range(self.nrows):
                            for j in range(self.ncols):
                                new_entries[i][j] = self.base_ring(entries[i][j])                       
                    # Unless, the base_ring is special, in which case we use
                    # the underlying data to later construct a FLINT matrix.
                    else:
                        new_entries = [[None] * self.ncols] * self.nrows
                        for i in range(self.nrows):
                            for j in range(self.ncols):
                                new_entries[i][j] = self.base_ring(entries[i][j]).data

            if self.base_ring == ZZ:
                self._data = flint.fmpz_mat(new_entries)
            elif self.base_ring == QQ:
                self._data = flint.fmpq_mat(new_entries)
            else:
                self._data = _MatrixGenericData(
                    base_ring = self.base_ring,
                    nrows = self.nrows,
                    ncols = self.ncols,
                    entries = new_entries,
                )

    # Constructor helper function.
    def _zero_entries(self):
        entries = []
        if self._is_generic:
            zero_element = self.base_ring(0)
        else:
            zero_element = self.base_ring(0).data

        for i in range(self.nrows):
            i_row = []
            for j in range(self.ncols):
                i_row.append(zero_element)
            entries.append(i_row)
        return entries

    @property
    def data(self):
        return self._data

    def __add__(self, other):
        return other.__class__(
            base_ring = other.base_ring,
            data = self._data + other._data,
        )

    def __mul__(self, other):
        return other.__class__(
            base_ring = other.base_ring,
            data = self._data * other._data,
        )
    
    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self._data.__repr__()
