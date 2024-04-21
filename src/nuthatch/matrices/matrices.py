"""
MATRICES

A module for dense matrices.
"""

class Matrix:
    def __init__(
        self,
        *,
        ring,
        nrows=None,
        ncols=None,
        entries=None,
        flint_data=None,
    ):
        if (nrows is None or ncols is None) and entries is None:
            raise ValueError("Not enough information provided to construct Matrix.")
        if entries is None:
            self._entries = []
            zero_element = ring(0)
            for i in range(nrows):
                i_row = []
                for j in range(ncols):
                    i_row.append(zero_element)
                self._entries.append(i_row)
            self.nrows = nrows
            self.ncols = ncols
        else:
            # Check it's all well-posed.

        self._flint_type = None

        if ring in {ZZ}:
            self._flint_type = flint.fmpz_mat
            self._flint_data = blah

    def @classmethod
    def _from_flint(cls,flint_matrix):

    def __add__(self,other):
        if self._flint_type:
            return Matrix(self.ring...)
        else:
            ...
