"""
MATRICES

A module for dense matrices. We create a single abstract matrix class which
holds a `_data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.
"""


class _MatrixGeneric:
    def __init__(self):
        pass


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
            pass

        self._flint_type = None

        if ring in {ZZ}:
            self._flint_type = flint.fmpz_mat
            self._flint_data = blah

    @classmethod
    def ZZ():
        pass

    @classmethod
    def QQ():
        pass

    @classmethod
    def RR():
        pass

    @classmethod
    def CC():
        pass

    @classmethod
    def ZN():
        pass

    @classmethod
    def GF():
        pass

    @classmethod
    def _from_flint(cls, flint_matrix):
        pass

    def __add__(self, other):
        if self._flint_type:
            pass

        pass
