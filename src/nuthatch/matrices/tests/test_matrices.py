"""
TEST_MATRICES.

Tests for the Integer, IntegerPython, IntegerRing, and IntegerRingPython
classes.
"""

import pytest
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.rationals import QQ
from nuthatch.matrices.matrices import Matrix


class TestMatrix:
    """Tests for the Matrix class."""

    m = Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
    q = Matrix(base_ring=ZZ, entries=[[1, 0], [3, 0]])
    n = Matrix(base_ring=ZZ, entries=[[1, 2, 3], [4, 5, 6]])
    thin = Matrix(base_ring=ZZ, nrows=3, ncols=0)
    flat = Matrix(base_ring=ZZ, nrows=0, ncols=2)

    # Generics
    m_py = Matrix(base_ring=ZZ_py, entries=[[1, 2], [3, 4]])
    q_py = Matrix(base_ring=ZZ_py, entries=[[1, 0], [3, 0]])
    n_py = Matrix(base_ring=ZZ_py, entries=[[1, 2, 3], [4, 5, 6]])
    thin_py = Matrix(base_ring=ZZ_py, nrows=3, ncols=0)
    flat_py = Matrix(base_ring=ZZ_py, nrows=0, ncols=2)

    def test_construction_from_nuthatch(self):
        """Tests for construction from Nuthatch elements."""
        assert self.m == Matrix(base_ring=ZZ, entries=[[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]])

    def test_classes(self):
        assert isinstance(self.m, Matrix)

    def test_equality_among_the_classes(self):
        pass

    def test_sizes(self):
        """Tests nrows and ncols."""
        assert self.m.nrows == 2
        assert self.m.ncols == 2
        assert self.n.nrows == 2
        assert self.n.ncols == 3
        assert self.thin.nrows == 3
        assert self.thin.ncols == 0
        assert self.flat.nrows == 0
        assert self.flat.ncols == 2

        assert self.m_py.nrows == 2
        assert self.m_py.ncols == 2
        assert self.n_py.nrows == 2
        assert self.n_py.ncols == 3
        assert self.thin_py.nrows == 3
        assert self.thin_py.ncols == 0
        assert self.flat_py.nrows == 0
        assert self.flat_py.ncols == 2

    def test_add(self):
        """Tests __add__."""
        assert self.m + self.m == Matrix(base_ring=ZZ, entries=[[2, 4], [6, 8]])
        assert self.m_py + self.m_py == Matrix(
            base_ring=ZZ_py, entries=[[2, 4], [6, 8]]
        )
        with pytest.raises(ValueError):
            self.m + self.n
        with pytest.raises(ValueError):
            self.m_py + self.n_py
        with pytest.raises(TypeError):
            self.m + self.m_py
        with pytest.raises(ValueError):
            self.m_py + self.m

    def test_mul(self):
        """Tests __mul__."""
        x = self.m * self.n
        x_py = self.m_py * self.n_py

        assert x == Matrix(base_ring=ZZ, entries=[[9, 12, 15], [19, 26, 33]])
        assert x_py == Matrix(base_ring=ZZ_py, entries=[[9, 12, 15], [19, 26, 33]])

    def test_empty_mul(self):
        """Tests multiplication with empty matrices of various sizes."""

        y = self.flat * self.m
        y_py = self.flat_py * self.m_py

        assert y.nrows == 0
        assert y.ncols == 2
        assert y_py.nrows == 0
        assert y_py.ncols == 2

        assert y == Matrix(base_ring=ZZ, nrows=0, ncols=2)
        assert y_py == Matrix(base_ring=ZZ_py, nrows=0, ncols=2)

        z = self.n * self.thin
        z_py = self.n_py * self.thin_py

        assert z.nrows == 2
        assert z.ncols == 0
        assert z_py.nrows == 2
        assert z_py.ncols == 0

        assert z == Matrix(base_ring=ZZ, nrows=2, ncols=0)
        assert z_py == Matrix(base_ring=ZZ_py, nrows=2, ncols=0)

        w = Matrix(base_ring=ZZ, nrows=0, ncols=0)
        w_py = Matrix(base_ring=ZZ_py, nrows=0, ncols=0)

        assert w * w == w
        assert w_py * w_py == w_py

        assert w.nrows == 0
        assert w.ncols == 0
        assert w_py.nrows == 0
        assert w_py.ncols == 0

    def test_sub(self):
        """Tests __sub__."""
        assert self.m - self.q == Matrix(base_ring=ZZ, entries=[[0, 2], [0, 4]])
        assert self.m_py - self.q_py == Matrix(
            base_ring=ZZ_py, entries=[[0, 2], [0, 4]]
        )
        with pytest.raises(ValueError):
            self.m - self.n
        with pytest.raises(ValueError):
            self.m_py - self.n_py
        with pytest.raises(TypeError):
            self.m - self.m_py
        with pytest.raises(ValueError):
            self.m_py - self.m

    def test_rational(self):
        a = Matrix(base_ring=QQ, entries=[[QQ(1, 2), QQ(3, 5)], [QQ(46, 3), QQ(23, 2)]])
        b = Matrix(
            base_ring=QQ, entries=[[QQ(189, 20), QQ(36, 5)], [QQ(184), QQ(2829, 20)]]
        )
        c = Matrix(base_ring=QQ, entries=[[QQ(1), QQ(6, 5)], [QQ(92, 3), QQ(23)]])
        assert a * a == b
        assert a + a == c


class TestGenericMatrices:
    """Tests for generic matrices."""

    pass
