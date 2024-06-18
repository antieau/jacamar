"""
TEST_MATRICES.

Tests for the Integer, IntegerPython, IntegerRing, and IntegerRingPython
classes.
"""

import pytest
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.rationals import QQ
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.complexes import CC
from nuthatch.rings.polynomials import PolynomialRing
import numpy as np
from nuthatch.matrices.matrices import Matrix, _MatrixGenericData, generate, random
import time
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
        assert self.m._is_generic == False

    def test_classes(self):
        """Tests classes."""
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
        # x = self.m * self.n
        # x_py = self.m_py * self.n_py
        # assert x == Matrix(base_ring=ZZ, entries=[[9, 12, 15], [19, 26, 33]])
        # assert x_py == Matrix(base_ring=ZZ_py, entries=[[9, 12, 15], [19, 26, 33]])
        # s1 = 0
        # s2 = 0
        # s3 = 17
        # a = Matrix(base_ring=RR, entries=np.ones((s1, s1),dtype=int).tolist())
        # b = Matrix(base_ring=RR, entries=np.ones((s2, s2),dtype=int).tolist())
        # c = Matrix(base_ring=RR, entries=np.ones((s3, s3),dtype=int).tolist())
        # assert Matrix(base_ring=RR, entries=(-1*np.ones((s1, s1),dtype=int)).tolist()) == a*RR(-1)
        # assert Matrix(base_ring=ZZ, entries=(-1*np.ones((s1, s1),dtype=int)).tolist()) == a*ZZ(-1)
        # assert a*a == Matrix(base_ring=RR, entries=(s1*np.ones((s1, s1),dtype=int)).tolist())
        # assert b*b == Matrix(base_ring=RR, entries=(s2*np.ones((s2, s2),dtype=int)).tolist())
        # assert c*c == Matrix(base_ring=RR, entries=(s3*np.ones((s3, s3),dtype=int)).tolist())
        s = 2000
        a = random(ZZ, 100, s, s)
        b = random(ZZ, 100, s, s)
        t1 = time.time()
        a * b
        print(time.time() - t1)
        assert 1 == 0

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
        """Tests QQ with matrices."""
        a = Matrix(base_ring=QQ, entries=[[QQ(1, 2), QQ(3, 5)], [QQ(46, 3), QQ(23, 2)]])
        b = Matrix(
            base_ring=QQ, entries=[[QQ(189, 20), QQ(36, 5)], [QQ(184), QQ(2829, 20)]]
        )
        c = Matrix(base_ring=QQ, entries=[[QQ(1), QQ(6, 5)], [QQ(92, 3), QQ(23)]])
        assert a * a == b
        assert a + a == c

    def test_index(self):
        """Tests modified __get_index__."""

        a = Matrix(base_ring=RR, entries=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        b = Matrix(base_ring=RR, entries=[[5]])
        c = Matrix(base_ring=RR, entries=[[3], [6], [9]])
        d = Matrix(base_ring=RR, entries=[[4, 5]])
        assert a[:, :] == a
        assert a[1, 1] == b
        assert a[:, 2] == c
        assert a[1, 0:2] == d

    def test_size(self):
        """Tests .size()"""
        a = Matrix(base_ring=RR, entries=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert a.size() == (3, 3)

    def test_det(self):
        """Tests flint determiant method."""
        a = Matrix(base_ring=RR, entries=[[1, 2], [2, 1]])
        assert a.det() == RR(-3)

    def test_transpose(self):
        """Tests flint transpose method."""
        a = Matrix(base_ring=RR, entries=[[1, 2], [3, 4]])
        b = Matrix(base_ring=RR, entries=[[1, 3], [2, 4]])
        p = random(RR, 2, 5, 4, 1, 2)
        assert p.transpose()[2, 3] == p[3, 2]
        assert a.transpose() == b
        assert a.T() == b

    def test_np_construction(self):
        """Tests numpy construction for RR_py and ZZ_py."""
        a = Matrix(base_ring=RR_py, entries=[[1, 2], [3, 4]])
        print(a*a)
        s = 1000
        p = random(RR_py, 10, s, s)
        print(p)
        # f = random(RR, 10, s, s)
        # t1 = time.time()
        # p*p
        # print(time.time()-t1)
        # print(p*p)
        # t1 = time.time()
        # print(time.time()-t1)

        assert 1 == 0




class TestGenericMatrices:
    """Tests for generic matrices."""
    z = PolynomialRing(base_ring=ZZ, ngens=3, prefix="x", packed=False)
    x0 = z.gens[0]
    x1 = z.gens[1]
    x2 = z.gens[2]

    r = PolynomialRing(base_ring=RR, ngens=3, prefix="x", packed=False)
    y0 = r.gens[0]
    y1 = r.gens[1]
    y2 = r.gens[2]


    # ZZ tests
    def test_generic_creation(self):
        """Test creation of a generic ZZ matrix of polynomials."""
        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert f == ZZ(2) * self.x1 * self.x2 + ZZ(9) * self.x0 ** ZZ(4)
        assert a[1, 1] == Matrix(base_ring=PolynomialRing, ncols=1, nrows=1, entries=[[f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=1, nrows=1, entries=[[f]]))

    def test_generic_mult(self):
        """Tests __mul__ of a generic ZZ matrix."""
        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        s = 127
        a = generate(f, s, s)
        b = generate(ZZ(s)*f*f, s, s)
        assert a * a == b
    
    def test_generic_add(self):
        """Tests __add__ of a generic ZZ matrix."""
        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a + a == Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f+f, f+f], [f+f, f+f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f+f, f+f], [f+f, f+f]]))
    
    def test_generic_sub(self):
        """Tests __sub__ of a generic ZZ matrix."""
        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a + a - a == a

    def test_straussen_mult(self):
        """Tests __mul__ (strassen alogorithm) of a generic ZZ matrix."""

        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        s = 32
        a = generate(f, s, s)
        b = random(RR, 10, s, s, 1, 15)
        t1 = time.time()

        b * b
        print(f'Test ({s}x{s}): {round(time.time()-t1, 2)}s')
        assert 1 == 0



    def test_kmb_mult(self):
        """Tests __mul__ (Kauers-Moosbauer alogorithm) of a generic ZZ matrix."""
        f = self.z({(1, 1, 2, 1): ZZ(2), (0, 4): ZZ(9)})
        s = 5
        a = generate(f, s, s)
        b = generate(ZZ(s)*f*f, s, s)
        print(len((a *a).data.entries))
        assert a * a == b
        
    # RR tests
    def test_generic_creation_RR(self):
        """Test creation of a generic RR matrix of polynomials."""
        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a[1, 1] == Matrix(base_ring=PolynomialRing, ncols=1, nrows=1, entries=[[f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=1, nrows=1, entries=[[f]]))

    def test_generic_mult_RR(self):
        """Tests __mul__ of a generic ZZ matrix."""
        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a * a == Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]]))
    
    def test_generic_add_RR(self):
        """Tests __add__ of a generic ZZ matrix."""
        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a + a == Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f+f, f+f], [f+f, f+f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f+f, f+f], [f+f, f+f]]))

    def test_straussen_mult_RR(self):
        """Tests __mul__ (strassen alogorithm) of a generic ZZ matrix."""

        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a * a == Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]]))

    def test_kmb_mult_RR(self):
        """Tests __mul__ (Kauers-Moosbauer alogorithm) of a generic ZZ matrix."""

        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[f, f], [f, f]]))
        assert a * a == Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[RR(2)*f*f, RR(2)*f*f], [RR(2)*f*f, RR(2)*f*f]]))
    
    def test_generic_det(self):
        """Tests generic implemnetation of the deternimant method for polynomials."""
        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})
        a = generate(f, 2, 2)
        assert a.det() == f * (f**ZZ(2) - f**ZZ(2)) - f * (f**ZZ(2) - f**ZZ(2)) + f * (f**ZZ(2) - f**ZZ(2))
 
    def test_generate(self):
        """Tests generate fuinction, which creates a matrix of size (r, c) with a specific entry."""
        f = self.r({(1, 1, 2, 1): RR(2), (0, 4): RR(9)})

    
    def test_random(self):
        """Tests random matrix generation function."""
        f = random(RR, 10, 5, 5, 1, 3)
        a = random(RR, 10, 10, 10)
        b = random(RR, 2, 10, 10)
        assert a * b, f