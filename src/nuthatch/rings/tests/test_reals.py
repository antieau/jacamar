"""
TEST_REALS

Tests for the RealNumber and RealRing classes.
"""

import pytest
import flint
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.polynomials import (
    _Monomial,
    _Polynomial,
    PolynomialRing,
)
from nuthatch.matrices.matrices import Matrix


class TestRealNumber:
    """Tests for the RealNumber class."""

    def test_no_args(self):
        """Tests RR()."""
        assert RR() == RR(0)

    def test_truediv(self):
        """Tests that x / y."""
        x = RR(1)
        y = RR(3)
        q = x / y
        assert str(q) == "[0.333333333333333 +/- 3.71e-16]"

    def test_itruediv_notimplemented(self):
        """Tests x /= y."""
        x = RR(1)
        y = RR(3)
        x /= y
        assert str(x) == "[0.333333333333333 +/- 3.71e-16]"

    def test_str(self):
        """Tests __str__."""
        assert str(RR(10.25)) == "10.2500000000000"
        assert str(RR) == "The ring of real numbers (via flint.arb)."

    def test_repr(self):
        """Tests __repr__."""
        assert RR.__repr__() == "The ring of real numbers (via flint.arb)."

    def test_mult(self):
        """Test x * y."""
        assert RR(4.0) * RR(2.25) == RR(9.0)

    def test_add(self):
        """Test x + y."""
        assert RR(4.2) + RR(3.5) == RR(7.7)

    def test_sub(self):
        """Test x - y."""
        assert RR(2.1) - RR(3.1) == RR(-1.0)

    def test_add_with_integer(self):
        """Tests __add__ with an integer."""
        assert ZZ(1) + RR(6.765) == RR(7.765)
        with pytest.raises(TypeError):
            RR(6.765) + ZZ(1)

    def test_mult_with_integer(self):
        """Tests ___mul__ with an integer."""
        assert ZZ(1) * RR(6.765) == RR(6.765)
        with pytest.raises(TypeError):
            RR(6.765) * ZZ(1)

    # def test_mult_with_matrix(self):
    #     """Tests ___mul__ with a matrix."""
    #     mat = Matrix(base_ring=ZZ, entries=[[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]])
    #     assert RR(2.0) * mat == Matrix(base_ring=ZZ, entries=[[ZZ(2), ZZ(4)], [ZZ(6), ZZ(8)]])

    def test_polynomial(self):
        """Tests creation of real polynomials."""
        r = PolynomialRing(base_ring=RR,ngens=2,prefix='x')
        x0,x1=r.gens
        f=x0+x1
        assert f**ZZ(2)==x0**ZZ(2) + ZZ(2)*x0*x1+x1**ZZ(2)

        r = PolynomialRing(base_ring=RR_py,ngens=2,prefix='x')
        x0,x1=r.gens
        f=x0+x1
        assert f**ZZ_py(2)==x0**ZZ_py(2) + ZZ_py(2)*x0*x1+x1**ZZ_py(2)
