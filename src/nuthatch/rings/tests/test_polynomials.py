"""
TEST_POLYNOMIALS

Tests for the MonomialData, PolynomialData, Polynomial, and PolynomialRing classes.
"""

import pytest
import flint
from nuthatch.rings.polynomials import (
    PackedMonomialData,
    SparseMonomialData,
    PolynomialData,
    Polynomial,
    PolynomialRing,
)
from nuthatch.rings.integers import ZZ


class TestMonomial:
    """Tests for the MonomialData class."""

    def test_empty(self):
        """Edge tests for the empty monomial."""
        assert SparseMonomialData(()) == SparseMonomialData(())
        m = SparseMonomialData((1, 2, 4, 7))
        assert m * SparseMonomialData(()) == m

    def test_eq(self):
        """Tests the __eq__ method."""
        assert SparseMonomialData((1, 2, 4, 7)) == SparseMonomialData((1, 2, 4, 7))
        assert SparseMonomialData((1, 2, 4, 7)) != SparseMonomialData((1, 2, 4, 8))

    def test_mul(self):
        """Tests the __mul__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        n = SparseMonomialData((0, 3, 3, 3, 4, 1))
        assert m * n == SparseMonomialData((0, 3, 1, 2, 3, 3, 4, 8))
        assert m * SparseMonomialData(()) == m

    def test_hash(self):
        """Tests the __hash__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        assert hash(m) == hash((1, 2, 4, 7))

    def test_str_and_repr(self):
        """Tests the __str__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        assert str(m) == str((1, 2, 4, 7))
        assert repr(m) == str(m)



class TestPolynomialDataClass:
    """Tests for the PolynomialData class."""

    d = {SparseMonomialData((0, 1)): flint.fmpz(1), SparseMonomialData((1, 1)): flint.fmpz(1)}
    p = PolynomialData(ZZ, d)

    def test_mul(self):
        """Tests for the __mul__ method."""
        # (x+y)^2 == x^2 + 2xy + y^2
        assert self.p * self.p == PolynomialData(
            ZZ,
            {
                SparseMonomialData((0, 2)): flint.fmpz(1),
                SparseMonomialData((0, 1, 1, 1)): flint.fmpz(2),
                SparseMonomialData((1, 2)): flint.fmpz(1),
            },
        )

    def test_zero(self):
        """Tests the zero function."""
        assert self.p.is_zero() == False
        assert PolynomialData(ZZ,{}).is_zero()
        assert PolynomialData(ZZ,{SparseMonomialData(0,1): ZZ.zero.data}).is_zero()


class TestPolynomialRing:
    """Tests for the PolynomialRing class."""

    r = PolynomialRing(base_ring=ZZ, ngens=3, prefix="x", packed=False)
    x0 = r.gens[0]
    x1 = r.gens[1]
    x2 = r.gens[2]

    def test_gens_data(self):
        """Tests the correct creation of the generators."""
        assert self.r.gens[0].data == PolynomialData(
            ZZ, {SparseMonomialData((0, 1)): flint.fmpz(1)}
        )

    def test_square(self):
        """Tests that (x0+x1)**2 == x0**2 + 2*x0x1 + x1**2."""
        assert (self.x0 + self.x1) ** ZZ(2) == (
            self.x0 ** ZZ(2) + ZZ(2) * self.x0 * self.x1 + self.x1 ** ZZ(2)
        )

    def test_constructor(self):
        """Tests dictionary constructor."""
        f = self.r({(1,1,2,1):ZZ(2),(0,4):ZZ(9)})
        assert f == ZZ(2) * self.x1 * self.x2 + ZZ(9) * self.x0**ZZ(4)


class TestCalls:
    """Test evaluation of polynomials."""
    s=PolynomialRing(base_ring=ZZ,ngens=4,prefix='x', packed=False)
    x0,x1,x2,x3 = s.gens
    f = x0+ZZ(2)*x1*x2+x3

    def test_evaluation_at_coorindates(self):
        """Tests evaluation at tuples of base_ring elements."""
        assert self.f(ZZ(1), ZZ(2), ZZ(3), ZZ(4)) == ZZ(17)

    def test_evaluation_at_monomial(self):
        """Tests evaluation at a SparseMonomialData instance."""
        assert self.f(SparseMonomialData((1, 1, 2, 1))) == ZZ(2)
        assert self.f(SparseMonomialData(1, 1, 2, 1)) == ZZ(2)

    def test_evaluation_at_implicit_monomial(self):
        """
        Tests evaluation at something that can be coerced into a
        SparseMonomialData instance.
        """
        assert self.f(1,1,2,1) == ZZ(0)
