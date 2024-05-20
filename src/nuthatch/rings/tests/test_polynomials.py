"""
TEST_POLYNOMIALS

Tests for the _Monomial, _Polynomial, Polynomial, and PolynomialRing classes.
"""

import pytest
import flint
from nuthatch.rings.polynomials import (
    _Monomial,
    _Polynomial,
    Polynomial,
    PolynomialRing,
)
from nuthatch.rings.integers import ZZ


class TestMonomial:
    """Tests for the _Monomial class."""

    def test_empty(self):
        """Edge tests for the empty monomial."""
        assert _Monomial(()) == _Monomial(())
        assert _Monomial(()).in_ith_variable(1) == 0
        m = _Monomial((1, 2, 4, 7))
        assert m * _Monomial(()) == m

    def test_eq(self):
        """Tests the __eq__ method."""
        assert _Monomial((1, 2, 4, 7)) == _Monomial((1, 2, 4, 7))
        assert _Monomial((1, 2, 4, 7)) != _Monomial((1, 2, 4, 8))

    def test_in_ith_variable(self):
        """Tests the `in_ith_variable` method."""
        m = _Monomial((1, 2, 4, 7))
        assert m.in_ith_variable(1) == 2
        assert m.in_ith_variable(4) == 7
        assert m.in_ith_variable(0) == 0
        assert m.in_ith_variable(2) == 0
        assert m.in_ith_variable(5) == 0

    def test_mul(self):
        """Tests the __mul__ method."""
        m = _Monomial((1, 2, 4, 7))
        n = _Monomial((0, 3, 3, 3, 4, 1))
        assert m * n == _Monomial((0, 3, 1, 2, 3, 3, 4, 8))
        assert m * _Monomial(()) == m


class TestPolynomialDataClass:
    """Tests for the _Polynomial class."""

    d = {_Monomial((0, 1)): flint.fmpz(1), _Monomial((1, 1)): flint.fmpz(1)}
    p = _Polynomial(ZZ, d)

    def test_mul(self):
        """Tests for the __mul__ method."""
        # (x+y)^2 == x^2 + 2xy + y^2
        assert self.p * self.p == _Polynomial(
            ZZ,
            {
                _Monomial((0, 2)): flint.fmpz(1),
                _Monomial((0, 1, 1, 1)): flint.fmpz(2),
                _Monomial((1, 2)): flint.fmpz(1),
            },
        )


class TestPolynomialRing:
    """Tests for the PolynomialRing class."""

    r = PolynomialRing(base_ring=ZZ, ngens=3, prefix="x")
    x0 = r.gens[0]
    x1 = r.gens[1]
    x2 = r.gens[2]

    def test_gens_data(self):
        """Tests the correct creation of the generators."""
        assert self.r.gens[0].data == _Polynomial(
            ZZ, {_Monomial((0, 1)): flint.fmpz(1)}
        )

    def test_square(self):
        """Tests that (x0+x1)**2 == x0**2 + 2*x0x1 + x1**2."""
        assert (self.x0 + self.x1) ** ZZ(2) == (
            self.x0 ** ZZ(2) + ZZ(2) * self.x0 * self.x1 + self.x1 ** ZZ(2)
        )
