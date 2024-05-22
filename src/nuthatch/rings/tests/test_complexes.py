"""
TEST_COMPLEXES

Tests for the RealNumber and RealRing classes.
"""

import pytest
import flint
from nuthatch.rings.complexes import CC
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.polynomials import (
    _Monomial,
    _Polynomial,
    Polynomial,
    PolynomialRing,
)


class TestComplexNumber:
    """Tests for the ComplexNumber class."""

    def test_no_args(self):
        """Tests CC()."""
        assert CC() == CC(0)

    def test_truediv(self):
        """Tests that x / y."""
        x = CC(1)
        y = CC(3)
        q = x / y
        assert str(q) == "[0.333333333333333 +/- 3.71e-16]"

    def test_itruediv_notimplemented(self):
        """Tests x /= y."""
        x = CC(1)
        y = CC(3)
        x /= y
        assert str(x) == "[0.333333333333333 +/- 3.71e-16]"

    def test_str(self):
        """Tests __str__."""
        assert str(CC(2 + 3j)) == "2.00000000000000 + 3.00000000000000j"
        assert str(CC) == "The ring of complex numbers (via flint.acb)."

    def test_repr(self):
        """Tests __repr__."""
        assert CC.__repr__() == "The ring of complex numbers (via flint.acb)."

    def test_mult(self):
        """Test x * y."""
        assert CC(4 + 2j) * CC(2) == CC(8 + 4j)

    def test_add(self):
        """Test x + y."""
        assert CC(3 - 2j) + CC(1 + 1j) == CC(4 - 1j)

    def test_sub(self):
        """Test x - y."""
        assert CC(5 - 2j) - CC(2 + 3j) == CC(3 - 5j)

    def test_mult_with_polynomial(self):
        """Tests __mult__ with a polynomial."""
        poly = _Polynomial(
            CC, {_Monomial((0, 1)): flint.fmpz(1), _Monomial((1, 1)): flint.fmpz(1)}
        )

    def test_add_with_integer(self):
        """Tests __add__ with an integer."""
        assert ZZ(1) + CC(6.765) == CC(7.765)
        with pytest.raises(TypeError):
            CC(6.765) + ZZ(1)

    def test_mult_with_integer(self):
        """Tests __mult__ with an integer."""
        assert ZZ(1) * CC(6.765) == CC(6.765)
        with pytest.raises(TypeError):
            CC(6.765) * ZZ(1)

    def test_add_with_real(self):
        """Tests __add__ with an integer."""
        assert RR(1.5) + CC(2 + 1j) == CC(3.5 + 1j)
        with pytest.raises(TypeError):
            CC(2 + 1j) + RR(1.5)

    def test_mult_with_real(self):
        """Tests __mult__ with an integer."""
        assert RR(2) * CC(1 + 1j) == CC(2 + 2j)
        with pytest.raises(TypeError):
            CC(1 + 1j) + RR(2)
