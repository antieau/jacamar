"""
TEST_COMPLEXES

Tests for the RealNumber and RealRing classes.
"""

import pytest
import flint
from nuthatch.rings.complexes import (
    CC,
    sin,
    cos,
    tan,
    asin,
    imag,
    real,
    acos,
    atan,
    csc,
    sec,
    cot,
    abs_lower,
    abs_upper,
    exp,
    log,
    root,
    sqrt,
)
from nuthatch.rings.reals import RR, RR_py
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.polynomials import (
    SparseMonomialData,
    PolynomialData,
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
        poly = PolynomialData(
            CC,
            {
                SparseMonomialData((0, 1)): flint.fmpz(1),
                SparseMonomialData((1, 1)): flint.fmpz(1),
            },
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


class TestComplexNumberFunctions:
    """Tests flint CC functions."""

    def test_sin(self):
        """Tests the sin global function."""
        assert str(sin(CC(1))) == "[0.841470984807897 +/- 6.08e-16]"

    def test_cos(self):
        """Tests the cos global function."""
        assert str(cos(CC(1))) == str(flint.acb(1).cos())

    def test_tan(self):
        """Tests the tan function."""
        assert str(tan(CC(1))) == str(flint.acb(1).tan())

    def test_asin(self):
        """Tests the asin global function."""
        assert str(asin(CC(1))) == str(flint.acb(1).asin())

    def test_acos(self):
        """Tests the acos function."""
        assert str(acos(CC(1))) == str(flint.acb(1).acos())

    def test_atan(self):
        """Tests the atan function."""
        assert str(atan(CC(1))) == str(flint.acb(1).atan())

    def test_csc(self):
        """Tests the csc global function."""
        assert str(csc(CC(1))) == str(flint.acb(1).csc())

    def test_sec(self):
        """Tests the sec function."""
        assert str(sec(CC(1))) == str(flint.acb(1).sec())

    def test_cot(self):
        """Tests the cot function."""
        assert str(cot(CC(1))) == str(flint.acb(1).cot())

    def test_abs_upper(self):
        """Tests the abs_upper global function."""
        assert str(abs_upper(CC("1.5+/-1"))) == str(flint.acb("1.5+/-1").abs_upper())

    def test_abs_lower(self):
        """Tests the abs_lower function."""
        assert str(abs_lower(CC("1.5+/-1"))) == str(flint.acb("1.5+/-1").abs_lower())

    def test_exp(self):
        """Tests the exp global function."""
        assert str(exp(CC(2 + 1j))) == str(flint.acb(2 + 1j).exp())

    def test_log(self):
        """Tests the log function."""
        assert str(log(CC(3.3))) == str(flint.acb(3.3).log())

    def test_root(self):
        """Tests the root function."""
        assert str(root(CC(7 + 2j), 3)) == str(flint.acb(7 + 2j).root(3))

    def test_sqrt(self):
        """Tests the sqrt function."""
        assert str(sqrt(CC(4))) == str(flint.acb(4).sqrt())

    def test_imag(self):
        """Tests the sqrt function."""
        assert str(imag(CC(4 + 3j))) == str(flint.acb(4 + 3j).imag)
        assert type(imag(CC(4 + 3j))) == type(RR(3))

    def test_real(self):
        """Tests the sqrt function."""
        assert str(real(CC(4 + 3j))) == str(flint.acb(4 + 3j).real)
        assert type(real(CC(4 + 3j))) == type(RR(4))
