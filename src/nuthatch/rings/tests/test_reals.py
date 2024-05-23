"""
TEST_REALS

Tests for the RealNumber and RealRing classes.
"""

import pytest
import flint
from nuthatch.rings.reals import (
    RR, RR_py, sin, cos, tan, asin, 
    acos, atan, csc, sec, cot, abs_lower,
    abs_upper, exp, factorial, log, root, sqrt
)
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.polynomials import (
    _Monomial,
    _Polynomial,
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
        """Tests __mult__ with an integer."""
        assert ZZ(1) * RR(6.765) == RR(6.765)
        with pytest.raises(TypeError):
            RR(6.765) * ZZ(1)

    # def test_mult_with_matrix(self):
    #     """Tests __mult__ with a matrix."""
    #     mat = Matrix(base_ring=ZZ, entries=[[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]])
    #     assert RR(2.0) * mat == Matrix(base_ring=ZZ, entries=[[ZZ(2), ZZ(4)], [ZZ(6), ZZ(8)]])

    # def test_mult_with_polynomial(self):
    #     """Tests __mult__ with a polynomial."""
    #     poly = _Polynomial(RR, {_Monomial((0, 1)): flint.fmpz(1), _Monomial((1, 1)): flint.fmpz(1)})
    #     assert RR(2.5) * poly == _Polynomial(RR, {_Monomial((0, 1)): flint.arb(2.5), _Monomial((1, 1)): flint.arb(2.5)})


class TestRealNumberFunctions:
    """Tests flint RR functions."""
    def test_sin(self):
        """Tests the sin global function."""
        assert str(sin(RR(1))) == '[0.841470984807897 +/- 6.08e-16]'

    def test_cos(self):
        """Tests the cos global function."""
        assert str(cos(RR(1))) == str(flint.arb(1).cos())

    def test_tan(self):
        """Tests the tan function."""
        assert str(tan(RR(1))) == str(flint.arb(1).tan())

    def test_asin(self):
        """Tests the asin global function."""
        assert str(asin(RR(1))) == str(flint.arb(1).asin())

    def test_acos(self):
        """Tests the acos function."""
        assert str(acos(RR(1))) == str(flint.arb(1).acos())

    def test_atan(self):
        """Tests the atan function."""
        assert str(atan(RR(1))) == str(flint.arb(1).atan())

    def test_csc(self):
        """Tests the csc global function."""
        assert str(csc(RR(1))) == str(flint.arb(1).csc())

    def test_sec(self):
        """Tests the sec function."""
        assert str(sec(RR(1))) == str(flint.arb(1).sec())

    def test_cot(self):
        """Tests the cot function."""
        assert str(cot(RR(1))) == str(flint.arb(1).cot())
        
    def test_abs_upper(self):
        """Tests the abs_upper global function."""
        assert str(abs_upper(RR("1.5+/-1"))) == str(flint.arb("1.5+/-1").abs_upper())

    def test_abs_lower(self):
        """Tests the abs_lower function."""
        assert str(abs_lower(RR("1.5+/-1"))) == str(flint.arb("1.5+/-1").abs_lower())

    def test_factorial(self):
        """Tests the factorial function."""
        assert factorial(RR(5))==RR(120)
        
    def test_exp(self):
        """Tests the exp global function."""
        assert str(exp(RR(2))) == str(flint.arb(2).exp())

    def test_log(self):
        """Tests the log function."""
        assert str(log(RR(3.3))) == str(flint.arb(3.3).log())

    def test_root(self):
        """Tests the root function."""
        assert str(root(RR(7), 3)) == str(flint.arb(7).root(3))

    def test_sqrt(self):
        """Tests the sqrt function."""
        assert str(sqrt(RR(4))) == str(flint.arb(4).sqrt())
