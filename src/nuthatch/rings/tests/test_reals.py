"""
TEST_RATIONALS

Tests for the RealNumber and RealRing classes.
"""

import pytest
from nuthatch.rings.reals import RR


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

a = TestRealNumber()
a.test_no_args()
a.test_truediv()
a.test_itruediv_notimplemented()
a.test_repr()
a.test_str()