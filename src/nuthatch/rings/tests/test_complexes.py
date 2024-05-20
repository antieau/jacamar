"""
TEST_RATIONALS

Tests for the RealNumber and RealRing classes.
"""

import pytest
from nuthatch.rings.complexes import CC


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
        assert str(q) == "[0.33333333333333333333333333333333333333333333333333 +/- 3.78e-51]"

    def test_itruediv_notimplemented(self):
        """Tests x /= y."""
        x = CC(1)
        y = CC(3)
        x /= y
        assert str(x) == "[0.33333333333333333333333333333333333333333333333333 +/- 3.78e-51]"

    def test_str(self):
        """Tests __str__."""
        assert str(CC(2+3j)) == "2.00000000000000 + 3.00000000000000j"
        assert str(CC) == "The ring of complex numbers (via flint.acb)."

    def test_repr(self):
        """Tests __repr__."""
        assert CC.__repr__() == "The ring of complex numbers (via flint.acb)."

a = TestComplexNumber()
a.test_no_args()
a.test_truediv()
a.test_itruediv_notimplemented()
a.test_repr()
a.test_str()