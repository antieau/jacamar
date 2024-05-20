"""
TEST_RATIONALS.

Tests for the Rational and RationalRing classes.
"""

import pytest
from nuthatch.rings.integers import ZZ
from nuthatch.rings.rationals import QQ


class TestRational:
    """Tests for the Rational class."""

    def test_no_args(self):
        """Tests QQ()."""
        assert QQ() == QQ(0)

    def test_bad_args(self):
        """Tests that everything fails with too many arguments."""
        with pytest.raises(TypeError):
            QQ(1, 2, 3)

    def test_truediv(self):
        """Tests that x / y."""
        x = QQ(7)
        y = QQ(11)
        q = x / y
        assert str(q) == "7/11"

    def test_itruediv_notimplemented(self):
        """Tests x /= y."""
        x = QQ(7)
        y = QQ(11)
        x /= y
        assert str(x) == "7/11"

    def test_str(self):
        """Tests __str__."""
        assert str(QQ(7, 11)) == "7/11"
        assert str(QQ) == "The ring of Rational numbers (via flint.fmpq)."

    def test_repr(self):
        """Tests __repr__."""
        assert QQ(7, 11).__repr__() == "7/11"
        assert QQ.__repr__() == "The ring of Rational numbers (via flint.fmpq)."

    def test_add_with_integer(self):
        """Tests __add__ with an integer."""
        assert ZZ(1) + QQ(1, 2) == QQ(3, 2)
        with pytest.raises(TypeError):
            QQ(1, 2) + ZZ(1)
