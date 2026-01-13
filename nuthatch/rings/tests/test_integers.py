"""
TEST_INTEGERS

Tests for the Integer, IntegerPython, IntegerRing, and IntegerRingPython
classes.
"""

import pytest
from flint import fmpz
from flint.utils.flint_exceptions import DomainError
from nuthatch.rings.integers import ZZ, ZZ_py, Integer, IntegerPython, gcd, factor, factorial


class TestInteger:
    """Tests for the Integer class."""

    def test_rings(self):
        assert ZZ(6027945939101000).ring == ZZ
        assert ZZ_py(6027945939101000).ring == ZZ_py

    def test_classes(self):
        assert Integer.data_class == fmpz
        assert IntegerPython.data_class == int

    def test_equality_among_the_classes(self):
        """Tests that ZZ(n) == ZZ_py(n) == Integer(n) == Integer_py(n)."""
        with pytest.raises(AttributeError):
            fmpz(159728757) == ZZ(159728757)
        assert ZZ(159728757) == ZZ_py(159728757)
        assert ZZ(159728757) == Integer(ZZ, 159728757)
        assert ZZ_py(159728757) == IntegerPython(ZZ_py, 159728757)
        with pytest.raises(AttributeError):
            159728757 == ZZ(159728757)
        assert ZZ(159728757).data == fmpz(159728757)

    def test_add(self):
        """Tests __add__."""
        assert ZZ(61146125) + ZZ(98582632) == ZZ(159728757)
        assert ZZ_py(61146125) + ZZ_py(98582632) == ZZ_py(159728757)
        with pytest.raises(AttributeError):
            ZZ(7) + 7

    def test_iadd(self):
        """Tests __iadd__."""
        x = ZZ(61146125)
        x += ZZ(98582632)
        assert x == ZZ(159728757)

        x_py = ZZ_py(61146125)
        x_py += ZZ_py(98582632)
        assert x == ZZ_py(159728757)

    def test_sub(self):
        """Tests __sub__."""
        assert ZZ(61146125) - ZZ(98582632) == ZZ(-37436507)
        assert ZZ_py(61146125) - ZZ_py(98582632) == ZZ_py(-37436507)
        with pytest.raises(AttributeError):
            ZZ(7) - 7

    def test_isub(self):
        """Tests __isub__."""
        x = ZZ(61146125)
        x -= ZZ(98582632)
        assert x == ZZ(-37436507)

        x_py = ZZ_py(61146125)
        x_py -= ZZ_py(98582632)
        assert x == ZZ_py(-37436507)

    def test_neg(self):
        """Tests __neg__."""
        assert -ZZ(98582632) == ZZ(-98582632)
        assert -ZZ_py(98582632) == ZZ_py(-98582632)

    def test_mul(self):
        """Tests __mul__."""
        assert ZZ(61146125) * ZZ(98582632) == ZZ(6027945939101000)
        assert ZZ_py(61146125) * ZZ_py(98582632) == ZZ_py(6027945939101000)
        with pytest.raises(AttributeError):
            ZZ(7) * 7

    def test_imul(self):
        """Tests __imul__."""
        x = ZZ(61146125)
        x *= ZZ(98582632)
        assert x == ZZ(6027945939101000)

        x_py = ZZ_py(61146125)
        x_py *= ZZ_py(98582632)
        assert x == ZZ_py(6027945939101000)

    def test_truediv_notimplemented(self):
        """Tests that x/y returns NotImplemented."""
        with pytest.raises(TypeError):
            ZZ(7) / ZZ(7)
        with pytest.raises(TypeError):
            ZZ_py(7) / ZZ_py(7)

    def test_itruediv_notimplemented(self):
        """Tests that x /= y returns NotImplemented."""
        with pytest.raises((DomainError, TypeError)):
            x = ZZ(7)
            x /= ZZ(7)
        with pytest.raises(TypeError):
            x_py = ZZ_py(7)
            x_py /= ZZ_py(7)

    def test_divmod(self):
        """Tests integer division and remainder."""
        assert divmod(ZZ(98582632), ZZ(61146125)) == (ZZ(1), ZZ(37436507))
        assert divmod(ZZ_py(98582632), ZZ_py(61146125)) == (ZZ_py(1), ZZ_py(37436507))
        with pytest.raises(TypeError):
            ZZ(7) / 7

    def test_floordiv(self):
        """Tests the integer __floordiv__ method."""
        assert ZZ(7) // ZZ(2) == ZZ(3)
        assert ZZ_py(7) // ZZ_py(2) == ZZ_py(3)
        assert ZZ(7) // ZZ_py(2) == ZZ_py(3)
        with pytest.raises(AttributeError):
            ZZ(7) // 2

    def test_ifloordiv(self):
        """Tests __ifloordiv__."""
        x = ZZ(7)
        x //= ZZ(2)
        assert x == ZZ(3)

        x_py = ZZ_py(7)
        x_py //= ZZ_py(2)
        assert x == ZZ_py(3)

    def test_mod(self):
        """Tests __mod___."""
        assert ZZ(11) % ZZ(7) == ZZ(4)
        assert ZZ_py(11) % ZZ_py(7) == ZZ_py(4)
        with pytest.raises(AttributeError):
            ZZ(7) % 3

    def test_imod(self):
        """Tests __imod__."""
        x = ZZ(11)
        x %= ZZ(7)
        assert x == ZZ(4)

        x_py = ZZ_py(11)
        x_py %= ZZ_py(7)
        assert x == ZZ_py(4)

    def test_pow(self):
        """Tests powering."""
        with pytest.raises(AttributeError):
            ZZ(415) ** 6
        assert ZZ(415) ** ZZ(6) == ZZ(5108443333890625)
        assert ZZ_py(415) ** ZZ_py(6) == ZZ_py(5108443333890625)

    def test_ipow(self):
        """Tests __ipow__."""
        x = ZZ(415)
        x **= ZZ(6)
        assert x == ZZ(5108443333890625)

        x_py = ZZ_py(415)
        x_py **= ZZ_py(6)
        assert x == ZZ_py(5108443333890625)

    def test_mixing(self):
        """Test correct type of output."""
        assert isinstance(ZZ(7) + ZZ_py(9), IntegerPython)
        assert isinstance(ZZ(7) * ZZ_py(9), IntegerPython)
        assert isinstance(ZZ(7) - ZZ_py(9), IntegerPython)
        assert isinstance(ZZ(7) ** ZZ_py(9), Integer)
        assert isinstance(ZZ(11) % ZZ_py(4), Integer)
        assert isinstance(divmod(ZZ(11), ZZ_py(4))[0], Integer)
        assert isinstance(ZZ(7) // ZZ_py(2), Integer)

        assert isinstance(ZZ_py(7) + ZZ(9), Integer)
        assert isinstance(ZZ_py(7) * ZZ(9), Integer)
        assert isinstance(ZZ_py(7) - ZZ(9), Integer)
        assert isinstance(ZZ_py(7) ** ZZ(9), IntegerPython)
        assert isinstance(ZZ_py(11) % ZZ(4), IntegerPython)
        assert isinstance(divmod(ZZ_py(11), ZZ(4))[0], IntegerPython)
        assert isinstance(ZZ_py(7) // ZZ(2), IntegerPython)

    def test_prime(self):
        """Tests is_prime."""
        assert ZZ(104729).is_prime()
        assert not ZZ(104727).is_prime()

    def test_str(self):
        """Tests __str__ methods."""
        assert ZZ(104729).__str__() == "104729"
        assert ZZ_py(104729).__str__() == "104729"
        assert ZZ.__str__() == "The ring of Integers (via flint.fmpz)."
        assert ZZ_py.__str__() == "The ring of Integers (via Python's int)."

    def test_repr(self):
        """Tests __repr__ methods."""
        assert ZZ(104729).__repr__() == "104729"
        assert ZZ_py(104729).__repr__() == "104729"
        assert ZZ.__repr__() == "The ring of Integers (via flint.fmpz)."
        assert ZZ_py.__repr__() == "The ring of Integers (via Python's int)."

    def test_abs(self):
        """Tests __abs__ method."""
        assert abs(ZZ(-5)) == ZZ(5)
        assert abs(ZZ_py(-5)) == ZZ_py(5)

    def test_eq(self):
        """Tests the __eq__ method."""
        assert ZZ(5) == ZZ_py(5)
        assert not ZZ(5) == ZZ(6)
        assert not ZZ_py(5) == ZZ_py(6)
        with pytest.raises(AttributeError):
            ZZ(5) == 5
        with pytest.raises(AttributeError):
            5 == ZZ(5)

    def test_ge(self):
        """Tests the __ge__ method."""
        assert ZZ(5) >= ZZ(5)
        assert ZZ(5) >= ZZ(-13)
        assert ZZ(5) >= ZZ_py(0)
        assert not ZZ(5) >= ZZ(7)
        with pytest.raises(AttributeError):
            ZZ(5) >= 5

    def test_gt(self):
        """Tests the __gt__ method."""
        assert not ZZ(5) > ZZ(5)
        assert ZZ(5) > ZZ(-13)
        assert ZZ(5) > ZZ_py(0)
        assert not ZZ(5) > ZZ(7)
        with pytest.raises(AttributeError):
            ZZ(5) > 5

    def test_le(self):
        """Tests the __le__ method."""
        assert ZZ(5) <= ZZ(5)
        assert not ZZ(5) <= ZZ(-13)
        assert not ZZ(5) <= ZZ_py(0)
        assert ZZ(5) <= ZZ(7)
        with pytest.raises(AttributeError):
            ZZ(5) <= 5

    def test_lt(self):
        """Tests the __lt__ method."""
        assert not ZZ(5) < ZZ(5)
        assert not ZZ(5) < ZZ(-13)
        assert not ZZ(5) < ZZ_py(0)
        assert ZZ(5) < ZZ(7)
        with pytest.raises(AttributeError):
            ZZ(5) < 5

    def test_ne(self):
        """Tests the __ne__ method."""
        assert ZZ(5) != ZZ(6)
        assert ZZ_py(5) != ZZ_py(6)
        assert not ZZ(5) != ZZ_py(5)
        with pytest.raises(AttributeError):
            ZZ(5) != 5


class TestIntegerFunctions:
    """Tests FLINT ZZ functions."""

    def test_gcd(self):
        """Tests the gcd global function."""
        assert gcd(ZZ(15),ZZ(6))==ZZ(3)

    def test_factor(self):
        """Tests the factor function."""
        assert factor(ZZ(12)) == [(ZZ(2),ZZ(2)),(ZZ(3),ZZ(1))]

    def test_factorial(self):
        """Tests the factorial function."""
        assert factorial(ZZ(5))==ZZ(120)
