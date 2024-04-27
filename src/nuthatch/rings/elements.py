"""
ELEMENTS

Base classes for elements in commutative rings. The abstract class stores a
concrete `_data` attribute, on which it is assumed that the standard arithmetic
operations are implemented.
"""


class AbstractRingElement:
    def __init__(self, ring, data):
        self.ring = ring
        self._data = data

    # Arithmetic functions.
    def __add__(self, other):
        """Returns self + other with type that of other."""
        return other.__class__(self._data + other._data)

    def __mul__(self, other):
        """Returns self * other with type that of other."""
        return other.__class__(self._data * other._data)

    def __sub__(self, other):
        """Returns self - other with type that of other."""
        return other.__class__(self._data - other._data)

    def __neg__(self):
        """Returns - self."""
        return self.__class__(-self._data)

    def __truediv__(self, other):
        return NotImplemented

    def __floordiv__(self, other):
        """Returns self // other with type that of self."""
        return self.__class__(self._data // other._data)

    def __divmod__(self, other):
        """Returns (self // other, self % other) with types that of self."""
        d = divmod(self._data, other._data)
        return (self.__class__(d[0]), self.__class__(d[1]))

    def __mod__(self, other):
        """Returns self % other with type that of self."""
        return self.__class__(self._data % other._data)

    def __pow__(self, n):
        """Returns self ** n with type that of self."""
        return self.__class__(self._data**n._data)

    # Comparison functions.
    def __abs__(self):
        return self.__class__(self._data.__abs__())

    def __eq__(self, other):
        return self._data == other._data

    def __ge__(self, other):
        return self._data >= other._data

    def __gt__(self, other):
        return self._data > other._data

    def __le__(self, other):
        return self._data <= other._data

    def __lt__(self, other):
        return self._data < other._data

    def __ne__(self, other):
        return self._data != other._data
