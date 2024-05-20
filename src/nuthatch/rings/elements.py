"""
ELEMENTS

Base classes for elements in commutative rings. The abstract class stores a
concrete `data` attribute, on which it is assumed that the standard arithmetic
operations are implemented.

AUTHORS:
- Benjamin Antieau (2024).
"""


class AbstractRingElement:
    def __init__(self, ring, data):
        self.ring = ring
        # TODO: does this have unexpected copying behavior? It seems like if we
        # change .data in one place it will be changed elsewhere.
        if isinstance(data, self.ring.element_class):
            self.data = data.data
        elif isinstance(data, self.ring.element_class.data_class):
            self.data = data
        else:
            self.data = self.ring.element_class.data_class(data)

    # Arithmetic functions.
    def __add__(self, other):
        """Returns self + other with type that of other."""
        return other.__class__(other.ring, self.data + other.data)

    def __mul__(self, other):
        """Returns self * other with type that of other."""
        return other.__class__(other.ring, self.data * other.data)

    def __sub__(self, other):
        """Returns self - other with type that of other."""
        return other.__class__(other.ring, self.data - other.data)

    def __neg__(self):
        """Returns - self."""
        return self.__class__(self.ring, -self.data)

    def __truediv__(self, other):
        """Returns self / other with type that of self."""
        return NotImplemented

    def __floordiv__(self, other):
        """Returns self // other with type that of self."""
        return self.__class__(self.ring, self.data // other.data)

    def __divmod__(self, other):
        """Returns (self // other, self % other) with types that of self."""
        d = divmod(self.data, other.data)
        return (self.__class__(self.ring, d[0]), self.__class__(self.ring, d[1]))

    def __mod__(self, other):
        """Returns self % other with type that of self."""
        return self.__class__(self.ring, self.data % other.data)

    def __pow__(self, n):
        """Returns self ** n with type that of self."""
        return self.__class__(self.ring, self.data**n.data)

    # Comparison functions.
    def __abs__(self):
        return self.__class__(self.ring, self.data.__abs__())

    def __eq__(self, other):
        return self.data == other.data

    def __ge__(self, other):
        return self.data >= other.data

    def __gt__(self, other):
        return self.data > other.data

    def __le__(self, other):
        return self.data <= other.data

    def __lt__(self, other):
        return self.data < other.data

    def __ne__(self, other):
        return self.data != other.data
