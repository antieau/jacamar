"""
RATIONALS

Base class for QQ.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class Rational(AbstractRingElement):
    _data_class = flint.fmpq

    def __init__(self, n):
        AbstractRingElement.__init__(
            self,
            QQ,
            self._data_class(n),
        )

    def __truediv__(self, other):
        return self.__class__(self._data / other._data)

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self.__str__()


class RationalRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, Rational, exact=True)

    def __call__(self, n):
        return Rational(n)

    def __str__(self):
        return "The ring of Rational Numbers (via flint.fmpq)."

    def __repr__(self):
        return self.__str__()


# Create global version of the integer ring.
QQ = RationalRing()
