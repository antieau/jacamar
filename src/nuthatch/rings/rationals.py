"""
RATIONALS

Base class for QQ.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class Rational(AbstractRingElement):
    _data_class = flint.fmpq

    def __init__(self, *args):
        if not args:
            AbstractRingElement.__init__(
                self,
                QQ,
                0,
            )
        elif len(args) == 2:
            p, q = args
            AbstractRingElement.__init__(
                self,
                QQ,
                self._data_class(p, q),
            )

        elif len(args) == 1:
            p = args[0]
            AbstractRingElement.__init__(
                self,
                QQ,
                p,
            )

        else:
            raise TypeError(
                f"Rational() takes at most 2 arguments, but {len(args)} were given."
            )

    def __truediv__(self, other):
        """Returns self / other with type that of self."""
        return self.__class__(self._data / other._data)

    def __itruediv__(self, other):
        """Modifies self in place by self / other."""
        self._data /= other._data
        return self

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self.__str__()


class RationalRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, Rational, exact=True)

    def __call__(self, *args):
        return Rational(*args)

    def __str__(self):
        return "The ring of Rational numbers (via flint.fmpq)."

    def __repr__(self):
        return self.__str__()


# Create global version of the integer ring.
QQ = RationalRing()
