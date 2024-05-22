"""
RATIONALS

Base class for QQ.

AUTHORS:
- Benjamin Antieau (2024).
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class Rational(AbstractRingElement):
    data_class = flint.fmpq

    def __init__(self, *args):
        if len(args) == 3:
            ring, p, q = args
            AbstractRingElement.__init__(
                self,
                ring,
                self.data_class(p, q),
            )

        elif len(args) == 2:
            ring, p = args
            AbstractRingElement.__init__(
                self,
                ring,
                p,
            )

        elif len(args) == 1:
            ring = args[0]
            AbstractRingElement.__init__(
                self,
                ring,
                0,
            )

        else:
            raise TypeError(
                f"Rational() takes 1, 2, or 3 arguments, but {len(args)} were given."
            )

    def __truediv__(self, other):
        """Returns self / other with type that of self."""
        return self.__class__(self.ring, self.data / other.data)

    def __itruediv__(self, other):
        """Modifies self in place by self / other."""
        self.data /= other.data
        return self

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()


class RationalRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, Rational, exact=True)

    def __call__(self, *args):
        return Rational(self, *args)

    def __str__(self):
        return "The ring of Rational numbers (via flint.fmpq)."

    def __repr__(self):
        return self.__str__()


# Create global version of the rational ring.
QQ = RationalRing()


def p(x):
    """Returns numerator of x as ZZ."""
    return x.ring(x.data.p())

def q(x):
    """Returns denominator of x as ZZ."""
    return x.ring(x.data.q())