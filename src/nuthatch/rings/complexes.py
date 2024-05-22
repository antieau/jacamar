"""
COMPLEXES

Base class for CC.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class ComplexNumber(AbstractRingElement):
    """
    Base class for complex numbers, built on AbstractRingElement.
    """

    data_class = flint.acb  # type: ignore

    def __init__(self, *args):
        if len(args) == 3:
            ring, a, b = args
            AbstractRingElement.__init__(self, ring, self.data_class(a, b))

        elif len(args) == 2:
            ring, a = args
            AbstractRingElement.__init__(self, ring, self.data_class(a))

        elif len(args) == 1:
            ring = args[0]
            AbstractRingElement.__init__(self, ring, 0)

        else:
            raise TypeError(
                f"ComplexNumber() takes 1, 2, or 3 arguments, but {len(args)} were given."
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


class ComplexRing(AbstractRing):
    """
    Base class for complex rings, built on AbstractRing.
    """

    def __init__(self):
        AbstractRing.__init__(self, ComplexNumber, exact=True)

    def __call__(self, *args):
        return ComplexNumber(self, *args)

    def __str__(self):
        return "The ring of complex numbers (via flint.acb)."

    def __repr__(self):
        return self.__str__()


# Create a global version of the complex ring.
CC = ComplexRing()


# Flint.arb functions in nuthatch
def sin(x):
    return x.ring(x.data.sin())


def cos(x):
    return x.ring(x.data.cos())


def tan(x):
    return x.ring(x.data.tan())


def asin(x):
    return x.ring(x.data.asin())


def acos(x):
    return x.ring(x.data.acos())


def atan(x):
    return x.ring(x.data.atan())

def csc(x):
    return x.ring(x.data.csc())

def sec(x):
    return x.ring(x.data.sec())

def cot(x):
    return x.ring(x.data.cot())

def acsc(x):
    return x.ring(x.data.acsc())

def asec(x):
    return x.ring(x.data.asec())

def acot(x):
    return x.ring(x.data.acot())

def abs_lower(x):
    return x.ring(x.data.abs_lower())


def abs_upper(x):
    return x.ring(x.data.abs_upper())


def exp(x):
    return x.ring(x.data.exp())


def fac(x):
    return x.ring(x.data.fac())


def log(x):
    return x.ring(x.data.log())

def root(x, n):
    """nth root of x."""
    return x.ring(x.data.root(n))

def sqrt(x):
    return x.ring(x.data.sqrt())

def real(x):
    return x.ring(x.data.real())

def imag(x):
    return x.ring(x.data.imag())
