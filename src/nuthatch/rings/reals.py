"""
REALS

Base class for RR.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement  # type: ignore
from nuthatch.rings.rings import AbstractRing  # type: ignore


class RealNumber(AbstractRingElement):
    """
    Base class for real numbers, built on AbstractRingElement.
    """

    data_class = flint.arb  # type: ignore

    def __init__(self, ring, x=0):
        AbstractRingElement.__init__(self, ring, self.data_class(x))

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


class RealNumberPython(AbstractRingElement):
    """
    Base class for real numbers, built on AbstractRingElement. Interfaces with python.
    """

    data_class = float

    def __init__(self, ring, x=0):
        AbstractRingElement.__init__(self, ring, self.data_class(x))

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()


class RealRing(AbstractRing):
    """
    Base class for real rings, built on AbstractRing.
    """

    def __init__(self):
        AbstractRing.__init__(self, RealNumber, exact=True)

    def __call__(self, *args):
        return RealNumber(self, *args)

    def __str__(self):
        return "The ring of real numbers (via flint.arb)."

    def __repr__(self):
        return self.__str__()


class RealRingPython(AbstractRing):
    """
    Base class for real rings, built on AbstractRing. Interfaces with python.
    """

    def __init__(self):
        AbstractRing.__init__(self, RealNumberPython, exact=True)

    def __call__(self, *args):
        return RealNumberPython(self, *args)

    def __str__(self):
        return "The ring of real numbers (via flint.arb)."

    def __repr__(self):
        return self.__str__()


# Create a global version of the real ring.
RR = RealRing()
RR_py = RealRingPython()



"""
Select functions of RR's based on flint.arb attributes. (https://fredrikj.net/python-flint/arb.html)
"""


def sin(x):
    """
    Sine of x.
    """
    return x.ring(x.data.sin())


def cos(x):
    """
    Cosine of x.
    """
    return x.ring(x.data.cos())


def tan(x):
    """
    Tangent of x.
    """
    return x.ring(x.data.tan())


def asin(x):
    """
    Inverse sine of x.
    """
    return x.ring(x.data.asin())


def acos(x):
    """
    Inverse cosine of x.
    """
    return x.ring(x.data.acos())


def atan(x):
    """
    Inverse tanget of x.
    """
    return x.ring(x.data.atan())


def csc(x):
    """
    Cosecant of x.
    """
    return x.ring(x.data.csc())


def sec(x):
    """
    Secant of x.
    """
    return x.ring(x.data.sec())


def cot(x):
    """
    Cotanget of x.
    """
    return x.ring(x.data.cot())


def acsc(x):
    """
    Inverse cosecant of x.
    """
    return x.ring(x.data.acsc())


def asec(x):
    """
    Inverse secant of x.
    """
    return x.ring(x.data.asec())


def acot(x):
    """
    Inverse cotanget of x.
    """
    return x.ring(x.data.acot())


def abs_lower(x):
    """
    The absolute value of the lower bound of x.
    """
    return x.ring(x.data.abs_lower())


def abs_upper(x):
    """
    The absolute value of the upper bound of x.
    """
    return x.ring(x.data.abs_upper())


def exp(x):
    """
    The constat e raised to the power of x.
    """
    return x.ring(x.data.exp())


def fac(x):
    """
    Factorial, using the gamma function of x.
    """
    return x.ring(x.data.fac())


def log(x):
    """
    Natural logarithm of x.
    """
    return x.ring(x.data.log())


def root(x, n):
    """nth root of x."""
    return x.ring(x.data.root(n))


def sqrt(x):
    """
    Square root of x.
    """
    return x.ring(x.data.sqrt())
