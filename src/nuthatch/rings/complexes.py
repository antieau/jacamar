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



"""
Select functions of CC's based on flint.acb attributes. (https://fredrikj.net/python-flint/acb.html)
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

def real(x):
    """
    Returns the real component of x.
    """
    return x.ring(x.data.real())


def imag(x):
    """
    Returns the imaginary component of x.
    """
    return x.ring(x.data.imag())
