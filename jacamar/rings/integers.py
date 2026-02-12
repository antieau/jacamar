"""
INTEGERS

Base classes for ZZ and ZZ_py.

AUTHORS:
- Benjamin Antieau (2024).
"""

import flint
from jacamar.rings.elements import AbstractRingElement
from jacamar.rings.rings import AbstractRing


class Integer(AbstractRingElement):
    data_class = flint.fmpz

    def __init__(self, ring, n):
        AbstractRingElement.__init__(
            self,
            ring,
            n,
        )

    def is_prime(self):
        return bool(self.data.is_prime())

    def is_unit(self):
        return self.data == flint.fmpz(1) or self.data == flint.fmpz(-1)

    def inverse(self):
        if self.data == flint.fmpz(1):
            return self
        elif self.data == flint.fmpz(-1):
            return self
        else:
            raise ValueError("Is not invertible.")

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()


class IntegerPython(AbstractRingElement):
    data_class = int

    def __init__(self, ring, n):
        AbstractRingElement.__init__(
            self,
            ring,
            n,
        )

    def is_unit(self):
        return self.data == -1 or self.data == 1

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()


class IntegerRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, Integer, exact=True)
        self.one = self(1)
        self.zero = self(0)

    def __str__(self):
        return "The ring of Integers (via flint.fmpz)."

    def __repr__(self):
        return self.__str__()

    def combinations(self, n, k):
        return self.element_class(ring=self,n=flint.fmpz.bin_uiui(n,k))

    binomial = combinations


class IntegerRingPython(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, IntegerPython, exact=True)
        self.one = self(1)
        self.zero = self(0)

    def __str__(self):
        return "The ring of Integers (via Python's int)."

    def __repr__(self):
        return self.__str__()


# Create global version of the integer rings.
ZZ = IntegerRing()
ZZ_py = IntegerRingPython()


"""
Select functions of ZZ's based on flint.fmpz attributes. (https://fredrikj.net/python-flint/fmpz.html)
"""


def factorial(x):
    """
    Integer factorial of x.
    """
    return x.ring(x.data_class.fac_ui(x.data))


def gcd(x, y):
    """
    Greatest common factor of x and y.
    """
    return x.ring(x.data.gcd(y.data))


def factor(x):
    """
    Factors of x in polynomial form.
    """
    ret = []
    factors = x.data.factor()
    for f in factors:
        ret.append((x.ring(f[0]), x.ring(f[1])))
    return ret
