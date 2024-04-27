"""
INTEGERS

Base classes for ZZ and ZZ_py.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class Integer(AbstractRingElement):
    _data_class = flint.fmpz

    def __init__(self, n):
        AbstractRingElement.__init__(
            self,
            ZZ,
            self._data_class(n),
        )

    def is_prime(self):
        return bool(self._data.is_prime())

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self.__str__()


class IntegerPython(AbstractRingElement):
    _data_class = int

    def __init__(self, n):
        AbstractRingElement.__init__(
            self,
            ZZ_py,
            self._data_class(n),
        )

    def __str__(self):
        return self._data.__str__()

    def __repr__(self):
        return self.__str__()


class IntegerRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, Integer, exact=True)

    def __str__(self):
        return "The ring of Integers (via flint.fmpz)."

    def __repr__(self):
        return self.__str__()


class IntegerRingPython(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, IntegerPython, exact=True)

    def __str__(self):
        return "The ring of Integers (via Python's int)."

    def __repr__(self):
        return self.__str__()


# Create global version of the integer ring.
ZZ = IntegerRing()
ZZ_py = IntegerRingPython()
