"""
INTEGERS

Base classes for ZZ, ZZ_word, and ZZ_py.
"""

import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing

class Integer(AbstractRingElement):
    def __init__(self,n):
        self._fmpz = flint.fmpz(n)
        AbstractRingElement.__init__(self,ZZ)

    def __pow__(self,n):
        return self.__class__(self._fmpz**n)

class IntegerPython(AbstractRingElement):
    def __init__(self,n):
        self._int = n
        AbstractRingElement.__init__(self,ZZ_py)

    def __pow__(self,n):
        return self.__class__(self._int**n)

class IntegerRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self,Integer,exact=True)

class IntegerRingPython(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self,IntegerPython,exact=True)



# Create global version of the integer ring.
ZZ = IntegerRing()
ZZ_py = IntegerRingPython()
