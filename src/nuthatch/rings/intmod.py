"""
INTEGERS MODULO N

Base class for ZN.

AUTHORS:
- Joseph Grantham (2024).
"""

from flint import nmod # pylint: disable=no-name-in-module
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing
from nuthatch.rings.integers import ZZ

class IntegerModN(AbstractRingElement):
    """Element representing an integer modulo N"""
    data_class = nmod

    def __init__(self, ring, data):
        if isinstance(data, AbstractRingElement):
            data = data.data
        if isinstance(ring.modulus, AbstractRingElement):
            data = nmod(data, ring.modulus.data)
        AbstractRingElement.__init__(
            self,
            ring,
            data,
        )

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()

    def __truediv__(self, other):
        return self.__class__(self.ring, self.data / other.data)

    def __floordiv__(self, other):
        return NotImplemented

    def __divmod__(self, other):
        return NotImplemented

    def __mod__(self, other):
        return NotImplemented

    def __abs__(self):
        raise AttributeError

    def __ge__(self, other):
        return NotImplemented

    def __gt__(self, other):
        return NotImplemented

    def __le__(self, other):
        return NotImplemented

    def __lt__(self, other):
        return NotImplemented


class IntegerModNRing(AbstractRing):
    """Ring of integers modulo N"""
    def __init__(self, modulus):
        AbstractRing.__init__(self, IntegerModN, exact=True)
        if isinstance(modulus, AbstractRingElement):
            modulus = ZZ(modulus.data)
        else:
            modulus = ZZ(modulus)
        self.modulus = modulus
        self.one = self(1)
        self.zero = self(0)

    def __str__(self):
        return "The ring of Integers modulo " + str(self.modulus) + "."

    def __repr__(self):
        return self.__str__()
