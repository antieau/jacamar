"""
REALS

Base class for RR.
"""
import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing


class RealNumber(AbstractRingElement):
    data_class = flint.arb
    
    def __init__(self, ring, x=0):
        AbstractRingElement.__init__(
            self,
            ring,
            self.data_class(x)
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
   

class RealNumberPython(AbstractRingElement):
    data_class = float

    def __init__(self, ring, x=0):
        AbstractRingElement.__init__(
            self,
            ring,
            self.data_class(x)
        )

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.__str__()


class RealRing(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, RealNumber, exact=True)

    def __call__(self, *args):
        return RealNumber(self, *args)
    
    def __str__(self):
        return "The ring of real numbers (via flint.arb)."
    
    def __repr__(self):
        return self.__str__()


class RealRingPython(AbstractRing):
    def __init__(self):
        AbstractRing.__init__(self, RealNumberPython, exact=True)

    def __call__(self, *args):
        return RealNumberPython(self, *args)
    
    def __str__(self):
        return "The ring of real numbers (via flint.arb)."
    
    def __repr__(self):
        return self.__str__()
    

RR = RealRing()
RR_py = RealRingPython()
