"""
# JACAMAR

A tiny Python computer algebra system.

JACAMAR provides three basic sets of functionality:
- a high-level interface to FLINT, via the PythonFLINT project;
- generic code for matrices;
- generic code for multivariable Polynomials, LaurentPolynomials, and PowerSeries.
Where possible, we follow the conventions of SAGE for the creation of objects.



# Philosophy

JACAMAR favors explicit, user-imposed coercion when it makes sense. In
particular, Python `int` instances are not mathematical objects and must typically be
coerced mannualy into JACAMAR objects.

Several base commutative rings are provided. At the moment, these are ZZ, ZZ_py, and QQ.
The rings ZZ and QQ wrap FLINT classes `fmpz` and `fmpq`, while ZZ_py wraps Python's `int` class.
We will be adding RR and CC in the coming days and then later integers modulo
N, finite fields, and q-adics.

New commutative rings are constructed by building PolynomialRings, LaurentPolynomialRings,
or PowerSeries rings. These rings each refer to a class of element objects. For
example, `Polynomial` is the generic element class of a PolynomialRing. The class
`Polynomial` itself is a wrapper for a lower-level _Polynomial class.

However, for those rings provided by FLINT, the data class will be the
corresponding FLINT polynomial class.

The goal is to provide a very thin, familiar interface which the user can
easily abstract and customize for their own purposes.
"""

from jacamar.matrices.matrices import *

from jacamar.rings.elements import *
from jacamar.rings.rings import *
from jacamar.rings.morphisms import *
from jacamar.rings.integers import *
from jacamar.rings.intmod import *
from jacamar.rings.polynomials import *
from jacamar.rings.rationals import *
from jacamar.rings.series import *
from jacamar.rings.reals import *
from jacamar.rings.complexes import *

__version__ = "0.1.0"
