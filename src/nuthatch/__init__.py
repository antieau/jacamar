"""
# NUTHATCH

A tiny Python computer algebra system.

NUTHATCH provides three basic sets of functionality:
- a high-level interface to FLINT, via the PythonFLINT project;
- generic code for matrices;
- generic code for multivariable Polynomials, LaurentPolynomials, and PowerSeries.
Where possible, we follow the conventions of SAGE for the creation of objects.



# Philosophy

NUTHATCH favors explicit, user-imposed coercion when it makes sense. In
particular, Python `int` instances are not mathematical objects and must typically be
coerced mannualy into NUTHATCH objects.

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


from nuthatch.matrices.matrices import *

from nuthatch.rings.elements import *
from nuthatch.rings.rings import *
from nuthatch.rings.morphisms import *
from nuthatch.rings.integers import *
from nuthatch.rings.polynomials import *
from nuthatch.rings.rationals import *
from nuthatch.rings.series import *
