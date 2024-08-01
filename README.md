```
______________________________________________________________________________
    _     _    _     _  ______    _     _    __    ______      __      _     _
    /|   /     /    /     /       /    /     / |     /       /    )    /    / 
---/-| -/-----/----/-----/-------/___ /-----/__|----/-------/---------/___ /--
  /  | /     /    /     /       /    /     /   |   /       /         /    /   
_/___|/_____(____/_____/_______/____/_____/____|__/_______(____/____/____/____
```

# Nuthatch.

A lightweight framework for computer algebra using `FLINT` and `Python-FLINT`
at its core.



# What it is.

- `Nuthatch` is designed to provide a minimal layer of abstraction, giving a
`SAGE`-like feel to `Python-FLINT`. The goal is for parity to `Python-FLINT` in
terms of speed.
- Its audience is mathematicians and computer scientists who need to use
algebraic objects (commutative rings, polynomials, chain complexes, and categories) to
easily build new algebraic objects.
- A `Python` package, which plays well with the ecosystem of packages.
Eventually, this will include wheel builds and `pip install` capabilities.
- Built for parallelizability?
- Mathematical.
- Fundamental principle: no automatic coercions.



# Coercion model.

We more or less inherit some basic coercions from `Python-FLINT`. Where it
makes sense, binary operations like `__add__(self,other)` or
`__mul__(self,other)` return elements of type `other.__class__`. Thus, the
constructor for this class must be able to take the output of the underlying
data `self._data + other._data` and construct a new element with it. In
practice, this corresponds to pretending like all such operations correspond to
a left action of the ring of `self` on the ring or module of `other`.

We try for this to fail whenever something non-obvious is happening.
In particular, all arithmetic operations should fail when mixing `Nuthatch`
classes with bare `Python` classes.

Python integers are not automatically coerced. Thus, if `n` is an `int` and `x`
is an element in some ring, `x**n` will raise a `TypeError`.

In general, `__truediv__` is not defined for non-fields. Similarly, powering by
non-integer (Integer or IntegerPy) elements in polynomial rings raises a `TypeError`.
Powering by rational numbers or real numbers might eventually be well-defined
in certain Puiseux or power series rings.


# What it is not.

- `Nuthatch` is not meant to be a replacement for `SAGE` or `SymPy`. It has
a tiny fraction of the capabilities of `SAGE`, a tiny fraction of the lifetime
of tests that `SAGE` has passed, and none of the printing and display options.



# Remarks.

- We do not provide a universal element-parent approach as in `SAGE`. Rather,
instances of rings point to element classes via `self.element_class`
and instances of ring elements have `self.ring`. We use this as well for objects like polynomial rings
or matrices.
- We prefer to access immutable characteristics of elements, like `self.ring`
for an element, or `self.codomain` for some kind of morphism, as attributes
(or @properties) rather than as functions. 
- We do not automatically coerce numbers into Integers. A number such as `7` in
the Python prompt will always still produce an `int`, so it can be used as an
index in lists, etc.


# Tools.

- pygount
- ruff
- black
- pylint
- coverage
- python environments


# Style.


```
from art import *
tprint("NUTHATCH", font="bigchief")
```

# Example Use.

Number rings behave similarly to those in `Python-FLINT`. 
The real and complex rings in nuthatch are balls, with a number and a radius for precision. 
This radius can be declared and/or will be automatically updated based on any operations performed with the number.

Number ring basic usage:

Reals and complexes:

```
> RR(7.3) + RR("2.7 +/- 0.1") == RR(10.05)
True
> RR(1.5) + CC(2 + 1j)
CC(3.5 + 1j)
> RR(2) * CC(1 + 1j)
CC(2 + 2j)
```

Integers and rationals:

```
> ZZ(2) + ZZ(4)
ZZ(6)
> ZZ(3) ** ZZ(2)
ZZ(9)
> QQ((1, 2)) + QQ((1, 4))
QQ((3, 4))
> ZZ(1) + QQ((3, 5))
QQ((8, 5))
```

Polynomials in `Nuthatch`

Classes for generic polynomial rings.

Currently, there are two supported exponent types: PackedMonomialData and
SparseMonomialData.

The former are packed into a Python `int`. The largest power supported is given
by the module-level constant PACKING_BOUND.

The latter are designed based on the `ETuple` class from `SAGE`. In particular,
they are sparse. In a polynomial ring on x0, x1, x2, x3, a monomial
x0*x2^3*x3^7 would be encoded as (0,1,2,3,3,7), where the even indices indicate
the present variables and the odd indices indicate their powers. The empty
tuple () represents the monomial of 1.

Polynomial construction:

```
> poly_ring = PolynomialRing(base_ring=ZZ, ngens=3, prefix='x')
> x0, x1, x2 = poly_ring.gens
> p = x0 + ZZ(2) * x1 * x2 ** ZZ(2)
> p
"x0 + 2x1*x2^2"
> poly_ring(p)
"x0 + 2x1*x2^2"
> poly_ring(ZZ(4))
"4"
```

Polynomial operations:

```
> p((ZZ(1), ZZ(1), ZZ(1)))
ZZ(3)
> p + p
"2x0 + 4x1*x2^2"
> p * ZZ(2)
"2x0 + 4x1*x2^2"
> p ** ZZ(2)
"x0^2 + 4x0*x1*x2^2 + 4x1^2*x2^4"
```

Matrices in `Nuthatch`:

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.

Matrix construction:

```
> Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
[[ZZ(1), ZZ(2)],
[ZZ(3), ZZ(4)]]
> random(ZZ_py, 10, 1000, 1000)
1000 x 1000 matrix of random integers on the interval [0, 10)
> Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]]))
[["x0 + 2x1*x2^2", "x0 + 2x1*x2^2"],
["x0 + 2x1*x2^2", "x0 + 2x1*x2^2"]]
```

Matrix operations:

```
> m = Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
> m + m
[[2, 4],
[6, 8]]
> m - m
[[0, 0],
[0, 0]]
> m * m
[[5, 10],
[15, 18]]
> m * ZZ(2)
[[2, 4],
[6, 8]]
```

Matrix methods:

```
> m.determinant()
ZZ(-2)
> m.transpose()
[[1, 3],
[2, 4]]
> m.size()
(2, 2)
> m[1, 1]
ZZ(4)
```
