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

Matrices in `Nuthatch` can be constructed in different ways, depending on the type of object in the matrix.
If the object is a ZZ, QQ, RR, or CC, the matrix will use the `Python-FLINT` mat framework.
If the object is a python wrapper (ZZ_py, QQ_py, RR_py, or CC_py), the matrix will use the `numpy.array` framework.
If the object is something else (polynomial, series, etc...), the matrix will be generic and use `Nuthatch` methods.

Matrix construction:

```
> Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
[[ZZ(1), ZZ(2)],
[ZZ(3), ZZ(4)]]
> random(ZZ_py, 10, 1000, 1000)
1000 x 1000 matrix of random integers on the interval [0, 10)
> poly_ring = PolynomialRing(base_ring=ZZ, ngens=3, prefix='x')
> x0, x1, x2 = poly_ring.gens
> p = x0 + ZZ(2) * x1 * x2 ** ZZ(2)
> Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]]))
[[x0 + 2x1*x2^2, x0 + 2x1*x2^2],
[x0 + 2x1*x2^2, x0 + 2x1*x2^2]]
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
