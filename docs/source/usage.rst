
Usage
=====

.. _installation:

Installation
------------

To use Nuthatch, first install it using pip:

.. code-block:: console

   (.venv) $ pip install nuthatch

Overview
--------
`Nuthatch` provides support for the collection of mathematical constructions listed:

Basic operations
----------------

**Number Rings**

Number rings behave similarly to those in `Python-FLINT`. 
The real and complex rings in nuthatch are balls, with a number and a radius for precision. 
This radius can be declared and/or will be automatically updated based on any operations performed with the number.


Reals and complexes:


>>> RR(7.3) + RR("2.7 +/- 0.1") == RR(10.05)
True
>>> RR(1.5) + CC(2 + 1j)
CC(3.5 + 1j)
>>> RR(2) * CC(1 + 1j)
CC(2 + 2j)


Integers and rationals:


>>> ZZ(2) + ZZ(4)
ZZ(6)
>>> ZZ(3) ** ZZ(2)
ZZ(9)
>>> QQ((1, 2)) + QQ((1, 4))
QQ((3, 4))
>>> ZZ(1) + QQ((3, 5))
QQ((8, 5))


**Polynomials in `Nuthatch`**

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

>>> poly_ring = PolynomialRing(base_ring=ZZ, ngens=3, prefix='x')
>>> x0, x1, x2 = poly_ring.gens
>>> p = x0 + ZZ(2) * x1 * x2 ** ZZ(2)
>>> p
"x0 + 2x1*x2^2"
>>> poly_ring(p)
"x0 + 2x1*x2^2"
>>> poly_ring(ZZ(4))
"4"

Polynomial operations:

>>> p((ZZ(1), ZZ(1), ZZ(1)))
ZZ(3)
>>> p + p
"2x0 + 4x1*x2^2"
>>> p * ZZ(2)
"2x0 + 4x1*x2^2"
>>> p ** ZZ(2)
"x0^2 + 4x0*x1*x2^2 + 4x1^2*x2^4"


**Matrices in `Nuthatch`:**

A module for dense matrices. We create a single abstract matrix class which
holds a `data` attribute. We assume that this points to a class on which all
matrix operations can be performed and we overload matrices in that way. We
provide a generic implementation of these operations as `_MatrixGeneric`. Special
examples are provided by `FLINT`.

Matrix construction:

>>> Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
[[ZZ(1), ZZ(2)],
[ZZ(3), ZZ(4)]]
>>> random(ZZ_py, 10, 1000, 1000)
1000 x 1000 matrix of random integers on the interval [0, 10)
>>> Matrix(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]], data=_MatrixGenericData(base_ring=PolynomialRing, ncols=2, nrows=2, entries=[[p, p], [p, p]]))
[["x0 + 2x1*x2^2", "x0 + 2x1*x2^2"],
["x0 + 2x1*x2^2", "x0 + 2x1*x2^2"]]

Matrix operations:

>>> m = Matrix(base_ring=ZZ, entries=[[1, 2], [3, 4]])
>>> m + m
[[2, 4],
[6, 8]]
>>> m - m
[[0, 0],
[0, 0]]
>>> m * m
[[5, 10],
[15, 18]]
>>> m * ZZ(2)
[[2, 4],
[6, 8]]

Matrix methods:

>>> m.determinant()
ZZ(-2)
>>> m.transpose()
[[1, 3],
[2, 4]]
>>> m.size()
(2, 2)
>>> m[1, 1]
ZZ(4)

