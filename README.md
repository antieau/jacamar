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



# Release plan.

- `v0.1.0` - AbstractRings, ZZ, `ZZ_py`, `ZZ_word`, GF(p), PolynomialRings, Matrices
- `v0.2.0` - LaurentPolynomialRings, PowerSeriesRings
