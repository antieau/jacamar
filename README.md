```
░░█ ▄▀█ █▀▀ ▄▀█ █▀▄▀█ ▄▀█ █▀█ 
█▄█ █▀█ █▄▄ █▀█ █░▀░█ █▀█ █▀▄
```

# Jacamar.

A lightweight framework for computer algebra using `FLINT` and `Python-FLINT`.



# What it is.

- `Jacamar` is designed to provide a minimal layer of abstraction, giving a
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
In particular, all arithmetic operations should fail when mixing `Jacamar`
classes with bare `Python` classes.

Python integers are not automatically coerced. Thus, if `n` is an `int` and `x`
is an element in some ring, `x**n` will raise a `TypeError`.

In general, `__truediv__` is not defined for non-fields. Similarly, powering by
non-integer (Integer or IntegerPy) elements in polynomial rings raises a `TypeError`.
Powering by rational numbers or real numbers might eventually be well-defined
in certain Puiseux or power series rings.



# What it is not.

- `Jacamar` is not meant to be a replacement for `SAGE` or `SymPy`. It has
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
tprint("JACAMAR", font="tarty2")
```


