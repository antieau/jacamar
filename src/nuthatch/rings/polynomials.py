"""
POLYNOMIALS

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

AUTHORS:
- Benjamin Antieau (2024): initial version.
"""

import functools
import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.morphisms import AbstractRingMorphism

# The following constant controls the maximum allowed weight of a power x^n in a
# monomial.
PACKING_BOUND = 2**16


class MonomialData:
    def __init__(self):
        pass

    @classmethod
    def from_packed_integer(cls, n):
        return NotImplemented

    @classmethod
    def from_sparse_tuple(cls, t):
        return NotImplemented

    @classmethod
    def from_tuple(cls, t):
        return NotImplemented



class PackedMonomialData(MonomialData):
    """
    A class for packing monomial data into a single integer.
    """
    def __init__(self, n):
        self.weight = n
        MonomialData.__init__(self)

    @classmethod
    def from_packed_integer(cls, n):
        return cls(n)

    @classmethod
    def from_sparse_tuple(cls, t):
        x = 0
        for i in range(len(t)//2):
            x+=t[2*i+1]*(PACKING_BOUND**t[2*i])
        return cls(x)

    @classmethod
    def from_tuple(cls, t):
        x = 0
        for i,pwr in enumerate(t):
            x += pwr * (PACKING_BOUND**i)
        return cls(x)

    @functools.cached_property
    def degrees(self):
        """Returns the sparse tuple representation of self."""
        count = 0
        x = []
        n = self.weight
        while n > 0:
            q, r = n.__divmod__(PACKING_BOUND)
            if r > 0:
                x.extend((count,r))
            n = q
            count += 1
        return tuple(x)

    def __hash__(self):
        return self.weight.__hash__()

    def __eq__(self, other):
        return self.weight == other.weight

    def __mul__(self, other):
        """
        Multiplication of monomials.

        WARNING: this could result in nonsense if there is an overflow
        encountered. Ensure that no variable is powered to more than the
        module-level constant PACKING_BOUND.
        """
        return other.__class__(self.weight+other.weight)

    def __str__(self):
        return str(self.weight)

    def __repr__(self):
        return self.__str__()


class SparseMonomialData(MonomialData):
    """
    The class of a monomial in a PolynomialData.

    INPUT:
    - `t`     -- a tuple of Python integers
    """

    def __init__(self, *t):
        if len(t) == 1 & isinstance(t,tuple):
            self.degrees = t[0]
        else:
            self.degrees = tuple(t)

    @classmethod
    def from_packed_integer(cls, n):
        count = 0
        x = []
        while n > 0:
            q, r = n.__divmod__(PACKING_BOUND)
            if r > 0:
                x.extend((count,r))
            n = q
            count += 1
        return cls(tuple(x))

    @classmethod
    def from_sparse_tuple(cls, t):
        return cls(t)

    @classmethod
    def from_tuple(cls, t):
        x = []
        for i,pwr in enumerate(t):
            if pwr != 0:
                x.extend((i,pwr))
        return cls(tuple(x))

    def __hash__(self):
        return hash(self.degrees)

    def __str__(self):
        return str(self.degrees)

    def __repr__(self):
        return str(self)

    def __mul__(self, other):
        new_list = []
        self_index = 0
        other_index = 0
        self_len = len(self.degrees)
        other_len = len(other.degrees)
        while self_index < self_len and other_index < other_len:
            if self.degrees[self_index] < other.degrees[other_index]:
                new_list.extend(self.degrees[self_index : self_index + 2])
                self_index += 2
            elif self.degrees[self_index] == other.degrees[other_index]:
                x = self.degrees[self_index + 1] + other.degrees[other_index + 1]
                if x != 0:
                    new_list.append(self.degrees[self_index])
                    new_list.append(x)
                self_index += 2
                other_index += 2
            else:
                new_list.extend(other.degrees[other_index : other_index + 2])
                other_index += 2

        # Now that we've reached the end of at least one list, we add
        # everything from the other list.
        while self_index < self_len:
            new_list.extend(self.degrees[self_index : self_index + 2])
            self_index += 2

        while other_index < other_len:
            new_list.extend(other.degrees[other_index : other_index + 2])
            other_index += 2

        return other.__class__(tuple(new_list))

    def __eq__(self, other):
        return self.degrees == other.degrees
        


class PolynomialData:
    """
    The underlying data class of polynomials.
    The `data` input should be a dictionary `{m:c}` where `m` is a `MonomialData`
    and `c` is an element of the data_class of the element_class of the `base_ring`.
    Each operation returns a PolynomialData class with no terms with
    coefficient equal to zero if the inputs have this property.

    In practice, the methods in this class work only under the assumption that
    the keys are of a fixed type of MonomialData. It is not designed so that
    the types can be mixed.
    """

    def __init__(self, base_ring, monomial_dictionary):
        self.base_ring = base_ring
        self.monomial_dictionary = {}
        for m,c in monomial_dictionary.items():
            if c != self.base_ring.zero.data:
                self.monomial_dictionary[m] = c

    def is_zero(self):
        if len(self.monomial_dictionary) == 0:
            return True
        for _, c in self.monomial_dictionary.items():
            if c != self.base_ring.zero.data:
                return False
        return True

    def term_data(self):
        """Returns the iterator over the items {m:c} in self.monomial_dictionary."""
        return self.monomial_dictionary.items()

    def __str__(self):
        return str(self.monomial_dictionary)

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        if len(self.monomial_dictionary) < len(other.monomial_dictionary):
            new_dict = other.monomial_dictionary.copy()
            other_dict = self.monomial_dictionary

        else:
            new_dict = self.monomial_dictionary.copy()
            other_dict = other.monomial_dictionary

        for m, c in other_dict.items():
            try:
                new_dict[m] += c
            except KeyError:
                new_dict[m] = c
        return other.__class__(other.base_ring, new_dict)

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        new_dict = self.monomial_dictionary.copy()
        for m, c in new_dict.items():
            new_dict[m] = -c
        return self.__class__(self.base_ring, new_dict)

    def __mul__(self, other):
        if self.base_ring != other.base_ring:
            return NotImplemented
        new_dict = {}
        for m, c in self.monomial_dictionary.items():
            for n, d in other.monomial_dictionary.items():
                k = m * n
                e = c * d
                try:
                    new_dict[k] += e
                except KeyError:
                    new_dict[k] = e
        return other.__class__(other.base_ring, new_dict)

    def __rmul__(self, other):
        """
        We just try to multiply the coefficients of self with other.
        """
        new_dict = {}
        for m, c in self.monomial_dictionary.items():
            new_dict[m] = c * other
        return self.__class__(self.base_ring, new_dict)

    def __eq__(self, other):
        return (self.base_ring == other.base_ring) & (
            self.monomial_dictionary == other.monomial_dictionary
        )

    def __pow__(self, n):
        # TODO: this low level power function should not really be checking for
        # types.
        if isinstance(n, (flint.fmpz, int)):
            if n < 0:
                raise TypeError("Cannot power a polynomial by a negative integer.")
            elif n == 0:
                return self.ring.one.data

            apow = self
            while not (n & 1):
                # While even
                apow *= apow
                # Bitshift
                n >>= 1

            # Now multiply together the correct factors a^(2^i)
            res = apow
            n >>= 1
            while n:
                apow *= apow
                if n & 1:
                    res = apow * res
                n >>= 1
            return res
        else:
            raise TypeError(f"Cannot power polynomial by {n}.")

    def __call__(self, arg):
        """The PolynomialData class can be called on a monomial."""
        try:
            return self.monomial_dictionary[arg]
        except KeyError:
            return self.base_ring.zero.data

    def evaluate(self, args):
        """Returns f(a_1,...,a_n) where arguments=[a_1,...,a_n]."""
        x = self.base_ring.zero.data
        for key, value in self.monomial_dictionary.items():
            monomial_part = self.base_ring.one.data
            for i in range(len(key.degrees) // 2):
                monomial_part *= args[key.degrees[2*i]]**key.degrees[2*i+1]
            x += value * monomial_part
        return x


class Polynomial(AbstractRingElement):
    """
    The ring of coefficients is `self.base_ring` while the ambient polynomial
    ring is `self.ring`.
    """

    data_class = PolynomialData

    def __init__(self, ring, data):
        self.base_ring = ring.base_ring
        AbstractRingElement.__init__(
            self,
            ring,
            data,
        )

    def __call__(self, *args):
        """
        Evaluate self at x. If x is an instance of MonomialData, returns the
        coefficient. If x is a tuple of elements, then evaluate on these.
        """
        if isinstance(args[0], MonomialData):
            return self.base_ring(self.data(args[0]))

        if len(args) == self.ring.ngens:
            try:
                return self.base_ring(self.data.evaluate([x.data for x in args]))
            except AttributeError:
                pass

            try:
                m = self.ring._monomial_class.from_tuple(tuple(args))
                return self.base_ring(self.data(m))
            except ValueError:
                pass

        return NotImplemented

    def __str__(self):
        if len(self.data.monomial_dictionary) == 0:
            return "0"
        keys = list(self.data.monomial_dictionary.keys())

        # Parse the first key
        key = keys[0]
        value = self.data.monomial_dictionary[key]
        try:
            bl = value < self.base_ring.zero.data
        except TypeError:
            bl = bool(0)
        if bl:
            token = str(-self.base_ring(value))
        else:
            token = str(self.base_ring(value))
        if token.count(" ") > 0:
            if bl:
                term = "- (" + token + ")"
            else:
                term = "(" + token + ")"
        else:
            if bl:
                term = "- " + token
            else:
                term = token
        
        for j in range(len(key.degrees) // 2):
            term += "*"
            term += self.ring._prefix
            term += str(key.degrees[2 * j])
            term += f"^{key.degrees[2*j+1]}"
        return_string = term

        # Parse the latter keys
        for i in range(1, len(keys)):
            key = keys[i]
            value = self.data.monomial_dictionary[key]
            # Parse the key
            try:
                bl = value < self.base_ring.zero.data
            except TypeError:
                bl = bool(0)
            if bl:
                token = str(-self.base_ring(value))
            else:
                token = str(self.base_ring(value))
            if token.count(" ") > 0:
                if bl:
                    term = " - (" + token + ")"
                else:
                    term = " + (" + token + ")"
            else:
                if bl:
                    term = " - " + token
                else:
                    term = " + " + token
            for j in range(len(key.degrees) // 2):
                term += "*"
                term += self.ring._prefix
                term += str(key.degrees[2 * j])
                term += f"^{key.degrees[2*j+1]}"
            return_string += term
        return return_string

    def __repr__(self):
        return str(self)


class PolynomialRing(AbstractRing):
    """
    Class for generic (multivariable) polynomial rings.

    INPUT:
    - `base_ring`   -- the FLINT ring from which the coefficients are taken
    - `ngens`       -- the numberg of generators
    - `prefix`      -- a string such as 'x' from which the generators x0,x1,... are named
    - `weights`     -- a list of `ngens` weights; if `None`, these default to 1; they are Python integers
    """

    element_class = Polynomial

    def __init__(
        self,
        *,
        base_ring,
        ngens,
        prefix,
        weights=None,
        packed=True,
    ):
        self.base_ring = base_ring
        self.ngens = ngens
        if packed:
            self._monomial_class = PackedMonomialData
        else:
            self._monomial_class = SparseMonomialData
        self.packed = packed
        self._prefix = prefix
        self._names = [prefix + str(i) for i in range(ngens)]
        if weights is None:
            self.weights = [1 for i in range(ngens)]
        else:
            if len(weights) != ngens:
                raise ValueError(
                    f"The number {len(weights)} is not equal to the number of generators {ngens}."
                )
            self.weights = weights

        self.gens = []
        for i in range(self.ngens):
            self.gens.append(
                Polynomial(
                    self,
                    PolynomialData(
                        self.base_ring, {self._monomial_class.from_sparse_tuple((i, 1)): self.base_ring.one.data}
                    ),
                )
            )

        AbstractRing.__init__(self, Polynomial, exact=base_ring.exact)

        # Initialize one and zero since these are used so often.
        self.one = self(1)
        self.zero = self(0)


    def __str__(self):
        return f"Ring of polynomials in {self.ngens} variables {self._names} with weights {self.weights} over {self.base_ring}"

    def __repr__(self):
        return self.__str__()

    def __call__(self, data):
        if isinstance(data,int):
            return Polynomial(
                self,
                PolynomialData(
                    self.base_ring,
                    {self._monomial_class.from_sparse_tuple(()): self.base_ring(data).data},
                ),
            )
        if isinstance(data, dict):
            new_dict = {}
            for key, value in data.items():
                new_dict[self._monomial_class.from_sparse_tuple(key)] = value.data
            return Polynomial(self, PolynomialData(self.base_ring, new_dict))
        if isinstance(data, self.element_class):
            if data.ring == self:
            # Create a new element with the same data.
                return self.element_class(self, data.data)
        if isinstance(data, self.element_class.data_class):
            if data.base_ring == self.base_ring:
                return self.element_class(self,data)
        if isinstance(data, self.base_ring.element_class):
            if data.ring == self.base_ring:
                return Polynomial(
                    self,
                    PolynomialData(
                        self.base_ring,
                        {self._monomial_class.from_sparse_tuple(()): self.base_ring(data).data},
                    ),
                )
        raise TypeError(f"No known constructor for input data of type {type(data)}.")


class PolynomialRingMorphism(AbstractRingMorphism):
    """
    Homomorphisms of polynomial rings.
    """

    def __init__(
        self, *, domain, codomain, coefficient_morphism, action_on_generators
    ):
        if (
            domain.base_ring != coefficient_morphism.domain
            or codomain.base_ring != coefficient_morphism.codomain
        ):
            raise TypeError(
                "The coefficient homomorphism is not compatible with the given domain and codomain."
            )
        if len(action_on_generators) != domain.ngens:
            raise TypeError(
                "Cannot interpret {} as a homomorphism from {} to {}".format(
                    action_on_generators, domain, codomain
                )
            )
        if codomain.precision_cap > domain.precision_cap:
            raise ValueError(
                "Codomain precision cap is larger than domain precision cap."
            )
        for i in range(len(action_on_generators)):
            if action_on_generators[i].ring != codomain:
                raise TypeError(
                    "The elements of 'action_on_generators' are not members of the codomain."
                )

        self.coefficient_morphism = coefficient_morphism
        self.domain = domain
        self.codomain = codomain
        self.action_on_generators = action_on_generators

    @functools.cache
    def _call_on_generator_power(self, idx, e):
        return self.action_on_generators[idx]**ZZ(e)

    @functools.cache
    def _call_on_monomial(self, t):
        new_term = self.codomain.one
        for j in range(len(t.degrees)//2):
            new_term *= self._call_on_generator_power(t[2*j], t[2*j+1])
        return new_term

    def __call__(self, f):
        # TODO: this assumes that the input and output of
        # self.coefficient_morphism consist of elements of the data classes of
        # the base_rings.
        if isinstance(f,MonomialData):
            return self._call_on_monomial(f)

        x = self.codomain.zero
        for m,c in f.term_data():
            x += self.coefficient_morphism(self.domain(c)) * self._call_on_monomial(m)
        return x

    def __str__(self):
        return f"Homomorphism from {self.domain} to {self.codomain} defined by {self.action_on_generators} on generators."

    def __repr__(self):
        return self.__str__()
