"""
POLYNOMIALS

Classes for generic polynomial rings. The exponents are designed based on the
`ETuple` class from `SAGE`. In particular, they are sparse. In a polynomial
ring on x0, x1, x2, x3, a monomial x0*x2^3*x3^7 would be encoded as
(0,1,2,3,3,7), where the even indices indicate the present variables and the
odd indices indicate their powers. The empty tuple () represents the monomial
of 1.

We attempt to trim away zero terms whenever possible.

AUTHORS:
- Benjamin Antieau (2024): initial version.
"""

import functools
import flint
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.rings import AbstractRing
from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.morphisms import AbstractRingMorphism


class MonomialData:
    """
    The class of a monomial in a PolynomialData.

    INPUT:
    - `t`     -- a tuple of Python integers

    TODO:
    - make this into C or Rust code
    """

    def __init__(self, *t):
        if len(t) == 1 & isinstance(t,tuple):
            self.exponent = t[0]
        else:
            self.exponent = tuple(t)

    def degree(self):
        """Sum over the odd entries in the exponent tuple."""
        return sum(self.exponent[1::2])

    def in_ith_variable(self, i):
        j = 0
        while j < len(self.exponent):
            if self.exponent[j] < i:
                j += 2
                continue
            elif self.exponent[j] == i:
                return self.exponent[j + 1]
            else:  # t[j] > i
                return 0
        return 0

    def __hash__(self):
        return hash(self.exponent)

    def __str__(self):
        return str(self.exponent)

    def __repr__(self):
        return str(self)

    def __mul__(self, other):
        new_list = []
        self_index = 0
        other_index = 0
        self_len = len(self.exponent)
        other_len = len(other.exponent)
        while self_index < self_len and other_index < other_len:
            if self.exponent[self_index] < other.exponent[other_index]:
                new_list.extend(self.exponent[self_index : self_index + 2])
                self_index += 2
            elif self.exponent[self_index] == other.exponent[other_index]:
                x = self.exponent[self_index + 1] + other.exponent[other_index + 1]
                if x != 0:
                    new_list.append(self.exponent[self_index])
                    new_list.append(x)
                self_index += 2
                other_index += 2
            else:
                new_list.extend(other.exponent[other_index : other_index + 2])
                other_index += 2

        # Now that we've reached the end of at least one list, we add
        # everything from the other list.
        while self_index < self_len:
            new_list.extend(self.exponent[self_index : self_index + 2])
            self_index += 2

        while other_index < other_len:
            new_list.extend(other.exponent[other_index : other_index + 2])
            other_index += 2

        return other.__class__(tuple(new_list))

    def __eq__(self, other):
        return self.exponent == other.exponent
        


class PolynomialData:
    """
    The underlying data class of polynomials.
    The `data` input should be a dictionary `{m:c}` where `m` is a `MonomialData`
    and `c` is an element of the data_class of the element_class of the `base_ring`.
    """

    def __init__(self, base_ring, monomial_dictionary):
        self.base_ring = base_ring
        self.monomial_dictionary = monomial_dictionary

    def is_zero(self):
        if len(self.monomial_dictionary) == 0:
            return True
        for _, c in self.monomial_dictionary.items():
            if c != self.base_ring.zero.data:
                return False
        return True

    def term_data(self):
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
                if new_dict[m] == other.base_ring.zero.data:
                    del new_dict[m]
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
        new_dict = {}
        for m, c in self.monomial_dictionary.items():
            for n, d in other.monomial_dictionary.items():
                k = m * n
                e = c * d
                try:
                    new_dict[k] += e
                except KeyError:
                    new_dict[k] = e
        final_dict = {}
        for m, c in new_dict.items():
            if c != other.base_ring.zero.data:
                final_dict[m] = c
        return other.__class__(other.base_ring, final_dict)

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
            try:
                return self.monomial_dictionary[arg]
            except KeyError:
                return self.base_ring.zero.data

    def evaluate(self, args):
        """Returns f(a_1,...,a_n) where arguments=[a_1,...,a_n]."""
        x = self.base_ring.zero.data
        for key, value in self.monomial_dictionary.items():
            monomial_part = self.base_ring.one.data
            for i in range(len(key.exponent) // 2):
                monomial_part *= args[key.exponent[2*i]]**key.exponent[2*i+1]
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

        return self.base_ring(self.data(MonomialData(tuple(args))))

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
            term = f"- {-self.data.monomial_dictionary[keys[0]]}"
        else:
            term = f"{self.data.monomial_dictionary[keys[0]]}"
            for j in range(len(key.exponent) // 2):
                term += "*"
                term += self.ring._prefix
                term += str(key.exponent[2 * j])
                term += f"^{key.exponent[2*j+1]}"
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
                term = f" - {-value}"
            else:
                term = f" + {value}"
            for j in range(len(key.exponent) // 2):
                term += "*"
                term += self.ring._prefix
                term += str(key.exponent[2 * j])
                term += f"^{key.exponent[2*j+1]}"
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
    ):
        self.base_ring = base_ring
        self.ngens = ngens
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
                        self.base_ring, {MonomialData((i, 1)): self.base_ring.one.data}
                    ),
                )
            )

        # Initialize one and zero since these are used so often.
        self.one = self(1)
        self.zero = self(0)

        AbstractRing.__init__(self, Polynomial, exact=base_ring.exact)

    def __str__(self):
        return f"Ring of polynomials in {self.ngens} variables {self._names} with weights {self.weights} over {self.base_ring}"

    def __repr__(self):
        return self.__str__()

    def __call__(self, data):
        if isinstance(data, int):
            if data == 0:
                return Polynomial(self, PolynomialData(self.base_ring, {}))
            else:
                return Polynomial(
                    self,
                    PolynomialData(
                        self.base_ring,
                        {MonomialData(()): self.base_ring.element_class.data_class(data)},
                    ),
                )
        if isinstance(data, dict):
            new_dict = {}
            for key, value in data.items():
                new_dict[MonomialData(key)] = value.data
            return Polynomial(self, PolynomialData(self.base_ring, new_dict))
        raise TypeError("No known constructor for input data.")


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
        for j in range(len(t.exponent)//2):
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
