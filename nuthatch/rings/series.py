"""
SERIES

This is a module for absolutely capped precision power series in multiple variables with arbitrary
positive integer weights.

To model k[[x_1,...,x_n]] we include it into k[y_1,...,y_n][[T]] x_i maps to y_i*T^w_i
where w_i is the weight of x_i. Under this inclusion, the power series ring is closed under all
arithmetic operations in the ambient ring.

Exposed classes:
 - PowerSeriesRing,
 - PowerSeriesRingHomomorphism, and
 - Series.

AUTHORS:
- Benjamin Antieau (2024): initial version.

"""

import functools
import itertools
from nuthatch.rings.elements import AbstractRingElement
from nuthatch.rings.integers import ZZ
from nuthatch.rings.rings import AbstractRing
from nuthatch.rings.polynomials import (
    PolynomialRing,
    PolynomialData,
    Polynomial,
    MonomialData,
)
from nuthatch.rings.morphisms import AbstractRingMorphism


class SeriesData:
    """
    Data class for elements of weighted power series ring.
    """

    def __init__(self, base_ring, term_list, precision):
        """
        The argument term_list is an ordered list [(N,coefficient)] representing g*T^N
        where N is a Python integer and g is a PolynomialData of degree N.
        """
        self.base_ring = base_ring
        self.term_list = term_list
        self.precision = precision

    def __str__(self):
        return str(self.term_list)

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        """
        Add self to other, taking other.precision as the result's precision
        cap.
        """
        term_dict = dict(other.term_list)
        for deg, coeff in self.term_list:
            if deg < other.precision:
                if deg in term_dict:
                    term_dict[deg] += coeff
                else:
                    term_dict[deg] = coeff
        sorted_degrees = list(term_dict)
        sorted_degrees.sort()
        newterm_list = []
        for deg in sorted_degrees:
            if not term_dict[deg].is_zero():
                newterm_list.append((deg, term_dict[deg]))
        return other.__class__(self.base_ring, newterm_list, other.precision)

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return self.__class__(
            self.base_ring,
            [(deg, -coeff) for deg, coeff in self.term_list],
            self.precision,
        )

    def __eq__(self, other):
        return (self.term_list == other.term_list) & (self.precision == other.precision)

    def __ne__(self, other):
        return not self == other

    def is_unit(self):
        if (
            self.term_list[0][0] != 0
            or not self.term_list[0][1].constant_coefficient().is_unit()
        ):
            return False
        else:
            return True

    def __rmul__(self, other):
        # TODO: this could result in zero terms.
        return self.__class__(
            self.base_ring,
            [(deg, other * coefficient) for deg, coefficient in self.term_list],
            self.precision,
        )

    def __mul__(self, other):
        """
        Multiplication assumes that the input term_lists are sorted by degree.

        TODO: this can be parallelized in an obvious way.
        """
        term_dict = {}
        for deg, coeff in self.term_list:
            for deg1, coeff1 in other.term_list:
                if deg + deg1 < other.precision:
                    if deg + deg1 in term_dict:
                        term_dict[deg + deg1] += coeff * coeff1
                    else:
                        term_dict[deg + deg1] = coeff * coeff1
                else:
                    break
        sorted_degrees = list(term_dict)
        sorted_degrees.sort()
        newterm_list = []
        for deg in sorted_degrees:
            if not term_dict[deg].is_zero():
                newterm_list.append((deg, term_dict[deg]))
        return self.__class__(self.base_ring, newterm_list, other.precision)

    def __pow__(self, n):
        # This is probably not very pythonic. But, it is the only
        # way I can figure out how to get self**n to work. The only other option
        # would be to implement __rpow__(self,n) and apply it n**self,
        # which would make for garbage code.
        # From sage/arith/power.pyx
        apow = self
        while not n & 1:  # While even...
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

    def _naive_power(self, n):
        pass

    def is_homogeneous(self):
        """
        Constants and zero count as homogeneous.
        """
        non_zero_terms = 0
        for _, coeff in self.term_list:
            if not coeff.is_zero():
                non_zero_terms += 1
            if non_zero_terms > 1:
                return False
        return True

    def degree(self):
        degree = -1
        for deg, coeff in self.term_list:
            if not coeff.is_zero():
                degree = deg
        return degree


class Series(AbstractRingElement):
    data_class = SeriesData

    def __init__(self, ring, data):
        self.base_ring = ring.base_ring
        AbstractRingElement.__init__(
            self,
            ring,
            data,
        )

    def term_data(self, a=None, b=None):
        """
        An iterator returning the terms of the series, packaged as tuples (m,c)
        where m is a MonomialData and c is an element of the base_ring.

        If a and b are provided, return only those terms of total degree
        between a and b.
        """
        iterators = []
        for n, g in self.data.term_list:
            if a <= n <= b:
                iterators.append(g.term_data())
        return itertools.chain.from_iterable(iterators)

    def constant_coefficient(self):
        if self.data.term_list[0][0] != 0:
            return self.base_ring.zero
        else:
            for m, c in self.data.term_list[0][1].monomial_dictionary.items():
                return self.base_ring(c)

    def is_unit(self):
        if self.data.term_list[0][0] == 0:
            for m, c in self.data.term_list[0][1].monomial_dictionary.items():
                if self.base_ring(c).is_unit():
                    return True
        return False

    def inverse(self):
        newton_approximation_inverse = self.ring(self.constant_coefficient())
        N = 1
        while N < self.ring.precision_cap:
            N = 2 * N
            newton_approximation_inverse = (
                newton_approximation_inverse
                + newton_approximation_inverse
                * (self.ring.one - self * newton_approximation_inverse)
            )
        return newton_approximation_inverse

    def __rmul__(self, other):
        return other.__class__(other.ring, self.data * other.data)

    def __str__(self):
        """
        This will return a string representing the underlying element
        in whatever order the coefficient_dictionary has stored its keys in. In particular
        it is not guaranteed to be in lexicographic ordering or anything else useful.
        """
        if len(self.data.term_list) == 0:
            return f"0 + O(F^{self.ring.precision_cap})"
        else:
            return_str = ""
            for _, coef in self.data.term_list:
                return_str += str(self.ring._polynomial_ring(coef)) + " + "
            return return_str + f"O(F^{self.ring.precision_cap})"

    def __repr__(self):
        return self.__str__()


class PowerSeriesRing(AbstractRing):
    """
    Weighted power series ring.
    """

    element_class = Series

    def __init__(
        self,
        *,
        base_ring,
        ngens,
        prefix,
        precision_cap,
        weights=None,
    ):
        """
        Returns a new PowerSeriesRing.

        base_ring - a commutative ring
        precision_cap - an integer giving the precision for the prefix-adic coefficients
        ngens - the number of variables
        weights - a tuple of *positive* Python integer weights, one for each variable
        prefix - a string giving the prefix for the variables
        names - a list of variable names; this overrides prefix if given

        """
        self.base_ring = base_ring
        self.precision_cap = precision_cap
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
            for i in weights:
                if i <= 0:
                    raise ValueError("Weights must be positive.")
            self.weights = weights

        # Initialize one and zero since these are used so often.
        AbstractRing.__init__(self, Series, exact=base_ring.exact)

        self._polynomial_ring = PolynomialRing(
            base_ring=self.base_ring,
            ngens=self.ngens,
            prefix=self._prefix,
            packed=False,
        )

        self.one = self(1)
        self.zero = self(0)

        self.gens = []
        for generator in self._polynomial_ring.gens:
            # We have to coerce the SAGE integer into a Python integer.
            self.gens.append(self._unflatten(generator))

    def __call__(self, data):
        if self._polynomial_ring._special:
            self._polynomial_ring = self._polynomial_ring.to_generic()
            # if isinstance(data, int):
            #     if data == 0:
            #         return Series(self, SeriesData(self.base_ring, [], self.precision_cap))
            #     else:
            #         return Series(
            #             self,
            #             SeriesData(
            #                 self.base_ring,
            #                 [
            #                     (
            #                         0,
            #                         self._polynomial_ring.element_class.data_class(
            #                             data,
            #                             self._polynomial_ring.ctx
            #                         ),
            #                     ),
            #                 ],
            #                 self.precision_cap,
            #             ),
            #         )

            # if isinstance(data, self.base_ring.element_class):
            #     return Series(
            #         self,
            #         SeriesData(
            #             self.base_ring,
            #             [
            #                 (
            #                     0,
            #                     self._polynomial_ring.element_class.data_class(
            #                         self.base_ring,
            #                         {
            #                             self._polynomial_ring._monomial_class.from_sparse_tuple(
            #                                 tuple()
            #                             ): data.data
            #                         },
            #                     ),
            #                 ),
            #             ],
            #             self.precision_cap,
            #         ),
            #     )

            # if isinstance(data, self.base_ring.element_class.data_class):
            #     return Series(
            #         self,
            #         SeriesData(
            #             self.base_ring,
            #             [
            #                 (
            #                     0,
            #                     self._polynomial_ring.element_class.data_class(
            #                         self.base_ring,
            #                         {
            #                             self._polynomial_ring._monomial_class.from_sparse_tuple(
            #                                 tuple()
            #                             ): data
            #                         },
            #                     ),
            #                 ),
            #             ],
            #             self.precision_cap,
            #         ),
            #     )

            # return Series(self, data)

        if isinstance(data, int):
            if data == 0:
                return Series(self, SeriesData(self.base_ring, [], self.precision_cap))
            else:
                return Series(
                    self,
                    SeriesData(
                        self.base_ring,
                        [
                            (
                                0,
                                self._polynomial_ring.element_class.data_class(
                                    self.base_ring,
                                    {
                                        self._polynomial_ring._monomial_class.from_sparse_tuple(
                                            tuple()
                                        ): self.base_ring(
                                            data
                                        ).data
                                    },
                                ),
                            ),
                        ],
                        self.precision_cap,
                    ),
                )

        if isinstance(data, self.base_ring.element_class):
            return Series(
                self,
                SeriesData(
                    self.base_ring,
                    [
                        (
                            0,
                            self._polynomial_ring.element_class.data_class(
                                self.base_ring,
                                {
                                    self._polynomial_ring._monomial_class.from_sparse_tuple(
                                        tuple()
                                    ): data.data
                                },
                            ),
                        ),
                    ],
                    self.precision_cap,
                ),
            )

        if isinstance(data, self.base_ring.element_class.data_class):
            return Series(
                self,
                SeriesData(
                    self.base_ring,
                    [
                        (
                            0,
                            self._polynomial_ring.element_class.data_class(
                                self.base_ring,
                                {
                                    self._polynomial_ring._monomial_class.from_sparse_tuple(
                                        tuple()
                                    ): data
                                },
                            ),
                        ),
                    ],
                    self.precision_cap,
                ),
            )

        return Series(self, data)

    def _flatten(self, element):
        """
        Returns the element as an element in the underlying polynomial ring.
        """
        if element.ring != self:
            raise TypeError("Input must be a member of self.")
        out = self._polynomial_ring.zero
        for _, coeff in element.term_list:
            out += coeff
        return out

    def _unflatten_data(self, flat_polynomial_data):
        """
        Takes a PolynomialData instance and returns a SeriesData instance.
        """
        out_degrees = {}
        for m, c in flat_polynomial_data.monomial_dictionary.items():
            deg = 0
            for i in range(len(m.degrees) // 2):
                deg += self.weights[2 * i] * m.degrees[2 * i + 1]
            if deg < self.precision_cap:
                if deg in out_degrees:
                    out_degrees[deg] += self._polynomial_ring.element_class.data_class(
                        self.base_ring, {m: c}
                    )
                else:
                    out_degrees[deg] = self._polynomial_ring.element_class.data_class(
                        self.base_ring, {m: c}
                    )
        out_degrees_list = list(out_degrees.keys())
        out_degrees_list.sort()
        term_list = []
        for deg in out_degrees_list:
            term_list.append((deg, out_degrees[deg]))
        return self.element_class.data_class(
            self.base_ring, term_list, self.precision_cap
        )

    def _unflatten(self, flat_polynomial):
        """
        Takes an element of self._polynomial_ring and unpacks it into a power series,
        ignoring terms of degree at or above the precision cap.
        """
        return self.element_class(self, self._unflatten_data(flat_polynomial.data))

    def __str__(self):
        return (
            f"Absolute precision power series ring in {self._names} with weights {self.weights} "
            f"and with absolutely capped precision {self.precision_cap} over {self.base_ring}"
        )

    def __repr__(self):
        return self.__str__()


class PowerSeriesRingMorphism(AbstractRingMorphism):
    """
    Homomorphisms of weighted power series ring.
    """

    def __init__(self, *, domain, codomain, coefficient_morphism, action_on_generators):
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
        return self.action_on_generators[idx] ** ZZ(e)

    @functools.cache
    def _call_on_monomial(self, t):
        new_term = self.codomain.one
        for j in range(len(t.degrees) // 2):
            new_term *= self._call_on_generator_power(
                t.degrees[2 * j], t.degrees[2 * j + 1]
            )
        return new_term

    def __call__(self, f):
        # TODO: this assumes that the input and output of
        # self.coefficient_morphism consist of elements of the data classes of
        # the base_rings.
        if isinstance(f, MonomialData):
            return self._call_on_monomial(f)

        x = self.codomain.zero
        for m, c in f.term_data():
            x += self.coefficient_morphism(self.domain(c)) * self._call_on_monomial(m)
        return x

    def __str__(self):
        return f"Homomorphism from {self.domain} to {self.codomain} defined by {self.action_on_generators} on generators."

    def __repr__(self):
        return self.__str__()
