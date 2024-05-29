"""
TEST_POLYNOMIALS

Tests for the MonomialData, PolynomialData, Polynomial, and PolynomialRing classes.
"""

import pytest
import flint
from nuthatch.rings.polynomials import (
    MonomialData,
    PackedMonomialData,
    SparseMonomialData,
    PolynomialData,
    Polynomial,
    PolynomialRing,
    PolynomialRingMorphism,
)
from nuthatch.rings.integers import ZZ
from nuthatch.rings.rationals import QQ



class TestMonomialData:
    """Tests for the MonomialData class."""

    def test_from_packed_integer(self):
        """Tests that from_packed_int is not implemented."""
        assert MonomialData.from_packed_integer(5) == NotImplemented

    def test_from_sparse_tuple(self):
        """Tests that from_sparse_tuple is not implemented."""
        assert MonomialData.from_sparse_tuple((1, 2, 4, 7)) == NotImplemented

    def test_from_tuple(self):
        """Tests that from_tuple is not implemented."""
        assert MonomialData.from_tuple((0, 2, 0, 0, 7)) == NotImplemented



class TestSparseMonomial:
    """Tests for the SparseMonomialData class."""

    def test_empty(self):
        """Edge tests for the empty monomial."""
        assert SparseMonomialData(()) == SparseMonomialData(())
        m = SparseMonomialData((1, 2, 4, 7))
        assert m * SparseMonomialData(()) == m

    def test_eq(self):
        """Tests the __eq__ method."""
        assert SparseMonomialData((1, 2, 4, 7)) == SparseMonomialData((1, 2, 4, 7))
        assert SparseMonomialData((1, 2, 4, 7)) != SparseMonomialData((1, 2, 4, 8))

    def test_mul(self):
        """Tests the __mul__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        n = SparseMonomialData((0, 3, 3, 3, 4, 1))
        assert m * n == SparseMonomialData((0, 3, 1, 2, 3, 3, 4, 8))
        assert m * SparseMonomialData(()) == m

    def test_hash(self):
        """Tests the __hash__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        assert hash(m) == hash((1, 2, 4, 7))

    def test_str_and_repr(self):
        """Tests the __str__ method."""
        m = SparseMonomialData((1, 2, 4, 7))
        assert str(m) == str((1, 2, 4, 7))
        assert repr(m) == str(m)

    def test_class_methods(self):
        """Tests the class method constructors."""
        assert SparseMonomialData((1, 2, 4, 7)) == SparseMonomialData.from_sparse_tuple((1, 2, 4, 7))
        assert SparseMonomialData((1, 2, 4, 7)) == SparseMonomialData.from_tuple((0, 2, 0, 0, 7))
        assert SparseMonomialData((1, 2, 4, 7)) == SparseMonomialData.from_packed_integer(129127208515966992384)


class TestPackedMonomial:
    """Tests for the PackedMonomialData class."""

    def test_empty(self):
        """Edge tests for the empty monomial."""
        assert PackedMonomialData(0) == PackedMonomialData(0)
        m = PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        assert m * PackedMonomialData(0) == m

    def test_eq(self):
        """Tests the __eq__ method."""
        assert PackedMonomialData.from_sparse_tuple((1, 2, 4, 7)) == PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        assert PackedMonomialData.from_sparse_tuple((1, 2, 4, 7)) != PackedMonomialData.from_sparse_tuple((1, 2, 4, 8))

    def test_mul(self):
        """Tests the __mul__ method."""
        m = PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        n = PackedMonomialData.from_sparse_tuple((0, 3, 3, 3, 4, 1))
        assert m * n == PackedMonomialData.from_sparse_tuple((0, 3, 1, 2, 3, 3, 4, 8))

    def test_hash(self):
        """Tests the __hash__ method."""
        m = PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        # This could fail if the constant PACKING_BOUND from
        # nuthatch.rings.polynomials is changed.
        assert hash(m) == hash(129127208515966992384)

    def test_str_and_repr(self):
        """Tests the __str__ method."""
        m = PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        # This could fail if the constant PACKING_BOUND from
        # nuthatch.rings.polynomials is changed.
        assert str(m) == "129127208515966992384"
        assert repr(m) == str(m)

    def test_class_methods(self):
        """Tests the class method constructors."""
        assert PackedMonomialData(129127208515966992384) == PackedMonomialData.from_sparse_tuple((1, 2, 4, 7))
        assert PackedMonomialData(129127208515966992384) == PackedMonomialData.from_tuple((0, 2, 0, 0, 7))
        assert PackedMonomialData(129127208515966992384) == PackedMonomialData.from_packed_integer(129127208515966992384)


class TestPolynomialDataClass:
    """Tests for the PolynomialData class."""

    d = {SparseMonomialData((0, 1)): flint.fmpz(1), SparseMonomialData((1, 1)): flint.fmpz(1)}
    p = PolynomialData(ZZ, d)

    def test_mul(self):
        """Tests for the __mul__ method."""
        # (x+y)^2 == x^2 + 2xy + y^2
        assert self.p * self.p == PolynomialData(
            ZZ,
            {
                SparseMonomialData((0, 2)): flint.fmpz(1),
                SparseMonomialData((0, 1, 1, 1)): flint.fmpz(2),
                SparseMonomialData((1, 2)): flint.fmpz(1),
            },
        )

    def test_is_zero(self):
        """Tests the is_zero method."""
        assert not self.p.is_zero()
        assert PolynomialData(ZZ,{}).is_zero()
        assert PolynomialData(ZZ,{SparseMonomialData(0,1): ZZ.zero.data}).is_zero()
        q = PolynomialData(ZZ, {PackedMonomialData(23512): flint.fmpz(0), PackedMonomialData(515167141362451626134541345): flint.fmpz(0)})
        assert q.is_zero()
        assert not PolynomialData(ZZ, {PackedMonomialData(5): flint.fmpz(5)}).is_zero()

    def test_term_data(self):
        """Tests the term_data method."""
        assert self.d.items() == self.p.term_data()

    def test_printing(self):
        """Tests the printing methods."""
        assert str(self.p) == "{(0, 1): 1, (1, 1): 1}"
        assert repr(self.p) == str(self.p)


class TestPolynomialRing:
    """Tests for the PolynomialRing class."""

    r = PolynomialRing(base_ring=ZZ, ngens=3, prefix="x", packed=False)
    x0 = r.gens[0]
    x1 = r.gens[1]
    x2 = r.gens[2]
    a = r(-5) + x0 + x1 + x2

    def test_gens_data(self):
        """Tests the correct creation of the generators."""
        assert self.r.gens[0].data == PolynomialData(
            ZZ, {SparseMonomialData((0, 1)): flint.fmpz(1)}
        )

    def test_square(self):
        """Tests that (x0+x1)**2 == x0**2 + 2*x0x1 + x1**2."""
        assert (self.x0 + self.x1) ** ZZ(2) == (
            self.x0 ** ZZ(2) + ZZ(2) * self.x0 * self.x1 + self.x1 ** ZZ(2)
        )

    def test_constructor_from_int(self):
        """Tests r(5)."""
        assert str(self.r(5)) == "5"

    def test_constructor_from_base_ring(self):
        """Tests r(ZZ(5))."""
        assert str(self.r(ZZ(5))) == "5"

    def test_constructor_from_element(self):
        assert self.r(self.a) == self.a

    def test_constructor_from_element_data(self):
        assert self.r(self.a.data) == self.a

    def test_constructor_from_dict(self):
        """Tests dictionary constructor."""
        f = self.r({(1,1,2,1):ZZ(2),(0,4):ZZ(9)})
        assert f == ZZ(2) * self.x1 * self.x2 + ZZ(9) * self.x0**ZZ(4)

    def test_packed_sparse_str(self):
        """Tests that sparse and packed versions print the same."""
        m = str(self.a*self.a)
        s = PolynomialRing(base_ring=ZZ, ngens=3, prefix="x", packed=True)
        x0, x1, x2 = s.gens
        b = s(-5) + x0 + x1 + x2
        m == str(b*b)

    def test_sub(self):
        """Tests the sub method."""
        assert self.a - self.a == self.r.zero

    def test_negative_power(self):
        """Tests that powering fails with negative input."""
        with pytest.raises(TypeError):
            self.a**ZZ(-5)

    def test_zero_power(self):
        """Tests that powering by zero returns one."""
        assert self.a**ZZ(0) == self.r.one

    def test_power_failure(self):
        """Tests that powering by rational or polynomial fails."""
        with pytest.raises(TypeError):
            self.a**QQ(1,2)
        with pytest.raises(TypeError):
            self.a**self.a


class TestLayers:
    """Tests polynomial rings over polynomial rings."""
    r = PolynomialRing(base_ring=ZZ,ngens=2,prefix='x', packed=True)
    s = PolynomialRing(base_ring=r,ngens=2,prefix='y', packed=True)
    x0,x1 = r.gens
    y0,y1 = s.gens
    a = r(1) + x0 + x1
    b = s(1) + y0 + y1

    def test_evaluation(self):
        """Tests evaluation."""
        print(self.a)
        print(self.s(self.a)(0, 0))
        assert self.s(self.a)(0, 0) == self.a

    def test_multiplication(self):
        with pytest.raises(TypeError):
            self.b * self.a
        assert self.a*self.b == self.s(self.a) + self.a*self.y0 + self.a*self.y1
        assert self.s(self.a) * self.b == self.a*self.b


class TestCalls:
    """Test evaluation of polynomials."""
    s=PolynomialRing(base_ring=ZZ,ngens=4,prefix='x', packed=False)
    x0,x1,x2,x3 = s.gens
    f = x0+ZZ(2)*x1*x2+x3

    def test_evaluation_at_coorindates(self):
        """Tests evaluation at tuples of base_ring elements."""
        assert self.f(ZZ(1), ZZ(2), ZZ(3), ZZ(4)) == ZZ(17)

    def test_evaluation_at_monomial(self):
        """Tests evaluation at a SparseMonomialData instance."""
        assert self.f(SparseMonomialData((1, 1, 2, 1))) == ZZ(2)
        assert self.f(SparseMonomialData(1, 1, 2, 1)) == ZZ(2)

    def test_evaluation_at_implicit_monomial(self):
        """
        Tests evaluation at something that can be coerced into a
        SparseMonomialData instance.
        """
        assert self.f(1,1,2,1) == ZZ(0)

    def test_evaluation_at_other(self):
        """
        Tests that evaluating at something else is NotImplemented.
        """
        assert self.f(QQ(1,2)) == NotImplemented
        assert self.f(1,1,2,1,1,2) == NotImplemented


class TestPrinting:
    """Test printing in various situations."""
    r = PolynomialRing(base_ring=ZZ, ngens=2, prefix='x')
    x0,x1=r.gens
    f = r(-5)*x0 - ZZ(9)*x0**ZZ(5)*x1**ZZ(2)

    s = PolynomialRing(base_ring=r, ngens=2, prefix='y')
    y0,y1=s.gens
    g = s(x0+x1)*(y0+y1)

    h = s(x0) * (y0+y1)

    def test_integer_coefficients(self):
        """Tests printing with integer coefficients."""
        assert str(self.f) == "- 5*x0^1 - 9*x0^5*x1^2"
        assert repr(self.f) == str(self.f)

    def test_polynomial_coefficients(self):
        """Tests printing with polynomial coefficients."""
        assert str(self.g) == "(1*x0^1 + 1*x1^1)*y0^1 + (1*x0^1 + 1*x1^1)*y1^1"

    def test_simple_polynomial_coefficients(self):
        """Tests printing with simple polynomial coefficients."""
        assert str(self.h) == "1*x0^1*y0^1 + 1*x0^1*y1^1"

    def test_printing_ring(self):
        """Tests the printing of polynomial rings."""
        assert str(self.r) == "Ring of polynomials in 2 variables ['x0', 'x1'] with weights [1, 1] over The ring of Integers (via flint.fmpz)."
        assert str(self.s) == "Ring of polynomials in 2 variables ['y0', 'y1'] with weights [1, 1] over Ring of polynomials in 2 variables ['x0', 'x1'] with weights [1, 1] over The ring of Integers (via flint.fmpz)."
        assert repr(self.r) == str(self.r)

    def test_printing_zero(self):
        """Tests printing zero."""
        assert str(self.r.zero) == "0"


class TestWeights:
    """Tests weights."""
    r = PolynomialRing(base_ring=ZZ, ngens=2, prefix='x', weights=[2,3])

    def test_weights(self):
        """Tests that weights are correctly initialized."""
        assert self.r.weights == [2,3]

    def test_bad_weights(self):
        """
        Tests that a ValueError is raised if the number of weights is not ngens.
        """
        with pytest.raises(ValueError):
            s = PolynomialRing(base_ring=ZZ, ngens=2, prefix='x', weights=[1,2,3,4])


class TestMorphism:
    """Tests PolynomialRingMorphisms."""
    r = PolynomialRing(base_ring=ZZ, ngens=2, prefix='x')
    x0,x1=r.gens

    s = PolynomialRing(base_ring=ZZ, ngens=2, prefix='y')
    y0,y1=s.gens

    f0 = ZZ(2)*y0**ZZ(3) + y0*y1
    f1 = y0+y1

    f = PolynomialRingMorphism(domain = r, codomain = s, coefficient_morphism = ZZ.identity_morphism(), action_on_generators = [f0,f1])

    def test_morphism(self):
        """Tests evaluation of morphisms."""
        assert self.f(self.x0) == self.f0
        assert self.f(self.x1) == self.f1
        assert self.f(PackedMonomialData.from_sparse_tuple((0,1,1,1))) == self.f0*self.f1

    def test_base_ring_check(self):
        """Tests the detection of incompatible base ring morphism."""
        t = PolynomialRing(base_ring=QQ, ngens=2, prefix='z')
        z0,z1=t.gens
        with pytest.raises(TypeError):
            PolynomialRingMorphism(domain=self.r, codomain=t, coefficient_morphism=ZZ.identity_morphism(), action_on_generators=[z0,z1])

    def test_bad_ngens(self):
        """Tests the detection of the wrong number of generators."""
        with pytest.raises(TypeError):
            PolynomialRingMorphism(domain=self.r, codomain=self.s, coefficient_morphism=ZZ.identity_morphism(), action_on_generators=[self.y0,self.y1,self.f0])

    def test_bad_action(self):
        """Tests ring detection."""
        with pytest.raises(TypeError):
            PolynomialRingMorphism(domain = self.r, codomain = self.s, coefficient_morphism = ZZ.identity_morphism(), action_on_generators = [self.x0,self.x1])

    def test_printing(self):
        """Tests printing."""
        assert repr(self.f) == str(self.f)
        assert str(self.f) == "Homomorphism from Ring of polynomials in 2 variables ['x0', 'x1'] with weights [1, 1] over The ring of Integers (via flint.fmpz). to Ring of polynomials in 2 variables ['y0', 'y1'] with weights [1, 1] over The ring of Integers (via flint.fmpz). defined by [2*y0^3 + 1*y0^1*y1^1, 1*y0^1 + 1*y1^1] on generators."
