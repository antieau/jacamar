"""
TEST_SERIES

Tests for the SeriesData, Series, and PowerSeriesRing classes.
"""

import pytest
import flint
from jacamar.rings.polynomials import (
    PackedMonomialData,
    SparseMonomialData,
    PolynomialData,
    Polynomial,
    PolynomialRing,
)
from jacamar.rings.series import (
    SeriesData,
    Series,
    PowerSeriesRing,
)
from jacamar.rings.integers import ZZ


class TestSeriesData:
    """Tests for the SeriesData class."""

    a = PolynomialData(
        ZZ,
        {
            SparseMonomialData((1, 1, 3, 4)): ZZ(5).data,
            SparseMonomialData((0, 2, 5, 3)): ZZ(99).data,
        },
    )
    b = PolynomialData(
        ZZ,
        {
            SparseMonomialData((1, 4, 3, 4)): ZZ(5).data,
            SparseMonomialData((0, 3, 5, 5)): ZZ(99).data,
        },
    )
    f = SeriesData(ZZ, [(5, a)], 20)
    g = SeriesData(ZZ, [(8, b)], 20)
    h = SeriesData(ZZ, [(5, a), (8, b)], 20)

    def test_printing(self):
        """Tests printing."""
        assert str(self.f) == str([(5, self.a)])
        assert repr(self.f) == str(self.f)

    def test_add(self):
        """Tests __add__."""
        assert self.f + self.g == self.h

    def test_sub(self):
        """Tests __sub__."""
        assert self.f - self.g == SeriesData(ZZ, [(5, self.a), (8, -self.b)], 20)

    def test_ne(self):
        """Tests __ne__."""
        assert not self.f + self.g != self.h
        assert self.f != self.h

    def test_mul(self):
        """Tests __mul__."""
        assert self.f * self.g == SeriesData(ZZ, [(13, self.a * self.b)], 20)
        assert self.g * self.g * self.g == SeriesData(ZZ, [], 20)

        x = PolynomialData(ZZ, {SparseMonomialData((0, 1)): ZZ(1).data})
        y = PolynomialData(ZZ, {SparseMonomialData((1, 1)): ZZ(1).data})
        w = SeriesData(
            ZZ,
            [(0, PolynomialData(ZZ, {SparseMonomialData(()): ZZ(1).data})), (1, x)],
            20,
        )
        z = SeriesData(
            ZZ,
            [(0, PolynomialData(ZZ, {SparseMonomialData(()): ZZ(1).data})), (1, y)],
            20,
        )
        f = w + z
        assert f**2 == SeriesData(
            ZZ,
            [
                (0, PolynomialData(ZZ, {SparseMonomialData(()): ZZ(4).data})),
                (1, ZZ(4).data * (x + y)),
                (2, x * x + ZZ(2).data * x * y + y * y),
            ],
            20,
        )


class TestSeries:
    """Tests for the Series class."""

    s = PowerSeriesRing(base_ring=ZZ, ngens=4, prefix="x", precision_cap=5)
    x0, x1, x2, x3 = s.gens

    def test_mul(self):
        """Tests the __mul__ and __rmul__ methods."""
        assert (self.x0 + self.x1) ** ZZ(2) == self.x0 ** ZZ(2) + ZZ(
            2
        ) * self.x0 * self.x1 + self.x1 ** ZZ(2)
