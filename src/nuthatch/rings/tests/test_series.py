"""
TEST_SERIES

Tests for the _Series, Series, and PowerSeriesRing classes.
"""

import pytest
import flint
from nuthatch.rings.series import (
    _Series,
    Series,
    PowerSeriesRing,
)
from nuthatch.rings.integers import ZZ


class TestSeries:
    """Tests for the Series class."""
    s=PowerSeriesRing(base_ring=ZZ,ngens=4,prefix='x',precision_cap=5)
    x0,x1,x2,x3 = s.gens

    def test_mul(self):
        """Tests the __mul__ and __rmul__ methods."""
        assert (self.x0+self.x1)**ZZ(2) == self.x0**ZZ(2) + ZZ(2)*self.x0*self.x1 + self.x1**ZZ(2)
