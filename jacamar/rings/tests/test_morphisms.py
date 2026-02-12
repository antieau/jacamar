"""
TEST_MORPHISMS

Tests for the IdentityRingMorphism class.
"""

from jacamar.rings.integers import ZZ, ZZ_py
from jacamar.rings.morphisms import AbstractRingMorphism
from jacamar.rings.rationals import QQ


class TestIdentity:
    """Tests for the Integer class."""

    def test_identity_zz(self):
        i = ZZ.identity_morphism()
        assert i(ZZ(5)) == ZZ(5)

    def test_identity_zzpy(self):
        i = ZZ_py.identity_morphism()
        assert i(ZZ_py(5)) == ZZ_py(5)

    def test_identity_qq(self):
        i = QQ.identity_morphism()
        assert i(QQ(5)) == QQ(5)
        assert i(QQ(3, 5)) == QQ(12, 20)


class TestGeneric:
    """Tests for the AbstractRingMorphism class."""

    def test_not_implemented(self):
        """Tests that __call__ is not implemented."""
        assert AbstractRingMorphism(ZZ, ZZ)(ZZ(4)) == NotImplemented
