"""
TEST_MORPHISMS

Tests for the IdentityRingMorphism class.
"""

from nuthatch.rings.integers import ZZ, ZZ_py
from nuthatch.rings.rationals import QQ

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
        assert i(QQ(3,5)) == QQ(12,20)
