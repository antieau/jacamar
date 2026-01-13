"""
RINGS

A module for commutative rings.

AUTHORS:
- Benjamin Antieau (2024).
"""

from nuthatch.rings.morphisms import IdentityRingMorphism


class AbstractRing:
    def __init__(self, element_class, exact=True):
        self.element_class = element_class

        # Whether or not the ring is exact.
        self.exact = exact

        # The following should be overridden for each concrete class.
        self.zero = None
        self.one = None

    def identity_morphism(self):
        return IdentityRingMorphism(self)

    def __call__(self, x):
        return self.element_class(self, x)
