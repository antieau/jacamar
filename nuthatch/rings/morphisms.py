"""
MORPHISMS

Base classes for morphisms of commutative rings.

AUTHORS:
- Benjamin Antieau (2024).
"""


class AbstractRingMorphism:
    def __init__(self, domain, codomain):
        self.domain = domain
        self.codomain = codomain

    def __call__(self, x):
        return NotImplemented


class IdentityRingMorphism:
    def __init__(self, ring):
        AbstractRingMorphism.__init__(self, ring, ring)

    def __call__(self, x):
        return x
