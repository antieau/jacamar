"""
MORPHISMS

Base classes for morphisms of commutative rings.
"""


class AbstractRingMorphism:
    def __init__(self, domain, codomain):
        self.domain = domain
        self.codomain = codomain
