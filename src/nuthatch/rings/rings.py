"""
RINGS

A module for commutative rings.
"""


class AbstractRing:
    def __init__(self, element_class, exact=True):
        self.element_class = element_class
        self.exact = exact

    def __call__(self, x):
        return self.element_class(self, x)
