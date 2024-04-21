"""
ELEMENTS

Base classes for elements in commutative rings.
"""

class AbstractRingElement:
    def __init__(self,ring):
        self.ring = ring
