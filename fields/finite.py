#!/usr/bin/python3
# -*- coding: utf-8 -*-

from rings.quotients import QuotientRing
from rings.integers import IntegersRing


class FiniteField(QuotientRing):
    """
    The finite field with p^k elements. Support only k=1.
    """

    def __init__(self, p):
        super().__init__(IntegersRing(), p)

    def characteristic(self):
        return self.modulus

    def size(self):
        return self.characteristic()