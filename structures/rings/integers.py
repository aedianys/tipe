#!/usr/bin/python3
# -*- coding: utf-8 -*-

from structures.rings import CommutativeRing, CommutativeRingElement


class IntegersRing(CommutativeRing):
    """
    The ring of integers.
    """

    def __init__(self):
        pass

    def __eq__(self, other):
        return self.__class__() == other.__class__()

    def __call__(self, element):
        return Integers(element)

    def zero(self):
        return 0

    def one(self):
        return 1


class Integers(CommutativeRingElement, int):
    def ring(self):
        return IntegersRing()