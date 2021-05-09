#!/usr/bin/python3
# -*- coding: utf-8 -*-

from rings import CommutativeRing, CommutativeRingElement
from maths.math_lib import euclide


class QuotientRing(CommutativeRing):
    def __init__(self, ring, modulus):
        self.modulus = modulus
        self.__ring = ring

    def __eq__(self, other):
        return (self.modulus, self.__ring) == (other.modulus, other.ring())

    def __call__(self, element):
        return QuotientRingElement(self, element)

    def zero(self):
        return QuotientRingElement(self, 0)

    def one(self):
        return QuotientRingElement(self, 1)
    
    def ring(self):
        return self.__ring


class QuotientRingElement(CommutativeRingElement):
    """
    A ring of residue classes of ring elements modulo a fixed ring element; the
    residue classes (quotient ring elements) support infix notation for ring
    operations, as well as mixed argument types.
    """

    def __init__(self, quotient_ring, representative):
        """
        Construct a new residue class representative modulo quotient_ring.modulus
        """

        self.__quotient_ring = quotient_ring

        if isinstance(representative, self.__class__):
            if self.__quotient_ring == representative.__quotient_ring:
                self.__remainder = representative.__remainder
            else:
                m = self.__quotient_ring.modulus
                self.__remainder = representative.__remainder % m
        else:
            m = self.__quotient_ring.modulus
            self.__remainder = representative % m

    def __str__(self):
        return f"{self.__remainder}"

    def remainder(self):
        return self.__remainder

    def modulus(self):
        return self.__quotient_ring.modulus

    def __bool__(self):
        return bool(self.__remainder)

    def __eq__(self, other):
        try:
            other = self.__quotient_ring(other)
            return self.__remainder == other.remainder()

        except TypeError:
            return False

    def __add__(self, other):
        other = self.__quotient_ring(other)
        return self.__quotient_ring(self.__remainder + other.remainder())

    def __neg__(self):
        return self.__quotient_ring(-self.__remainder)

    def __mul__(self, other):
        try:
            other = self.__quotient_ring(other)
            return self.__quotient_ring(self.__remainder * other.remainder())
        except:
            return other.__mul__(self)

    def __truediv__(self, other):
        if not other:
            raise ZeroDivisionError

        other = self.__quotient_ring(other)
        return self.__mul__(other.inverse())

    def __rtruediv__(self, other):
        return self.inverse() * other

    def inverse(self):
        if not self:
            raise ZeroDivisionError

        gcd, inverse, _ = euclide(self.__remainder, self.modulus())
        R = self.__quotient_ring.ring()

        if gcd == R.one():
            return self.__quotient_ring(inverse)
        else:
            raise ZeroDivisionError(gcd)

    def __int__(self):
        return int(self.__remainder)

    def ring(self):
        return self.__quotient_ring