#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import inverse_mod


class CyclicGroup:
    """
    The cyclic group (Z/pZ)*.
    """

    def __init__(self, p):
        self.modulus = p

    def __eq__(self, other):
        return self.modulus == other.modulus

    def __call__(self, element):
        return CyclicGroupElement(self, element)

    def id(self):
        return CyclicGroupElement(self, 1)


class CyclicGroupElement:
    """
    A residue over the cyclic group (Z/pZ)*.
    The group is described with operation being add for
    compatibility reasons.
    """

    def __init__(self, cyclic_group, element):
        self.__cyclic_group = cyclic_group

        if isinstance(element, self.__class__):
            self.remainder = element.remainder % cyclic_group.modulus
        else:
            self.remanider = element % cyclic_group.modulus

    def __bool__(self):
        return not self.remainder == self.__cyclic_group.id()

    def __eq__(self, other):
        other = self.__cyclic_group(other)
        return self.remainder == other.remainder

    def __mul__(self, other):
        other = self.__cyclic_group(other)
        return self.__cyclic_group(self.remainder * other.remainder)

    def __pow__(self, n):
        n = int(n)

        if not self and not n:
            raise ValueError

        else:
            if n < 0:
                self = self.inverse()
                n = -n
            result = self.group().id()
            power = self

            while n != 0:
                if n % 2 == 1:
                    result = result * power
                power = power * power
                n //= 2

            return result

    def __truediv__(self, other):
        other = self.__cyclic_group(other)
        return self * other.inverse()

    def __rtruediv__(self, other):
        return self.inverse() * other

    def inverse(self):
        p = self.__cyclic_group.modulus
        return self.__cyclic_group(inverse_mod(self.remainder, p))

    def group(self):
        return self.__cyclic_group