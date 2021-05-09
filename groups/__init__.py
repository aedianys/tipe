#!/usr/bin/python3
# -*- coding: utf-8 -*-


class AbelianGroup:
    """
    A base class for abelian groups that provides default
    operator overloading.
    """

    def __neq__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        raise NotImplementedError

    def __call__(self, element):
        raise NotImplementedError

    def id(self):
        raise NotImplementedError


class AbelianGroupElement:
    """
    A base class for elements of abelian groups that provides default
    operator overloading.
    """

    def __neq__(self, other):
        return not self.__eq__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return -self.__sub__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, n):
        n = int(n)

        if n < 0:
            self = -self
            n = -n

        result = self.group().id()
        power = self

        while n != 0:
            if n % 2 == 1:
                result = result + power
            power = power + power
            n //= 2

        return result

    def __bool__(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError

    def group(self):
        raise NotImplementedError