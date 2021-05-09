#!/usr/bin/python3
# -*- coding: utf-8 -*-


class CommutativeRing:
    """
    Template class for commutative rings.
    """

    def __neq__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        raise NotImplementedError

    def __call__(self, element):
        raise NotImplementedError

    def zero(self):
        raise NotImplementedError

    def one(self):
        raise NotImplementedError


class CommutativeRingElement:
    """
    A base class for elements of commutative rings that provides default
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

    def __rdivmod__(self, other):
        return divmod(self.ring()(other), self)

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __rfloordiv__(self, other):
        return divmod(other, self)[0]

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __rmod__(self, other):
        return divmod(other, self)[1]

    def __pow__(self, n):
        n = int(n)
        assert n >= 0

        if not self and not n:
            raise ValueError

        else:
            result = self.ring().one()
            power = self

            while n != 0:
                if n % 2 == 1:
                    result = result * power
                power = power * power
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

    def __mul__(self, other):
        raise NotImplementedError

    def __divmod__(self, other):
        raise NotImplementedError

    def ring(self):
        raise NotImplementedError