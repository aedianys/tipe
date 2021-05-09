#!/usr/bin/python3
# -*- coding: utf-8 -*-

from rings import CommutativeRingElement


class FieldElement(CommutativeRingElement):
    """
    A base class for field elements that provides default operator overloading.
    """

    def __truediv__(self, other):
        other = self.field()(other)
        return self.__mul__(other.inverse())

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, n):
        n = int(n)
        if n < 0:
            return self.inverse().__pow__(-n)
        else:
            return super().__pow__(n)

    def inverse(self):
        raise NotImplementedError
    
    def ring(self):
        raise NotImplementedError
    
    def field(self):
        return self.ring()