#!/usr/bin/python3
# -*- coding: utf-8 -*-

from structures.rings import CommutativeRing
from structures.fields import FieldElement


class FractionField(CommutativeRing):
    """
    A field of formal fractions.
    """

    def __init__(self, field):
        self.field = field

    def __eq__(self, other):
        return self.field == other.field

    def __call__(self, numerator, denominator=None):
        return FractionFieldElement(self, numerator, denominator)

    def zero(self):
        return self(self.field.zero())

    def one(self):
        return self(self.field.one())

    def ring(self):
        return self.field

    def field(self):
        return self.field


class FractionFieldElement(FieldElement):
    """
    Formal fractions.
    """

    def __init__(self, fraction_field, numerator, denominator=None):
        self.__fraction_field = fraction_field

        if isinstance(numerator, self.__class__):
            # Copy an instance
            self.__numerator = numerator.__numerator
            self.__denominator = numerator.__denominator

        else:
            if denominator is None:
                denominator = self.__fraction_field.ring().one()
            if not denominator:
                raise ZeroDivisionError
            self.__numerator = self.__fraction_field.ring()(numerator)
            self.__denominator = self.__fraction_field.ring()(denominator)

    def numbers(self):
        return self.__numerator, self.__denominator

    def numerator(self):
        return self.__numerator

    def denominator(self):
        return self.__denominator

    def __str__(self):
        if self.__denominator == self.__fraction_field.ring().one():
            return f"{self.__numerator}"
        else:
            return f"({self.__numerator}) / ({self.__denominator})"

    def __bool__(self):
        return bool(self.__numerator)

    def __eq__(self, other):
        return (
            self.__numerator * other.__denominator
            == other.__numerator * self.__denominator
        )

    def __add__(self, other):
        other = self.ring()(other)
        numerator = (
            self.__numerator * other.__denominator
            + self.__denominator * other.__numerator
        )
        denominator = self.__denominator * other.__denominator
        return self.__fraction_field(numerator, denominator)

    def __neg__(self):
        return self.__fraction_field(-self.__numerator, self.__denominator)

    def __mul__(self, other):
        other = self.ring()(other)
        numerator = self.__numerator * other.__numerator
        denominator = self.__denominator * other.__denominator
        return self.__fraction_field(numerator, denominator)

    def inverse(self):
        if not self:
            raise ZeroDivisionError

        return self.__fraction_field(self.__denominator, self.__numerator)

    def ring(self):
        return self.__fraction_field