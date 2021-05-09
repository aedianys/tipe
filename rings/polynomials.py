#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.polynomial_multiplication import FPM
from rings import CommutativeRing, CommutativeRingElement
from functools import reduce


class Polynomials(CommutativeRing):
    """
    A ring of polynomials in one indeterminate.
    """

    def __init__(self, field):
        self.field = field

    def __eq__(self, other):
        return self.field == other.field

    def __call__(self, element_description, *further_coefficients):
        return Polynomial(self, element_description, *further_coefficients)

    def zero(self):
        return self(())

    def one(self):
        return self(1)

    def x(self):
        return self(0, 1)


class Polynomial(CommutativeRingElement):
    """
    A polynomial in one indeterminate.
    """

    def __init__(self, polynomials, element_description, *further_coefficients):
        # List of coefficients is in ascending order without leading zeros.
        # Example: x^2 + 2x + 5 = [5, 2, 1]
        if isinstance(element_description, self.__class__):  # copy
            self.__coefficients = element_description.__coefficients

        else:
            if type(element_description) in [list, tuple]:
                coefficients = element_description
            else:
                coefficients = [element_description] + list(further_coefficients)

            self.__coefficients = [polynomials.field(c) for c in coefficients]

        self.__strip_leading_zeroes()
        self.__polynomials = polynomials

    def __strip_leading_zeroes(self):
        while len(self.__coefficients) > 0 and not self.__coefficients[-1]:
            self.__coefficients.pop()

    def __str__(self):
        if not self:
            return f"{self.ring().field.zero()}"
        else:

            def monomial_str(c, k):
                if c == self.__polynomials.field.zero():
                    return ""
                elif c == self.__polynomials.field.one():
                    if k == 0:
                        return f"{c}"
                    elif k == 1:
                        return "x"
                    else:
                        return f"x^{k}"
                else:
                    if k == 0:
                        return f"{c}"
                    elif k == 1:
                        return f"{c}x"
                    else:
                        return f"{c}x^{k}"

            monomials = [monomial_str(c, k) for k, c in enumerate(self.__coefficients)]
            monomials_str = list(filter(lambda s: s != "", monomials))
            return reduce(lambda x, y: x + " + " + y, reversed(monomials_str))

    def coefficients(self):
        return self.__coefficients[:]

    def leading_coefficient(self):
        if self.__coefficients:
            return self.__coefficients[-1]
        else:
            return self.__polynomials.field.zero()

    def degree(self):
        if not self.__coefficients:
            return -(2 ** 30)
        return len(self.__coefficients) - 1

    def __bool__(self):
        return bool(self.__coefficients)

    def __eq__(self, other):
        try:
            other = self.__polynomials(other)
            if self.degree() != other.degree():
                return False

            for x, y in zip(self.__coefficients, other.coefficients()):
                if x != y:
                    return False

            return True

        except TypeError:
            return False

    def __add__(self, other):
        other = self.__polynomials(other)

        if not self:
            return other
        elif not other:
            return self
        else:
            zero = self.__polynomials.field.zero()
            max_length = max(self.degree(), other.degree())
            padded_self = self.__coefficients + [zero] * (max_length - self.degree())
            padded_other = other.coefficients() + [zero] * (max_length - other.degree())

            coefficient_sums = [x + y for x, y in zip(padded_self, padded_other)]
            return self.__polynomials(coefficient_sums)

    def __neg__(self):
        return self.__polynomials([-c for c in self.__coefficients])

    def __mul__(self, other):
        other = self.__polynomials(other)

        if not self or not other:
            return self.__polynomials.zero()
        else:
            r = (self.degree() + other.degree()).bit_length()
            m = 1 << r
            complexity_convolution = 1 << (2*r)
            complexity_fpm = 9/2 * m * r + 5 * m
            if complexity_fpm > complexity_convolution:
                zero = self.__polynomials.field.zero()
                result = [zero] * (self.degree() + other.degree() + 2)

                for i, x in enumerate(self.__coefficients):
                    for j, y in enumerate(other.coefficients()):
                        result[i + j] += x * y

                return self.__polynomials(result)
            else:
                p = [int(c) for c in self.coefficients()]
                q = [int(c) for c in other.coefficients()]

                res = FPM(p, q)
                return self.__polynomials(res)

    def __divmod__(self, other):
        other = self.__polynomials(other)

        # Lists will be modified, so copy them
        if not other:
            raise ZeroDivisionError
        else:
            dividend = self.coefficients()
            divisor = other.coefficients()
            n = other.degree()

            zero = self.__polynomials.zero()
            quotient = [zero] * (self.degree() - n + 1)

            for k in reversed(range(0, len(quotient))):
                quotient[k] = dividend[n + k] / divisor[n]
                for j in range(k, n + k):
                    dividend[j] -= quotient[k] * divisor[j - k]

            remainder = dividend[0:n]

            return self.__polynomials(quotient), self.__polynomials(remainder)

    def __call__(self, point):
        """
        Return the polynomial function's value at @p point.
        """
        value = self.__polynomials.field.zero()
        point = self.__polynomials.field(point)

        for c in reversed(self.__coefficients):
            value *= point
            value += c

        return value

    def ring(self):
        return self.__polynomials