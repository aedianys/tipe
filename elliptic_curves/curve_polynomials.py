#!/usr/bin/python3
# -*- coding: utf-8 -*-

from structures.rings import CommutativeRing, CommutativeRingElement
from structures.rings.polynomials import Polynomials


class CurvePolynomials(CommutativeRing):
    """
    A ring of polynomials over an elliptic curve. It is in fact the quotient ring
    F[x,y]/(y^2 - x^3 - Ax - B) where F is the field used.
    """

    def __init__(self, curve):
        self.curve = curve
        self.field = curve.field
        self.polynomials = Polynomials(self.field)
        self.y2_reduction = self.polynomials([self.curve.B, self.curve.A, 0, 1])

    def __eq__(self, other):
        return self.curve, self.field == other.curve, other.field

    def __call__(self, x_factor, y_factor=None):
        return CurvePolynomial(self, x_factor, y_factor)

    def zero(self):
        return self(0, 0)

    def one(self):
        return self(1, 0)

    def x(self):
        return self((0, 1), 0)

    def y(self):
        return self(0, 1)


class CurvePolynomial(CommutativeRingElement):
    """
    A polynom over an elliptic curve, since this is in fact a quotient class,
    we can describe it with a polynomial whose degree in y is one or less.
    """

    def __init__(self, curve_polynomials, x_factor, y_factor=None):
        """
        Create a new polynomial from the given coefficient lists x_factor
        and y_factor.
        x_factor: polynomial in x only
        y_factor: polynomial in x (implicitly * y)
        """

        self.__curve_polynomials = curve_polynomials

        if isinstance(x_factor, self.__class__):
            self.__x_factor = x_factor.__x_factor
            self.__y_factor = x_factor.__y_factor
        else:
            self.__x_factor = self.__curve_polynomials.polynomials(x_factor)
            if y_factor is None:
                self.__y_factor = self.__curve_polynomials.polynomials.zero()
            else:
                self.__y_factor = self.__curve_polynomials.polynomials(y_factor)

    def x_factor(self):
        return self.__x_factor

    def y_factor(self):
        return self.__y_factor

    def __str__(self):
        if not self.__y_factor:
            return f"{self.__x_factor}"
        elif not self.__x_factor:
            if self.__y_factor.degree() == 0:
                return f"{self.__y_factor}y"
            else:
                return f"({self.__y_factor})y"
        else:
            if self.__y_factor.degree() == 0:
                return f"{self.__x_factor} + {self.__y_factor}y"
            else:
                return f"{self.__x_factor} + ({self.__y_factor})y"

    def __bool__(self):
        return bool(self.__x_factor) or bool(self.__y_factor)

    def __eq__(self, other):
        try:
            other = self.__curve_polynomials(other)
            return (
                self.__x_factor == other.x_factor()
                and self.__y_factor == other.y_factor()
            )
        except TypeError:
            return False

    def __add__(self, other):
        x = self.__x_factor + other.x_factor()
        y = self.__y_factor + other.y_factor()
        return self.__curve_polynomials(x, y)

    def __neg__(self):
        return self.__curve_polynomials(-self.__x_factor, -self.__y_factor)

    def __mul__(self, other):
        """
        Multiply two curve polynomials using the equality
        y^2 = x^3 + A * x^2 + B in order to get a degree in
        y of one or less.
        """
        other = self.__curve_polynomials(other)
        xx = self.__x_factor * other.x_factor()
        yy = self.__y_factor * other.y_factor()
        xy = (self.__x_factor + self.__y_factor) * (other.x_factor() + other.y_factor())
        y = xy - xx - yy
        x = xx + yy * self.__curve_polynomials.y2_reduction

        return self.__curve_polynomials(x, y)

    def __divmod__(self, other):
        """
        Supports only the following cases :
        - other is a polynomial in x alone
        - self and other are polynomials in y alone
        """
        other = self.__curve_polynomials(other)

        if not other:
            raise ZeroDivisionError

        elif not self:
            return self.__curve_polynomials.zero(), self.__curve_polynomials.zero()

        elif other.y_factor() and (self.__x_factor or other.x_factor()):
            raise NotImplementedError

        elif other.y_factor():
            qx, rx = divmod(self.__y_factor, other.y_factor())
            zero = self.polynomial_ring().zero()
            qy, ry = zero, zero
        else:
            qx, rx = divmod(self.__x_factor, other.x_factor())
            qy, ry = divmod(self.__y_factor, other.x_factor())

        q = self.__curve_polynomials(qx, qy)
        r = self.__curve_polynomials(rx, ry)

        return q, r

    def ring(self):
        return self.__curve_polynomials

    def polynomial_ring(self):
        return self.ring().polynomials