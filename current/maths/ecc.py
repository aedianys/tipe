#!/usr/bin/python3

from .math_lib import euclide


class CurveError(Exception):
    def __init__(self, factor, message=None) -> None:
        self.factor = factor
        self.message = message
        super().__init__(self.message)


class EllipticCurve(object):
    """Represents a single elliptic curve defined over a finite field.

    See here:
        http://en.wikipedia.org/wiki/Elliptic_curve

    p must be prime, since we use the modular inverse to compute point
    addition.

    """

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def eval(self, x):
        return (x ** 3 + self.a * x + self.b) % self.p

    def __eq__(self, C):
        return (self.a, self.b) == (C.a, C.b)

    def has_point(self, P):
        return (P.y ** 2 * P.z) % self.p == (
            P.x ** 3 + self.a * P.x * P.z ** 2 + self.b + P.z ** 2
        ) % self.p

    def __str__(self):
        return f"y^2 = x^3 + {self.a}x + {self.b}"

    def inf(self):
        return Point(self, 0, 1, 0)


class Point(object):
    """A point on a specific curve."""

    def __init__(self, curve, x, y, z):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        self.z = z % curve.p

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __getitem__(self, index):
        return [self.x, self.y, self.z][index]

    def __eq__(self, Q):
        return (self.curve, self.x, self.y, self.z) == (Q.curve, Q.x, Q.y, Q.z)

    def __neg__(self):
        return Point(self.curve, self.x, -self.y, self.z)

    def __add__(self, Q):
        assert self.curve == Q.curve

        if self.z == 0:
            return Q

        if Q.z == 0:
            return self

        xp, yp, xq, yq = self.x, self.y, Q.x, Q.y
        p = self.curve.p
        R = None

        # Vertical line
        if xp == xq and (yp + yq) % p == 0:
            R = self.curve.inf()

        # P == Q
        elif self == Q:
            if yp == 0:
                R = self.curve.inf()
            else:
                num = 3 * xp * xp + self.curve.a
                denom = 2 * yp

        else:
            num = yq - yp
            denom = xq - xp

        if R is None:
            d, inv, _ = euclide(denom, p)
            if d != 1:
                raise CurveError(d)
#                return Point(self.curve, 0, 0, d)
            else:
                m = inv * num
                xr = (m ** 2 - xp - xq) % p
                yr = (m * (xp - xr) - yp) % p
                R = Point(self.curve, xr, yr, 1)

        return R

    def __rmul__(self, n):
        Q = self
        R = self.curve.inf()

        while n > 0:
            if n & 1:
                R = R + Q
            Q = Q + Q
            n >>= 1
        return R

    def __mul__(self, n):
        return self.__rmul__(n)