#!/usr/bin/python3

from .math_lib import inverse_mod

class EllipticCurve(object):
    """Represents a single elliptic curve defined over a finite field.

    See here:
        http://en.wikipedia.org/wiki/Elliptic_curve
        http://jeremykun.com/2014/02/24/elliptic-curves-as-python-objects/

    p must be prime, since we use the modular inverse to compute point
    addition.

    """

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def eval(self, x):
        return (x**3 + self.a * x + self.b) % self.p

    def __eq__(self, C):
        return (self.a, self.b) == (C.a, C.b)

    def has_point(self, x, y):
        return (y ** 2) % self.p == (x ** 3 + self.a * x + self.b) % self.p

    def __str__(self):
        return f"y^2 = x^3 + {self.a}x + {self.b}"


class Point(object):
    """A point on a specific curve."""

    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p

        if not self.curve.has_point(x, y):
            raise ValueError(f"{self} is not on curve {self.curve}")

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __getitem__(self, index):
        return [self.x, self.y][index]

    def __eq__(self, Q):
        return (self.curve, self.x, self.y) == (Q.curve, Q.x, Q.y)

    def __neg__(self):
        return Point(self.curve, self.x, -self.y)

    def __add__(self, Q):
        """Add two points together.

        We need to take care of special cases:
         * Q is the infinity point (0)
         * P == Q
         * The line crossing P and Q is vertical.

        """
        assert self.curve == Q.curve

        # 0 + P = P
        if isinstance(Q, Inf):
            return self

        xp, yp, xq, yq = self.x, self.y, Q.x, Q.y
        p = self.curve.p
        m = None

        # P == Q
        R = 0
        if self == Q:
            if yp == 0:
                R = Inf(self.curve)
            else:
                m = (
                    (3 * xp * xp + self.curve.a) * mod_inverse(2 * yp, p)
                ) % p

        # Vertical line
        elif xp == xq:
            R = Inf(self.curve)

        # Common case
        else:
            m = ((yq - yp) * mod_inverse(xq - xp, p)) % p

        if m is not None:
            xr = (m ** 2 - xp - xq) % p
            yr = (m * (xp - xr) - yp) % p
            R = Point(self.curve, xr, yr)

        return R

    def __rmul__(self, n):
        Q = self
        R = Inf(self.curve)

        while n > 0:
            if n % 2 == 1:
                R = R + Q
            Q = Q + Q
            n = n // 2
        return R

    def __mul__(self, n):
        return self.__rmul__(n)


class Inf(Point):
    """The custom infinity point."""

    def __init__(self, curve):
        self.curve = curve

    def __eq__(self, Q):
        return isinstance(Q, Inf)

    def __neg__(self):
        """-0 = 0"""
        return self

    def __add__(self, Q):
        """P + 0 = P"""
        return Q
