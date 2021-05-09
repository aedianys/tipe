#!/usr/bin/python3
# -*- coding: utf-8 -*-

from groups import AbelianGroup, AbelianGroupElement


class EllipticCurve(AbelianGroup):
    """Represents a single elliptic curve defined over a finite field.
    p must be prime, since we use the modular inverse to compute point
    addition.
    """

    def __init__(self, field, A, B):
        self.field = field
        self.A = field(A)
        self.B = field(B)

    def __str__(self):
        return f"y^2 = x^3 + {self.A}x + {self.B}"

    def __eq__(self, other):
        return (self.A, self.B) == (other.A, other.B)

    def parameters(self):
        return self.A, self.B

    def __call__(self, x, y, z=None):
        if z is None:
            return Point(self, x, y, 1)
        else:
            return Point(self, x, y, z)

    def id(self):
        return Point(self, 0, 1, 0)

    def inf(self):
        return self.id()

    def eval(self, x):
        A, B = self.parameters()
        return x ** 3 + A * x + B

    def has_point(self, P):
        if not P:
            return True
        else:
            return P.y ** 2 == self.eval(P.x)

    def is_singular(self):
        return 4 * (self.A) ** 3 + 27 * (self.B) ** 2 == 0


class Point(AbelianGroupElement):
    """A point on a specific curve."""

    def __init__(self, curve, x, y, z):
        self.__curve = curve
        self.x = self.__curve.field(x)
        self.y = self.__curve.field(y)
        self.z = z

    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __getitem__(self, index):
        return [self.x, self.y, self.z][index]

    def __bool__(self):
        return bool(self.z)

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __neg__(self):
        return self.__curve(self.x, -self.y, self.z)

    def __add__(self, other):
        assert self.__curve == other.curve()

        if not self:
            return other

        elif not other:
            return self

        elif self == -other:
            return self.__curve.inf()

        else:
            if self.x == other.x:
                if not self.y:
                    return self.__curve.inf()  # only if singular
                else:
                    return self.__double__()
            else:
                return self.__generic_add__(other)

    def __double__(self):
        A, B = self.__curve.parameters()
        delta = (3 * self.x ** 2 + A) / (2 * self.y)
        x = delta ** 2 - (2 * self.x)
        y = delta * (self.x - x) - self.y
        return self.__curve(x, y)

    def __generic_add__(self, other):
        delta = (other.y - self.y) / (other.x - self.x)
        x = delta ** 2 - self.x - other.x
        y = delta * (self.x - x) - self.y
        return self.__curve(x, y)

    def group(self):
        return self.__curve

    def curve(self):
        return self.group()