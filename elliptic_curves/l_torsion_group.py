#!/usr/bin/python3
# -*- coding: utf-8 -*-

from structures.rings.quotients import QuotientRing
from structures.fields.fraction import FractionField
from elliptic_curves.elliptic_curve import EllipticCurve
from elliptic_curves.curve_polynomials import CurvePolynomials
from elliptic_curves.division_polynomials import DivisionPolynomialsList


class LTorsionGroups:
    """
    List of all l-torsion groups of an elliptic curve for odd l.
    """

    def __init__(self, curve):
        self.curve = curve
        self.curve_polynomials = CurvePolynomials(self.curve)
        self.division_polynomials = DivisionPolynomialsList(self.curve_polynomials)

    def __eq__(self, other):
        return self.curve == other.curve

    def __getitem__(self, index):
        return LTorsionGroup(self, index)


class LTorsionGroup:
    def __init__(self, torsion_groups, torsion):
        torsion = int(torsion)
        self.torsion_groups = torsion_groups
        if torsion < 3 or torsion % 2 == 0 or torsion_groups.curve.field(torsion) == 0:
            raise ValueError("only odd torsions greater than 1 are supported")

        self.__torsion = torsion
        self.__init_point()

    def representative(self):
        return self.__point

    def torsion(self):
        return self.__torsion

    def __init_point(self):
        """
        Create the representative of the l-torsion group.
        """
        R = self.torsion_groups.curve_polynomials
        psi = self.torsion_groups.division_polynomials[self.__torsion]

        S = QuotientRing(R, psi)
        T = FractionField(S)
        A, B = self.torsion_groups.curve.parameters()

        # Polynomials x and y on the curve interpreted
        # as elements of the field of fractions
        x = T(R.x())
        y = T(R.y())

        self.__point = EllipticCurve(T, A, B)(x, y)
