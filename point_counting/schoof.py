#!/usr/bin/python3
# -*- coding: utf-8 -*-

from math import ceil, sqrt
from maths.math_lib import gcd
from maths.primes import inverse_primorial
from maths.tonelli_shanks import square_root_mod
from maths.quotients import representative_in_range, solve_congruences
from structures.fields.finite import FiniteField
from structures.rings.polynomials import Polynomials
from elliptic_curves.l_torsion_group import LTorsionGroups


def frobenius(point, q):
    """Frobenius endomorphism"""
    return point.curve()(point.x ** q, point.y ** q)


def frobenius_trace_mod_2(curve):
    """Computes the trace of the Frobenius endomorphism mod 2"""
    R = Polynomials(curve.field)

    x = R.x()
    A, B = curve.parameters()

    curve_polynomial = x ** 3 + A * x + B
    field_polynomial = R([0] * curve.field.size() + [1]) - x

    if gcd(curve_polynomial, field_polynomial).degree() == 0:
        return FiniteField(2)(1)
    else:
        return FiniteField(2)(0)


def frobenius_trace_mod_l(torsion_group):
    l = torsion_group.torsion()
    torsion_quotient_ring = FiniteField(l)
    field_size = torsion_group.torsion_groups.curve.field.size()
    psi = torsion_group.torsion_groups.division_polynomials

    point = torsion_group.representative()
    frobenius_point = frobenius(point, field_size)
    frobenius2_point = frobenius(frobenius_point, field_size)

    valid_range = range(-(l - 1) // 2, (l + 1) // 2)
    factor = representative_in_range(torsion_quotient_ring(field_size), valid_range)

    curve = point.curve()
    # determinant_point = factor * point
    # alternative way of computing point multipication with division polynomials
    def multiplication(p, n):
        sign = 1 if n > 0 else -1
        k = abs(n)
        x = p.x - curve.field((psi[k - 1] * psi[k + 1]), (psi[k] ** 2))
        y = curve.field(sign * psi[2 * k], 2 * psi[k] ** 4)
        return curve(x, y)

    determinant_point = multiplication(point, factor)

    point_sum = frobenius2_point + determinant_point
    if not point_sum:
        return torsion_quotient_ring(0)

    elif not (frobenius2_point - determinant_point):
        w = square_root_mod(field_size, l)
        # root_point = w * frobenius_point
        root_point = multiplication(frobenius_point, w)
        if determinant_point.y == root_point.y:
            return torsion_quotient_ring(2 * w)
        else:
            return torsion_quotient_ring(-2 * w)

    trace_point = frobenius_point
    for trace_candidate in range(1, (l + 1) // 2):
        if point_sum.x == trace_point.x:
            if point_sum.y == trace_point.y:
                return torsion_quotient_ring(trace_candidate)
            else:
                return torsion_quotient_ring(-trace_candidate)
        else:
            trace_point += frobenius_point


def hasse_trace_range(field):
    """
    The interval given by Hasse's theorem in which the
    trace of the Frobenius endomorphism must be.
    """
    l = 2 * ceil(sqrt(field.size()))
    return range(-l, l + 1)


def frobenius_trace(curve):
    """
    Computes the trace of the Frobenius endomorphism
    for the given curve.
    """

    trace_congruences = []
    search_range = hasse_trace_range(curve.field)
    torsion_primes = inverse_primorial(len(search_range), curve.field.characteristic())

    if 2 in torsion_primes:
        trace_congruences.append(frobenius_trace_mod_2(curve))
        torsion_primes.remove(2)

    torsion_group = LTorsionGroups(curve)
    for prime in torsion_primes:
        trace_congruences.append(frobenius_trace_mod_l(torsion_group[prime]))

    trace_residue = solve_congruences(trace_congruences)
    return representative_in_range(trace_residue, search_range)


def schoof(curve):
    """
    Return the order of the curve using schoof's algorithm
    """
    p = curve.field.characteristic()
    return p + 1 - frobenius_trace(curve)