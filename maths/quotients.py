#!/usr/bin/python3
# -*- coding: utf-8 -*-

from functools import reduce
from maths.math_lib import euclide
from rings.quotients import QuotientRing


def solve_congruences(congruences):
    """
    Returns a solution to a list of congruences using
    the Chinese Remainder Theorem in form of quotients elements.
    """
    if not congruences:
        raise ValueError("cannot solve empty equation system")

    mul = lambda x, y: x * y
    common_ring = congruences[0].ring().ring
    common_modulus = reduce(mul, [c.modulus() for c in congruences])
    common_representative = 0

    for c in congruences:
        neutralizer = common_modulus // c.modulus()
        common_representative += (
            c.remainder() * neutralizer * euclide(neutralizer, c.modulus())[1]
        )

    quotient_ring = QuotientRing(common_ring, common_modulus)
    return quotient_ring(common_representative)


def representative_in_range(quotient_class, valid_range):
    """
    Return the unique representative of quotient_class in valid_range.
    """
    if len(valid_range) > quotient_class.modulus():
        raise ValueError("solution not unique")

    # range[0] is the first element, that is, the lower bound
    q, r = divmod(valid_range[0], quotient_class.modulus())
    shifted_range = range(r, r + len(valid_range))

    if quotient_class.remainder() in shifted_range:
        return q * quotient_class.modulus() + quotient_class.remainder()
    elif quotient_class.remainder() + quotient_class.modulus() in shifted_range:
        return (q + 1) * quotient_class.modulus() + quotient_class.remainder()
    else:
        raise ValueError("no solution")