#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import is_quadratic_residue


def naive_order(curve):
    """Naive point counting for curve over Z/pZ"""
    count = 1
    p = curve.field.characteristic()

    for x in range(p):
        y2 = curve.eval(x)
        if not y2:
            count += 1
        elif is_quadratic_residue(y2.remainder(), p):
            count += 2

    return count