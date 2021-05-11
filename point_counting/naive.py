#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import jacobi_symbol


def naive_order(curve):
    """Naive point counting for curve over Z/pZ"""
    count = 1
    p = curve.field.characteristic()

    for x in range(p):
        y2 = curve.eval(x)
        j = jacobi_symbol(y2.remainder(), p)
        if j == 0:
            count += 1
        elif j == 1:
            count += 2

    return count