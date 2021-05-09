#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import int_nthroot


def isprimepower(n: int):
    x = n
    power = 1

    while x >= 2:
        power += 1
        x = int_nthroot(n, power)
        if x ** power == n:
            return x, power

    return n, 1