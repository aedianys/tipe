#!/usr/bin/python3

from ..maths.math_lib import int_nthroot


def isprimepower(n:int):
    x = n
    m = 1
    while x >= 2:
        m += 1
        x = int_nthroot(n,m)
        if x**m == n:
            return x,m
    return n,1