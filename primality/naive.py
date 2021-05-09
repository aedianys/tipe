#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import int_sqrt


def is_prime_naive(n):
    if n < 5:
        return n == 2 or n == 3
    elif n % 2 == 0:
        return False
    else:
        bound = int_sqrt(n)
        factor = 3
        while factor <= bound:
            if n % factor == 0:
                return False
            factor += 2
        return True