#!/usr/bin/python3

from ..maths.math_lib import int_sqrt


def isprime_naive(n):
    if n < 5:
        return n == 2 or n == 3
    elif n % 2 == 0:
        return False
    else:
        r = int_sqrt(n)
        k = 3
        while k <= r:
            if n % k == 0:
                return False
            k += 2
        return True