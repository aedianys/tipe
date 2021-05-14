#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import valuation2, is_quadratic_residue


def square_root_mod(n, p):
    """
    Computes the square root r of n modulo a prime p
    if it exists, ie such that r² = n (mod p)
    """
    assert is_quadratic_residue(n, p)

    q, s = valuation2(p - 1)

    z = 2
    while is_quadratic_residue(z, p):
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)

    while True:
        if t == 0:
            return 0
        elif t == 1:
            return r
        else:
            x = pow(t, 2, p)
            i = 1
            while x != 1:
                x = pow(x, 2, p)
                i += 1

            b = pow(c, 2 ** (m - i - 1), p)
            m = i
            c = pow(c, 2, p)
            t = (t * c) % p
            r = (r * b) % p
