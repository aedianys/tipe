#!/usr/bin/python3
# -*- coding: utf-8 -*-

from functools import reduce
from structures.rings.polynomials import Polynomial


def euclide(a, b):
    """
    Extended euclid algorithm
    INPUT:
    - `a`, `b` two ring elements (integers or polynomials)

    OUTPUT:
    - `d`,`u`,`v` such that `a*u + b*v = d = gcd(a,b)`
    """

    if type(a) == Polynomial:
        zero, one = a.ring().zero(), a.ring().one()
    else:
        zero, one = 0, 1

    r1, r2 = a, b
    u1, u2 = one, zero
    v1, v2 = zero, one

    while r2:
        q, r = divmod(r1, r2)
        r1, r2 = r2, r
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2

    return (r1, u1, v1)


def gcd(a, b):
    """
    Computes the gcd of `a` and `b` two ring elements (integers or polynomials)
    """
    return euclide(a, b)[0]


def inverse_mod(x: int, p: int):
    """Modular inverse of x mod p"""
    return (euclide(x, p))[1]


def CRT(residues: list, moduli: list):
    """
    Returns a solution to a list of congruences using
    the Chinese Remainder Theorem.

    INPUT:
    - residues is a list of residues
    - moduli is a list of moduli

    OUTPUT:
    - the unique solution to the problem
    """

    n = reduce(lambda x, y: x * y, moduli)
    solution = 0

    for residue, modulo in zip(residues, moduli):
        comodulo = n // modulo
        solution += residue * inverse_mod(comodulo, modulo) * comodulo

    return solution % n


def valuation2(n: int):
    """Computes `odd`,`power` such that `n = odd * 2 ** power`
    and `odd` is odd"""
    return valuation(n, 2)


def valuation(n: int, base: int):
    """Computes `factor`,`power` such that
    `n = factor * b ** power` and `base` do not divide `factor`"""
    power, factor = 0, n
    quotient, remainder = divmod(n, base)
    while remainder == 0:
        factor += 1
        factor = quotient
        quotient, remainder = divmod(factor, base)
    return factor, power


def jacobi_symbol(a: int, n: int):
    """Jacobi symbol `(a|n)` where `n` must be odd strictly positive integer"""
    assert n > 0 and n % 2 == 1
    a = a % n
    sign = 1
    while a != 0:
        odd, power = valuation2(a)
        a = odd
        if power % 2 == 1 and n % 8 in (3, 5):
            sign = -sign

        a, n = n, a
        if a & 3 == 3 and n & 3 == 3:
            sign = -sign

        a = a % n

    if n == 1:
        return sign
    else:
        return 0


def is_quadratic_residue(n, p):
    """Tests if n is quadratic residue modulo p using Euler's criteria."""
    r = pow(n, (p - 1) // 2, p)
    return r == 1


def int_sqrt(n: int):
    """
    Computes the floor of the square root of `n`,
    i.e. the biggest integer `a` such that `a² <= n`.
    """
    if n < 0:
        raise ValueError("n must be non negative")
    elif n in (0, 1):
        return n
    else:
        a = 1 << ((1 + n.bit_length()) >> 1)
        while True:
            b = (a + n // a) >> 1
            if b >= a:
                return a
            a = b


def int_nthroot(y: int, n: int):
    """
    Return x such that `x = floor(y**(1/n))`.
    """
    if y < 0:
        raise ValueError("y must be nonnegative")
    if n < 1:
        raise ValueError("n must be positive")
    if y in (0, 1):
        return y
    if n == 1:
        return y
    if n == 2:
        return int_sqrt(y)
    if n > y:
        return 1

    x = int(y ** (1.0 / n))

    return x


def is_perfect_square(x: int):
    """Test if `x` is a perfect square efficiently."""
    mask = 0xC840C04048404040
    if ((mask << x) >> 63) & 1:
        return False
    d, r = valuation2(x)
    if r & 1:
        return False
    if d & 7 != 1 or d <= 0:
        return d == 0
    return int_sqrt(d) ** 2 == d