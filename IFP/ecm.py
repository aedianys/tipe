#!/usr/bin/python3
# -*- coding: utf-8 -*-

from elliptic_curves.elliptic_curve import EllipticCurve
from structures.fields.finite import FiniteField
from maths.math_lib import gcd
from maths.primes import primes_range
from math import log
from IFP.primepower import isprimepower
from primality.prime import is_prime
from random import randrange


def randomCurve(p):
    """Returns a random curve on field `GF(p)` and a random point on it."""
    x, y = randrange(1, p), randrange(1, p)
    a = randrange(1, p)
    b = (y ** 2 - x ** 3 - a * x) % p
    field = FiniteField(p)
    curve = EllipticCurve(field, a, b)
    point = curve(x, y)
    return curve, point


def ECMFactor(n, sieve, limit, times):
    """Finds a non-trivial factor of `n` using
    the elliptic-curve factorization method"""
    for _ in range(times):
        g = n
        while g == n:
            curve, point = randomCurve(n)
            g = gcd(n, curve.discriminant().remainder())
            
        if g > 1:
            return g

        for prime in sieve:
            prod = prime
            i = 0
            while prod < limit and i < 10:
                i += 1
                try:
                    point = prime * point
                except ZeroDivisionError as e:
                    return gcd(e.args[1], n)
                prod *= prime

    return n


def ECM(n, times=5):
    """Returns the list of prime factors of `n`
    using ECM factorization algorithm."""
    factors = []

    k = int(13 * log(n) ** 0.42)
    bound = max(10 * k ** 2, 100)
    sieve = primes_range(2, bound)

    for prime in sieve:
        q, r = divmod(n, prime)
        while r == 0:
            n = q
            factors.append(prime)
            q, r = divmod(n, prime)

    root, power = isprimepower(n)
    if power != 1:
        return [root] * power

    if is_prime(n):
        return [n]

    while n != 1:
        factor = ECMFactor(n, sieve, bound, times)
        factors += ECM(factor)
        n //= factor

    return sorted(factors)
