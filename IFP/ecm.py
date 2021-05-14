#!/usr/bin/python3
# -*- coding: utf-8 -*-

from elliptic_curves.elliptic_curve import EllipticCurve
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
    E = EllipticCurve(a, b, p)
    point = E(x, y)
    return E, point


def ECMFactor(n, sieve, limit, times):
    """Finds a prime factor of `n` using
    the elliptic-curve factorization method"""
    for _ in range(times):
        g = n
        while g == n:
            E, point = randomCurve(n)
            A, B = E.parameters()
            g = gcd(4 * A ** 3 + 27 * B ** 2, n)
        if g > 1:
            return g

        for prime in sieve:
            prod = prime
            while prod < limit:
                try:
                    point = prime * point
                except ZeroDivisionError as e:
                    return gcd(e, n)
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
        factors.append(factor)
        n //= factor

    return sorted(factors)
