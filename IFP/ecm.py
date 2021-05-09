#!/usr/bin/python3
# -*- coding: utf-8 -*-

from groups.elliptic_curve import EllipticCurve
from maths.math_lib import gcd
from maths.primes import primes_range
from math import log
from IFP.primepower import isprimepower
from primality.prime import is_prime
from random import randrange


def randomCurve(n):
    """Returns a random curve on field `GF(n)` and a random point on it."""
    x, y = randrange(1, n), randrange(1, n)
    a = randrange(1, n)
    b = (y ** 2 - x ** 3 - a * x) % n
    E = EllipticCurve(a, b, n)
    point = E(x, y)
    return E, point


def ECMFactor(n, sieve, limit, times):
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
