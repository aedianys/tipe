#!/usr/bin/python3

from ..maths.ecc import CurveError, EllipticCurve, Point
from ..maths.math_lib import gcd, primes
from .primepower import isprimepower
from ..prime import is_prime
from random import randrange

def randomCurve(n):
    """Returns a random curve on field `GF(n)` and a random point on it."""
    x, y = randrange(1,n), randrange(1,n)
    a = randrange(1,n)
    b = (y ** 2 - x ** 3 - a * x) % n
    E = EllipticCurve(a, b, n)
    P = Point(E, x, y, 1)
    return E, P


def ECMFactor(n, ps, limit, times):
    for _ in range(times):
        g = n
        while g == n:
            E, P = randomCurve(n)
            g = gcd(4 * E.a ** 3 + 27 * E.b ** 2, n)
        if g > 1:
            return g

        for p in ps:
            prod = p
            while prod < limit:
                try:
                    P = p * P
                except CurveError as e:
                    return gcd(e.factor,n)
                prod *= p

    return n


def ECM(n, limit=1000, times=5):
    factors = []
    ps = primes(limit)

    for p in (2, 3):
        q, r = divmod(n, p)
        while r == 0:
            n = q
            factors.append(p)
            q, r = divmod(n, p)

    m,r = isprimepower(n)
    if r != 1:
        return [m]*r

    if is_prime(n):
        return [n]

    while n != 1:
        d = ECMFactor(n, ps, limit, times)
        factors.append(d)
        n //= d

    return sorted(factors)
