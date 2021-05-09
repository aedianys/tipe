#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import gcd
from primality.prime import is_prime


def pollard_rho_search_floyd(f, n: int, x0=2):
    """Pollard's-rho algorithm to find a non-trivial
    factor of `n` using Floyd's cycle detection algorithm."""
    d = 1
    tortoise = x0
    hare = x0

    while d == 1:
        tortoise = f(tortoise)
        hare = f(f(hare))
        d = gcd(abs(tortoise - hare), n)

    return d


def pollard_rho_search_brent(f, n: int, x0=2):
    """Pollard's-rho algorithm to find a non-trivial
    factor of `n` using Brent's cycle detection algorithm."""
    power = 1
    lam = 0
    tortoise = x0
    hare = x0
    d = 1
    while d == 1:
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        hare = f(hare)
        d = gcd(abs(tortoise - hare), n)
        lam += 1

    return d


def pollard_brent_search(f, n: int, x0=2, size=100):
    """Pollard's-rho algorithm to find a non-trivial
    factor of `n` using Brent's cycle detection algorithm,
    and including the product speed improvment."""
    power = 1
    lam = 0
    tortoise = x0
    hare = x0
    d = 1

    while d == 1:
        prod = 1
        last_pos = hare
        for c in range(size):
            if power == lam:
                tortoise = hare
                power *= 2
                lam = 0
            hare = f(hare)
            lam += 1
            prod = (prod * abs(tortoise - hare)) % n
        d = gcd(prod, n)

    if not (is_prime(d)):
        return pollard_rho_search_brent(f, n, last_pos)
    else:
        return d


def pollard_rho(n: int):
    """Pollard's-rho algorithm to find the factorization of `n`
    using Brent's cycle detection algorithm, and including
    the product speed improvment."""
    f = lambda x: (x ** 2 + 1) % n
    factors = []

    while not (is_prime(n)):
        x0 = 2
        factor = pollard_brent_search(f, n)
        while factor == n:
            x0 += 1
            factor = pollard_brent_search(f, n)
        factors.append(factor)
        n //= factor

    if n != 1:
        factors.append(n)

    return factors