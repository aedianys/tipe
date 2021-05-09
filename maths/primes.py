#!/usr/bin/python3
# -*- coding: utf-8 -*-


def primes_range(lesser_bound, upper_bound):
    """List all primes between `lesser_bound` and `upper_bound`
    (included) using Eratosthenes' sieve"""
    tab = [True] * (upper_bound + 1)
    sieve = []
    for p in range(2, upper_bound + 1):
        if tab[p]:
            sieve.append(p)
            for i in range(p, upper_bound + 1, p):
                tab[i] = False

    return list(filter(lambda x: x >= lesser_bound, sieve))


def inverse_primorial(n, ignored=0):
    """
    Return a list of the first primes whose product is greater than, or equal
    to `n`, but do not use `ignored`.
    """
    primes = primes_range(2, n + 1)

    # Find the smallest product of primes that is at least n, but don't use
    # the ignored prime.
    product = 1
    i = 0
    for index, prime in enumerate(primes):
        if prime != ignored:
            product *= prime
            if product >= n:
                i = index
                break

    # Throw away excess primes
    primes = primes[: i + 1]
    if ignored in primes:
        primes.remove(ignored)

    # Try to remove unnecessary primes, largest first.
    for index, prime in enumerate(reversed(primes)):
        modified_product = product / prime
        if modified_product >= n:
            product = modified_product
            primes[-(index + 1)] = 0

    return list(filter(None, primes))