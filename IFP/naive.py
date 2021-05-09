#!/usr/bin/python3
# -*- coding: utf-8 -*-


def factor_naive(n):
    """Naive factorization of n"""
    factors = []

    for factor in range(2, n // 2):
        q, r = divmod(n, factor)
        power = 0
        while r == 0:
            power += 1
            n = q
            q, r = divmod(q, factor)
        if power != 0:
            factors.append((factor, power))

    if factors == []:
        factors = [(n, 1)]

    return factors