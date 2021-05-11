#!/usr/bin/python3
# -*- coding: utf-8 -*-

from structures.fields.finite import FiniteField


def pollard_rho(a, b, group, hash):
    """Pollard's-rho algorithm to solve discrete logarithm
    factor of `n` using Floyd's cycle detection algorithm.

        INPUT:
        - `a` must be a primitive element of `group`
        - `b` must be an element of `group`

        OUTPUT:
        - `x` such that ``a^x = b``
    """

    n = group.order()
    F = FiniteField(n)
    zero = F.zero()

    def f(t):
        x, ai, bi = t
        y = hash(x)
        if y % 3 == 0:
            return (x * b, ai, bi + 1)
        elif y % 3 == 1:
            return (x ** 2, 2 * ai, 2 * bi)
        else:
            return (a * x, ai + 1, bi)

    tortoise = f(group.id(), zero, zero)
    hare = f(f(group.id(), zero, zero))

    while tortoise[0] != hare[0]:
        tortoise = f(tortoise)
        hare = f(f(hare))

    return (hare[1] - tortoise[1]) / (tortoise[2] - hare[2])