#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import valuation2

STEPS = [
    0,
    341531,
    1050535501,
    350269456337,
    55245642489451,
    7999252175582851,
    585226005592931977,
    0x10000000000000000,
]

BASES = [
    [9345883071009581737],
    [336781006125, 9639812373923155],
    [4230279247111683200, 14694767155120705706, 16641139526367750375],
    [2, 141889084524735, 1199124725622454117, 11096072698276303650],
    [2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805],
    [
        2,
        123635709730000,
        9233062284813009,
        43835965440333360,
        761179012939631437,
        1263739024124850375,
    ],
    [2, 325, 9375, 28178, 450775, 9780504, 1795265022],
]


def is_pseudoprime_MillerRabin_base(n: int, a: int):
    """Test if `n` is a strong probable prime in base `a`"""
    d, s = valuation2(n - 1)
    x = pow(a, d, n)

    if x == 1 or x == n - 1:
        return True

    for i in range(1, s):
        x = x * x % n
        if x == n - 1:
            return True

    return False


def is_prime_MillerRabin(n: int):
    if n >= STEPS[-1]:
        raise NotImplementedError

    for i in range(6):
        if STEPS[i] <= n < STEPS[i + 1]:
            base = BASES[i]

    return all(is_prime_MillerRabin(n, a) for a in base)
