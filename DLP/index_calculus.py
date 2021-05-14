#!/usr/bin/python3
# -*- coding: utf-8 -*-

from IFP.factor import factor
from maths.gaussian_elimination import reduced_row_echelon_form, vectors_are_inde
from structures.fields.finite import FiniteField


def is_Bsmooth(B: int, n: int):
    """Determines if `n` is `B`-smooth and returns its factorization if so."""
    factors = list(factor(n))
    if len(factors) != 0 and factors[-1] <= B:
        return True, factors
    else:
        return False, factors


def index_calculus(g: int, h: int, p: int, r: int):
    """Computes discrete logarithm of `h` in `(Z/pZ)*` in base `g`
    using the index calculus probabilistic algorithm

        INPUT:
        - `g` must be a generator of Z/pZ
        - `h` an element of Z/pZ
        - `p` the modulus
        - `r` the number of maximum equations to solve

        OUTPUT:
        - `x` such that ``g^x = h``

    """
    relations = []
    G = FiniteField(p)
    g, h = G(g), G(h)

    for k in range(p):
        b, factors = is_Bsmooth(r, pow(g, k, p))
        if b:
            rel = [G.zero()] * (r + 1)
            for f in factors:
                rel[f] += 1
            rel[-1] = k

            if all([vectors_are_inde(rel, l) for l in relations]):
                relations.append(rel)
                if len(relations) >= r:
                    break

    logs = reduced_row_echelon_form(relations)
    s = 0
    y = h

    while True:
        b, factors = is_Bsmooth(r, y)
        if b:
            fact = [0] * r
            for f in factors:
                fact[f] += 1
            l = [fact[i] * logs[i][-1] for i in range(r)]
            x = sum(l) - s
            return x.remainder()

        y *= g
        s += 1