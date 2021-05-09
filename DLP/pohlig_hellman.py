#!/usr/bin/python3
# -*- coding: utf-8 -*-

from maths.math_lib import CRT
from DLP.babystep_giantstep import BSGS
from sage.all import factor


def pohlig_hellman_prime(g, h, group, p, e):
    inv = g ** (-1)

    power = [1] * (e + 1)  # power[k] = p ** k
    for k in range(1, e):
        power[k] = power[k - 1] * p

    y = g ** power[e - 1]

    xk = 0  # the residue accumulation

    for k in range(e):
        # compute the k-th digit in base pi of x mod pi^ni
        hk = ((inv ** xk) * h) ** power[e - 1 - k]
        dk = BSGS(y, hk, group, p)
        xk = xk + dk * power[k]

    return xk


def pohlig_hellman(g, h, group, order=None):
    """Pohlig-Hellman algorithm for computing discrete logarithms
    in a finite abelian group `group` whose order is a smooth integer.

        INPUT:
        - `g` must be a primitive element of `group`
        - `h` must be an element of `group`

        OUTPUT:
        - `x` such that ``g^x = h``

    See https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
    for details.
    """

    if order is None:
        order = group.order()

    decomposition = factor(order)
    residues = []
    moduli = []

    for pi, ei in decomposition:
        # compute the residue of x mod pi^ni
        modulus = pi ** ei
        ni = order // modulus
        gi, hi = g ** ni, h ** ni
        residue = pohlig_hellman_prime(gi, hi, group, pi, ei)
        residues.append(residue)
        moduli.append(modulus)

    x = CRT(residues, moduli)
    return x