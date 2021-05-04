#!/usr/bin/sage

from ..maths.math_lib import CRT
from .babystep_giantstep import BSGS
from sage.all import factor
from ..maths.groups import AbelianGroup


def PH(g, h, G:AbelianGroup):
    r"""Pohlig-Hellman algorithm for computing discrete logarithms
    in a finite abelian group `G` whose order is a smooth integer.

        INPUT:
        - `g` must be a primitive element of `G`
        - `h` must be an element of `G`

        OUTPUT:
        - `x` such that ``g^x = h``
    """

    decomposition = factor(G.order)
    residues = []
    moduli = []
    beta = G.inv(g)

    for pi, ni in decomposition:
        # compute the residue of x mod pi^ni
        n = G.order // pi
        yi = G.exp(g, n)
        zk = h
        bk = beta
        power = 1
        residue = 0

        for k in range(ni):
            # compute the k-th digit in base pi of x mod pi^ni
            wk = G.exp(zk, n)
            H = AbelianGroup(G.op, G.id, G.inv, pi)
            dk = BSGS(wk, yi, H)
            residue += dk * power
            zk = G.op(zk, G.exp(bk, dk))
            power *= pi
            bk = G.exp(bk, pi)
            n = n // pi

        residues.append(residue)
        moduli.append(power)

    x = CRT(residues, moduli)
    return x
