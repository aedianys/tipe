#!/usr/bin/sage

from .math_lib import inverse_mod, CRT
from sage.all import factor, discrete_log
from groups import AbelianGroup, ECCGroup

def PH(a, y, G):
    r"""Pohlig-Hellman algorithm for computing discrete logarithms
    in a finite abelian group G whose order is a smooth integer.

        INPUT:
        `a` must be a primitive element of G
        y must be an element of G

        OUTPUT:
        x such that a^x = y
    """

    decomposition = factor(G.order)
    residues = []
    moduli = []
    beta = G.inv(a)

    for pi,ni in decomposition:
        #compute the residue of x mod pi^ni
        n = G.order//pi
        yi = G.exp(a,n)
        z = y
        b = beta
        power = 1
        residue = 0

        for j in range(ni):
            #compute the j-th digit in base pi of x mod pi^ni
            w = G.exp(z,n)
            bj = discrete_log(w,yi,ord = pi, identity = G.id, inverse=G.inv, op=G.op)
            residue += bj*power
            z = G.op(z, G.exp(b,bj))
            power *= pi
            b = G.exp(b,pi)
            n = n // pi

        residues.append(residue)
        moduli.append(power)

    x = CRT(residues, moduli)
    return x

def ECC_PH(G,A,E):
    r"""Pohlig-Hellman algorithm for computing discrete logarithms
        over an elliptic curve E whose order is a smooth integer.

        INPUT:
        G a point of E
        A a point of the group generate by G

        OUTPUT:
        n such that A = n*G
    """
    return PH(G, A, ECCGroup(E,G))
