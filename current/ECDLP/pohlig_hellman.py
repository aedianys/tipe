#!/usr/bin/sage

from ..maths.groups import ECCGroup
from ..DLP.pohlig_hellman import PH


def ECC_PH(G, A, E):
    r"""Pohlig-Hellman algorithm for computing discrete logarithms
    over an elliptic curve `E` whose order is a smooth integer.

    INPUT:
    - `G` a point of `E`
    - `A` a point of the group generate by `G`

    OUTPUT:
    - `n` such that `A = n*G`
    """
    return PH(G, A, ECCGroup(E, G))
