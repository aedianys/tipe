#!/usr/bin/python3

from math import ceil, sqrt

def BSGS(a, b, G):
    r"""Baby-step Giant-step algorithm for computing discrete logarithms
    in a finite abelian group `G`.

        INPUT:
        - `a` must be a primitive element of `G`
        - `b` must be an element of `G`

        OUTPUT:
        - `x` such that ``a^x = b``
    """

    dic = {}
    m = ceil(sqrt(G.order))

    aj = G.id
    for j in range(m):
        aj = G.op(aj,a)
        dic[aj] = j
    
    c = G.exp(G.inv(a),m)
    gam = b
    for i in range(m):
        try:
            _ = dic[gam]
            return m*i+j
        except:
            gam = G.op(gam,c)