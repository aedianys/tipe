#!/usr/bin/python3

from current.maths.math_lib import inverse_mod


def pollard_rho(a, b, G, hash):
    """Pollard's-rho algorithm to solve discrete logarithm
    factor of `n` using Floyd's cycle detection algorithm.

        INPUT:
        - `a` must be a primitive element of `G`
        - `b` must be an element of `G`

        OUTPUT:
        - `x` such that ``a^x = b``
    """

    n = G.order

    def f(t):
        x,ai,bi = t
        y = hash(x)
        if y % 3 == 0: return (G.op(x,b), ai, (bi + 1) % n)
        elif y % 3 == 1: return (G.op(x,x), (2 * ai) % n, (2 * bi) % n)
        else : return (G.op(a,x), (ai + 1) % n, bi)

    tortoise = f(G.id,0,0)
    hare = f(f(G.id,0,0))

    while tortoise[0] != hare[0]:
        tortoise = f(tortoise)
        hare = f(f(hare))
    
    r = tortoise[2] - hare[2]
    return (inverse_mod(r) * (hare[1] - tortoise[1])) % n