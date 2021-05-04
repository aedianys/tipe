#!/usr/bin/sage

from math import log
from sage.all import factor, next_prime
from ..maths.math_lib import inverse_mod, primes, reduced_row_echelon_form, array_is_inde,


def is_Bsmooth(B:int, n:int):
    """Determines if `n` is `B`-smooth and returns its factorization if so."""
    P = list(factor(n))
    if len(P) != 0 and P[-1] <= B: 
        return True, P
    else: return False, P


def index_calculus(g:int, h:int, p:int, r:int):
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
    sieve = primes(int(r*log(r)))
    while len(sieve) < r:
        sieve.append(next_prime(sieve[-1]))

    relations = []
    
    for k in range(p):
        b,P = is_Bsmooth(r, pow(g, k, p))
        if b:
            rel = [0] * (r+1)
            for x in P: rel[x] += 1
            rel[-1] = k

            if all([array_is_inde(rel,l) for l in relations]):
                relations.append(rel)
                if len(relations) >= r:
                    break
    
    add = lambda x,y : (x+y) % p
    mult = lambda x,y : (x*y) % p
    inv = lambda x : inverse_mod(x,p)

    logs = reduced_row_echelon_form(relations, inv, add, mult)
    s = 0
    y = h
    
    while True:
        b,P = is_Bsmooth(r, y)
        if b:
            fact = [0] * r
            for x in P: fact[x] += 1
            l = [fact[i]*logs[i][-1] for i in range(r)]
            return sum(l) - s
            
        y = (y*g) % p
        s += 1