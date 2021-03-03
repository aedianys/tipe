#!/usr/bin/python3

from .ecc import EllipticCurve, Point, Inf

class AbelianGroup(object):
    def __init__(self, op, id, inv, ord):
        self.ord = ord
        self.id = id
        self.inv = inv
        self.op = op

    def exp(x,n):
        y = identity
        while n != 0:
            if n % 2 == 1:
                y = op(y,x)
            x = op(x,x)
            n //= 2
        return y


def ZpZ(p):
    def op(x,y): return (x*y) % p
    return AbelianGroup(op, 1, mod_inverse, p-1)

def ECCGroup(E,G=None):
    def inv(P): return -P
    def op(P,Q): return P+Q

    if G == None:
        return AbelianGroup(op, Inf(E), inv, E.order())
    else:
        return AbelianGroup(op, Inf(E), inv, G.order())
