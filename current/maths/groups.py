#!/usr/bin/python3

from .math_lib import inverse_mod
from .ecc import Inf


class AbelianGroup(object):
    def __init__(self, op, id, inv, ord):
        self.ord = ord
        self.id = id
        self.inv = inv
        self.op = op

    def exp(self, x, n):
        y = self.id
        while n != 0:
            if n % 2 == 1:
                y = self.op(y, x)
            x = self.op(x, x)
            n //= 2
        return y


def ZpZ_mult(p):
    def op(x, y):
        return (x * y) % p
    
    def inv(x):
        return inverse_mod(x,p)

    return AbelianGroup(op, 1, inv, p - 1)


def ECCGroup(E, G=None):
    def inv(P):
        return -P

    def op(P, Q):
        return P + Q

    if G is None:
        return AbelianGroup(op, Inf(E), inv, E.order())
    else:
        return AbelianGroup(op, Inf(E), inv, G.order())
