#!/usr/bin/python3
# -*- coding: utf-8 -*-

from IFP.factor import factor

def element_order(x, ord):
    """
    Computes the order of element `x` of a group
    which has order `ord`.
    """
    factors = factor(ord)
    p = ord
    for (pi,ei) in factors:
        p = p // (pi**ei)
        y = x**p
        while not y:
            y = y**pi
            p *= pi
    
    return p
