#!/usr/bin/python3
# -*- coding: utf-8 -*-

from math import ceil, sqrt


def BSGS(a, b, group, order=None):
    """Baby-step Giant-step algorithm for computing discrete logarithms
    in a finite group `group`.

        INPUT:
        - `a` must be a primitive element of `group`
        - `b` must be an element of `group`

        OUTPUT:
        - `x` such that ``a^x = b``
    """
    if order is None:
        order = group.order()

    dic = {}
    m = ceil(sqrt(order))

    # Baby-step
    aj = group.id()
    for j in range(m):
        dic[aj] = j
        aj = aj * a

    # Giant-step
    c = a ** (-m)
    gam = b
    for i in range(m):
        try:
            j = dic[gam]
            return m * i + j
        except:
            gam *= c