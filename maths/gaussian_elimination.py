#!/usr/bin/python3
# -*- coding: utf-8 -*-


def vector_add(l1: list, l2: list):
    assert len(l1) == len(l2)
    return [a + b for a, b in zip(l1, l2)]


def vector_scal(x, l: list):
    return [x * a for a in l]


def vector_is_zero(l):
    return all([not e for e in l])


def vector_is_inde(l1: list, l2: list):
    """Tests if l1 and l2 are independant vectors."""
    assert len(l1) == len(l2)
    if vector_is_zero(l1):
        return not vector_is_zero(l2)
    else:
        i = 0
        while not l1[i]:
            i += 1

        if not l2[i]:
            return True
        else:
            l1_reduced = vector_scal(l1[i].inverse() * l2[i], l1)
            return l1_reduced != l2


def reduced_row_echelon_form(system: list):
    """Returns the reduced row echelon form of matrix `system`
    using Gaussian elimination"""
    p = len(system)
    q = len(system[0])
    for col in range(q):
        i = col
        while i < p and not system[i][col]:
            i += 1

        if i != p:
            system[i], system[col] = system[col], system[i]
            k = system[col][col].inverse()
            system[col] = vector_scal(k, system[col])

            for j in range(p):
                if j != i:
                    l = vector_scal(-system[j][col], system[col])
                    system[j] = vector_add(system[j], l)

    return system