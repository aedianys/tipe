#!/usr/bin/python3
# -*- coding: utf-8 -*-


def vector_add(v1: list, v2: list):
    """Vector addition"""
    assert len(v1) == len(v2)
    return [x + y for x, y in zip(v1, v2)]


def vector_scal(scalar, v: list):
    """Scalar multiplication of a vector"""
    return [scalar * x for x in v]


def vector_is_zero(v):
    """Tests if v is the zero vector."""
    return all([not e for e in v])


def vectors_are_inde(v1: list, v2: list):
    """Tests if v1 and v2 are independant vectors."""
    assert len(v1) == len(v2)

    if vector_is_zero(v1):
        return not vector_is_zero(v2)
    else:
        i = 0
        while not v1[i]:
            i += 1

        if not v2[i]:
            return True
        else:
            v1_reduced = vector_scal(v1[i].inverse() * v2[i], v1)
            return v1_reduced != v2


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