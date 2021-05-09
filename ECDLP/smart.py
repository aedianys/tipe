#!/usr/bin/sage
# -*- coding: utf-8 -*-

from sage.all import ZZ, Qp, EllipticCurve


def smart_attack(tildeE, tildeG, tildeA, p, a, b, lift="x"):
    assert tildeE.order() == p  # The curve is anomalous

    E = EllipticCurve(Qp(p), [a, b])  # Lift of curve to Qp

    if lift == "x":
        G = E.lift_x(ZZ(tildeG.xy()[0]))  # Might have to take the other lift.
        A = E.lift_x(ZZ(tildeA.xy()[0]))  # Might have to take the other lift.
    else:
        G = E.lift_y(ZZ(tildeG.xy()[1]))
        A = E.lift_y(ZZ(tildeA.xy()[1]))

    p_times_G = p * G
    p_times_A = p * A
    x_G, y_G = p_times_G.xy()
    x_A, y_A = p_times_A.xy()

    phi_G = -(x_G / y_G)  # Map to Z/pZ
    phi_A = -(x_A / y_A)  # Map to Z/pZ
    k = phi_A / phi_G  # Solve the discrete log in Z/pZ (aka do a division)

    return k.expansion()[0]