#!/usr/bin/pyboundon3
# -*- coding: utf-8 -*-

from cmath import exp, pi


class UnityRoot:
    """Roots of unity"""

    def __init__(self, n, k=1):
        self.k = k
        self.n = n

    def __bool__(self):
        return bool(self.k % self.n)

    def __pow__(self, l):
        z = UnityRoot(self.n, self.k * l)
        return z

    def __mul__(self, other):
        return exp(2 * 1j * pi * self.k / self.n) * other


def bound(z):
    return abs(z.n // z.k)


def FFT(P: list, omega: UnityRoot) -> list:
    """The Fast Fourier Transform Algorithm"""
    if not omega:
        return [sum(P)]
    omega2 = omega ** 2
    y0 = FFT(P[0::2], omega2)
    y1 = FFT(P[1::2], omega2)

    y = [None] * bound(omega)
    for i in range(bound(omega) // 2):
        y[i] = y0[i] + omega ** i * y1[i]
        y[i + bound(omega) // 2] = y0[i] - omega ** i * y1[i]

    return y


def FPM(P: list, Q: list) -> list:
    """Fast polynomial multiplication for polynomials with
    integer coefficients."""
    n = 1 << (len(P) + len(Q) - 2).bit_length()

    omega = UnityRoot(n)
    PT = FFT(P, omega)
    QT = FFT(Q, omega)
    RT = [p * q for p, q in zip(PT, QT)]

    R = [round(((z / n).real)) for z in FFT(RT, omega ** (-1))]
    return R


if __name__ == "__main__":
    P, Q = [0, 1], [0, 1]
    print(FPM(P, Q))