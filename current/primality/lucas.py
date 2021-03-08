#!/usr/bin/python3

from ..maths.math_lib import valuation2, jacobi_symbol


def ispseudoprime_Lucas(n: int, P: int, Q: int, jacob=-1):
    """Test if `n` is a strong Lucas probable prime using parameters `D`, `P` and `Q`.
    `jacob` is the jacobi symbol `(D/n)`
    """
    delta = n - jacob
    D = P * P - 4 * Q
    d, s = valuation2(delta)
    U, V, Qk = 1, P, Q
    digits = list(map(int, str(bin(d))[2:]))[1:]
    strong_prime = False

    if P == 1 and Q == -1:  # makes calculations easier
        m = -1
        for digit in digits:
            U = (U * V) % n
            V = (V * V - 2 * m) % n
            m = 1

            if digit == 1:
                U, V = P * U + V, D * U + P * V
                U = ((U + (n if U & 1 else 0)) >> 1) % n
                V = ((V + (n if V & 1 else 0)) >> 1) % n
                m *= -1

        if U == 0 or V == 0:
            strong_prime = True

        m = -1
        for r in range(1, s):
            V = (V * V - 2 * m) % n
            m = 1
            if V == 0:
                strong_prime = True

        if jacob == -1:
            V = (V * V - 2 * m) % n
            if (V + 2) % n != 0:
                return False

    else:
        for digit in digits:
            U = (U * V) % n
            V = (V * V - 2 * Qk) % n
            Qk = (Qk * Qk) % n
            if digit == 1:
                U, V = P * U + V, D * U + P * V
                U = ((U + (n if U & 1 else 0)) >> 1) % n
                V = ((V + (n if V & 1 else 0)) >> 1) % n
                Qk = (Qk * Q) % n

        if U == 0 or V == 0:
            strong_prime = True

        for r in range(1, s):
            V, Qk = (V * V - 2 * Qk) % n, (Qk * Qk) % n
            if V == 0:
                strong_prime = True

        if not strong_prime:
            return False

        if jacob == -1:
            V = (V * V - 2 * Qk) % n
            if (V - 2 * Q) % n != 0:
                return False

            if abs(Q) != 1:
                if Qk != (Q * jacobi_symbol(Q, n)) % n:
                    return False

    return True
