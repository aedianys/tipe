#!/usr/bin/python3
# -*- coding: utf-8 -*-

from primality.BPSW import is_prime_BPSW
from primality.naive import is_prime_naive
from primality.miller_rabin import is_prime_MillerRabin
from sage.all import is_prime as is_prime_sage


def is_prime(n: int, alg="sage"):
    """Returns if `n` is prime using different algorithms :
    - `'sage'` -> sage's probabilistic test (optimized)
    - `'BPSW'` -> Baillie-PSW probabilistic test
    - `'MR'` -> Miller-Rabin deterministic test for integers `n <= 2⁶⁴`
    - `'naive'` -> Naive deterministic test
    """
    if alg == "sage":
        return is_prime_sage(n)
    if alg == "BPSW":
        return is_prime_BPSW(n)
    elif alg == "MR":
        return is_prime_MillerRabin(n)
    elif alg == "naive":
        return is_prime_naive(n)
    else:
        raise NotImplementedError