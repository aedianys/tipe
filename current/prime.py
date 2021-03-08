#!/usr/bin/python3

from .primality.BPSW import isprime_BPSW
from .primality.naive import isprime_naive
from .primality.miller_rabin import isprime_MillerRabin


def is_prime(n: int, alg="BPSW"):
    """Returns if `n` is prime using different algorithms :
    - `'BPSW'` -> Baillie-PSW probabilistic test
    - `'MR'` -> Miller-Rabin deterministic test for integers `n <= 2⁶⁴`
    - `'naive'` -> Naive deterministic test
    """
    if alg == "BPSW":
        return isprime_BPSW(n)
    elif alg == "MR":
        return isprime_MillerRabin(n)
    elif alg == "naive":
        return isprime_naive(n)
    else:
        raise NotImplementedError