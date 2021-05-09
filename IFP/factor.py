#!/usr/bin/python3
# -*- coding: utf-8 -*-

from IFP.ecm import ECM
from IFP.pollard_rho import pollard_rho
from sage.all import factor as factor_sage


def factor(n: int, alg="sage"):
    """Returns the factorization of `n` using different algorithms :
    - `'sage'` -> sage's factorization algorithm (optimized)
    - `'ECM'` -> Lenstra's factorization algorithm on elliptic curves
    - `'pollard-rho'` -> Pollard's-rho algorithm with Brent's improvments
    """
    if alg == "sage":
        return factor_sage(n)
    elif alg == "ECM":
        return ECM(n)
    elif alg == "pollard-rho":
        return pollard_rho(n)
    else:
        raise NotImplementedError