#!/usr/bin/python3

from .IFP.ecm import ECM
from .IFP.pollard_rho import pollard_rho

def factor(n: int, alg="ECM"):
    """Returns the factorization of `n` using different algorithms :
    - `'ECM'` -> Lenstra's factorization algorithm on elliptic curves
    - `'pollard-rho'` -> Pollard's-rho algorithm with Brent's improvments
    """
    if alg == "ECM":
        return ECM(n)
    elif alg == "pollard-rho":
        return pollard_rho(n)
    else:
        raise NotImplementedError