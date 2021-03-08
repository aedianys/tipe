#!/usr/bin/python3

from .IFP.ecm import ECM

def factor(n: int, alg="ECM"):
    """Returns the factorization of `n` using Lenstra's factorization
    algorithm on elliptic curves"""
    if alg == "ECM":
        return ECM(n)
    else:
        raise NotImplementedError