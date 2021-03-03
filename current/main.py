#!/usr/bin/python3

from .primality.baillie_PSW import is_bailliePSW_prime

def is_prime(n:int, alg = 'PSW'):
    if alg == 'PSW':
        return is_bailliePSW_prime(n)
    else:
        raise NotImplementedError