#!/usr/bin/python3

from ..maths.math_lib import jacobi_symbol, is_perfect_square
from .miller_rabin import is_MillerRabin_prime
from .lucas import is_lucas_strong_prime

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]

def D_choice(n):
    D = 5
    while True:
       j = jacobi_symbol(D,n)
       if j == -1:
           break
       D = -D+2 if D < 0 else -D-2
       if abs(D) > 11:
           if is_perfect_square(n):
               return False
    return D

def is_bailliePSW_prime(n:int):
    """Pseudo-primality test using the baillie-PSW primality test"""
    if n in SMALL_PRIMES:
        return True

    if any((n % p == 0) for p in SMALL_PRIMES) or n in (0,1):
        return False

    if not is_MillerRabin_prime(n,2):
        return False
    
    D = D_choice(n)
    if not D:
        return False
    
    P,Q = 1, (1-D)//4
    if not is_lucas_strong_prime(n,P,Q):
        return False
    
    return True 

if __name__ == "__main__":
    print(is_bailliePSW_prime(
            4547337172376300111955330758342147474062293202868155909489))
    print(is_bailliePSW_prime(1511670205808207663))