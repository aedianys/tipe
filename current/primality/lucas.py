#!/usr/bin/python3

from ..maths.math_lib import valuation2

def is_lucas_strong_prime(n:int, P:int, Q:int, jacob = -1):
    """Test if n is a strong Lucas probable prime using parameters D, P and Q"""
    delta = n - jacob
    D = P*P - 4*Q
    d,s = valuation2(delta)
    U,V,Qk = 1,P,Q
    digits = list(map(int, str(bin(d))[2:]))[1:]
    
    if P == 1 and Q == -1: #makes calculations easier
        m = -1
        for digit in digits:
            U = (U*V) % n
            V = (V*V - 2*m) % n
            m = 1

            if digit == 1:
                U, V = P*U + V, D*U + P*V
                U = ((U + (n if U & 1 else 0)) >> 1) % n
                V = ((V + (n if V & 1 else 0)) >> 1) % n
                m *= -1

        
        if U == 0 or V == 0:
            return True
        
        m = -1
        for r in range(1,s):
            V = (V*V - 2*m) % n
            m = 1
            if V == 0:
                return True
    else:
        for digit in digits:
            U = (U*V) % n
            V = (V*V - 2*Qk) % n
            Qk = (Qk*Qk) % n
            if digit == 1:
                U, V = P*U + V, D*U + P*V
                U = ((U + (n if U & 1 else 0)) >> 1) % n
                V = ((V + (n if V & 1 else 0)) >> 1) % n
                Qk = (Qk*Q) % n
        
        if U == 0 or V == 0:
            return True
        
        for r in range(1,s):
            V, Qk = (V*V - 2*Qk) % n, (Qk*Qk) % n
            if V == 0:
                return True
        
    return False
