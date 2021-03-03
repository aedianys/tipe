#!/usr/bin/python3

from functools import reduce

def euclide(r1, r2, u1=1, v1=0, u2=0, v2=1):
    """Extended euclid algorithm"""
    if r2 != 0:
        q = r1//r2
        return euclide(r2, r1-q*r2, u2, v2, u1-q*u2, v1-q*v2)
    else:
        return r1,u1,v1
    
def inverse_mod(x:int, p:int):
    """Modular inverse of x mod p"""
    return (euclide(p,x) % p)[2]

def CRT(residues:list, moduli:list):
    """Return a solution to a Chinese Remainder Theorem.
    
    INPUT:
    residues is a list of residues
    moduli is a list of moduli

    OUTPUT:
    the unique solution to the problem
    """

    N = reduce(lambda x,y: x*y, moduli, 1)
    x = 0

    for r,p in zip(residues,moduli):
        m = N // p
        x += r * euclide(m,p)[1] * m

    return x % N

def factor_naive(n:int):
    """Naive factorization of n"""
    factorization = []

    for d in range(2,n//2):
        q,r = divmod(n,d)
        k = 0
        while r == 0:
            k += 1
            n = q
            q,r = divmod(q,d)
        if k != 0:
            factorization.append((d,k))

    if factorization == []: factorization = [(n,1)]
    
    return factorization

def valuation2(x:int):
    """Computes y,k such that x = y * 2**k and y is odd"""
    k = 0
    while not x & 1:
        k += 1
        x >>= 1
    return x,k

def valuation(n:int,b:int):
    """Computes d,r such that n = d * b**r and b do not divide d"""
    if b == 2: return valuation2(n)
    r,d = 0,n
    q,s = divmod(n,b)
    while s == 0:
        r += 1
        d = q
        q,s = divmod(d,b)
    return d,r

def jacobi_symbol(a:int, n:int):
    """Jacobi symbol (a/n) where n must be odd strictly positive integer""" 
    assert(n > 0 and n & 1)
    a = a % n
    t = 1
    while a != 0:
        q,r = divmod(a,2)
        while r == 0:
            a = q
            s = n & 7
            if s == 3 or s == 5: t = -t
            q,r = divmod(a,2)

        a,n = n,a
        if a & 3 == 3 and n & 3 == 3: t = -t
        a = a % n

    if n == 1: return t
    else: return 0

def intsqrt(n):
    if n < 0:
        raise ValueError
    elif n < 2:
        return n
    else:
        a = 1 << ((1 + n.bit_length()) >> 1)
        while True:
            b = (a + n // a) >> 1
            if b >= a:
                return a
            a = b

def is_perfect_square(x:int):
    """Test if x is a perfect square"""
    mask = 0xC840C04048404040
    if ((mask << x) >> 63) & 1: return False
    d,r = valuation2(x)
    if r & 1: return False
    if d & 7 != 1 or d <= 0: return d == 0
    return intsqrt(d)**2 == d
 
def isprime_naive(n):
    if n < 5:
        return n == 2 or n == 3
    elif n%2 == 0:
        return False
    else:
        r = intsqrt(n)
        k = 3
        while k <= r:
            if n%k == 0:
                return False
            k += 2
        return True

if __name__ == "__main__":
    pass
