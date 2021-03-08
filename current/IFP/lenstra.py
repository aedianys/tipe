#!/usr/bin/python3

from random import randint
from math import gcd
from .primepower import isprimepower
from ..prime import is_prime


# Sieve of Eratosthenes
def primes(n):
    b = [True] * (n + 1)
    ps = []
    for p in range(2, n + 1):
        if b[p]:
            ps.append(p)
            for i in range(p, n + 1, p):
                b[i] = False
    return ps


# Finds modular inverse
# Returns inverse, unused helper and gcd
def modular_inv(a, b):
    if b == 0:
        return 1, 0, a
    q, r = divmod(a, b)
    x, y, g = modular_inv(b, r)
    return y, x - q * y, g


# Addition in Elliptic curve modulo m space
def elliptic_add(p, q, a, b, m):
    # If one point is infinity, return other one
    if p[2] == 0:
        return q
    if q[2] == 0:
        return p
        
    if p[0] == q[0]:
        if (p[1] + q[1]) % m == 0:
            return 0, 1, 0  # Infinity
        num = (3 * p[0] * p[0] + a)
        denom = (2 * p[1])
    else:
        num = (q[1] - p[1])
        denom = (q[0] - p[0])
    inv, _, g = modular_inv(denom, m)
    # Unable to find inverse, arithmetic breaks
    if g > 1:
        return 0, 0, denom  # Failure
    z = (num * inv * num * inv - p[0] - q[0]) % m
    return z, (num * inv * (p[0] - z) - p[1]) % m, 1


# Multiplication (repeated addition and doubling)
def elliptic_mul(k, p, a, b, m):
    r = (0, 1, 0)  # Infinity
    while k > 0:
        # p is failure, return it
        if p[2] > 1:
            return p
        if k % 2 == 1:
            r = elliptic_add(p, r, a, b, m)
        k = k // 2
        p = elliptic_add(p, p, a, b, m)
    return r


# Lenstra's algorithm for factoring
# Limit specifies the amount of work permitted
def lenstra(n, limit=100000, times=5):
    for _ in range(times):
        g = n
        while g == n:
            # Randomized x and y
            q = randint(0, n - 1), randint(0, n - 1), 1
            # Randomized curve coefficient a, computed b
            a = randint(0, n - 1)
            b = (q[1] * q[1] - q[0] * q[0] * q[0] - a * q[0]) % n
            g = gcd(4 * a * a * a + 27 * b * b, n)  # singularity check
        # If we got lucky, return lucky factor
        if g > 1:
            return g
        # increase k step by step until lcm(1, ..., limit)
        for p in primes(limit):
            pp = p
            while pp < limit:
                q = elliptic_mul(p, q, a, b, n)
                # Elliptic arithmetic breaks
                if q[2] > 1:
                    return gcd(q[2], n)
                pp = p * pp
    return n

def factor(n, limit=100000):
    factors = []

    for p in (2, 3):
        q, r = divmod(n, p)
        while r == 0:
            n = q
            factors.append(p)
            q, r = divmod(n, p)

    m,r = isprimepower(n)
    if r != 1:
        return [m]*r

    if is_prime(n):
        return [n]

    while n != 1:
        d = lenstra(n, limit)
        factors.append(d)
        n //= d

    return sorted(factors)


if __name__ == "__main__":
    print(lenstra(7*5*13*53*31*47))
    print(lenstra(1359561509*8169704801))