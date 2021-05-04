#!/usr/bin/python3

from functools import reduce


def euclide(a: int, b: int):
    """
    Extended euclid algorithm
    INPUT:
    - `a`, `b` two integers

    OUTPUT:
    - `d`,`u`,`v` such that `a*u + b*v = d = gcd(a,b)`
    """

    r1, r2 = a, b
    u1, u2 = 1, 0
    v1, v2 = 0, 1

    while r2:
        q, r = divmod(r1,r2)
        r1, r2 = r2, r
        u1, u2 = u2, u1 - q * u2
        v1, v2 = v2, v1 - q * v2

    return (r1,u1,v1)


def gcd(a: int, b: int):
    return euclide(a, b)[0]


def inverse_mod(x: int, p: int):
    """Modular inverse of x mod p"""
    return (euclide(x, p) % p)[1]


def CRT(residues: list, moduli: list):
    """
    Returns a solution to a Chinese Remainder Theorem.

    INPUT:
    - residues is a list of residues
    - moduli is a list of moduli

    OUTPUT:
    - the unique solution to the problem
    """

    N = reduce(lambda x, y: x * y, moduli, 1)
    x = 0

    for r, p in zip(residues, moduli):
        m = N // p
        x += r * euclide(m, p)[1] * m

    return x % N


def factor_naive(n: int):
    """Naive factorization of n"""
    factors = []

    for d in range(2, n // 2):
        q, r = divmod(n, d)
        k = 0
        while r == 0:
            k += 1
            n = q
            q, r = divmod(q, d)
        if k != 0:
            factors.append((d, k))

    if factors == []:
        factors = [(n, 1)]

    return factors


def valuation2(x: int):
    """Computes `y`,`k` such that `x = y * 2**k` and `y` is odd"""
    k = 0
    while not x & 1:
        k += 1
        x >>= 1
    return x, k


def valuation(n: int, b: int):
    """Computes `d`,`r` such that `n = d * b**r` and `b` do not divide `d`"""
    if b == 2:
        return valuation2(n)
    r, d = 0, n
    q, s = divmod(n, b)
    while s == 0:
        r += 1
        d = q
        q, s = divmod(d, b)
    return d, r


def jacobi_symbol(a: int, n: int):
    """Jacobi symbol `(a/n)` where `n` must be odd strictly positive integer"""
    assert n > 0 and n & 1
    a = a % n
    t = 1
    while a != 0:
        q, r = divmod(a, 2)
        while r == 0:
            a = q
            s = n & 7
            if s == 3 or s == 5:
                t = -t
            q, r = divmod(a, 2)

        a, n = n, a
        if a & 3 == 3 and n & 3 == 3:
            t = -t
        a = a % n

    if n == 1:
        return t
    else:
        return 0


def int_sqrt(n: int):
    """
    Computes the floor of the square root of `n`,
    i.e. the biggest integer `a` such that `a² <= n`.
    """
    if n < 0:
        raise ValueError("n must be non negative")
    elif n in (0, 1):
        return n
    else:
        a = 1 << ((1 + n.bit_length()) >> 1)
        while True:
            b = (a + n // a) >> 1
            if b >= a:
                return a
            a = b


def int_nthroot(y: int, n: int):
    """
    Return x such that `x = floor(y**(1/n))`.
    """
    if y < 0:
        raise ValueError("y must be nonnegative")
    if n < 1:
        raise ValueError("n must be positive")
    if y in (0, 1):
        return y
    if n == 1:
        return y
    if n == 2:
        return int_sqrt(y)
    if n > y:
        return 1
       
    x = int(y ** (1.0 / n))
    
    return x


def is_perfect_square(x: int):
    """Test if `x` is a perfect square."""
    mask = 0xC840C04048404040
    if ((mask << x) >> 63) & 1:
        return False
    d, r = valuation2(x)
    if r & 1:
        return False
    if d & 7 != 1 or d <= 0:
        return d == 0
    return int_sqrt(d) ** 2 == d


def factorial(n: int):
    assert n >= 0
    f = 1
    while n != 0:
        f *= n
        n -= 1
    return f


def primes(n):
    """List all primes lesser than `n` using Eratosthenes' sieve"""
    tab = [True] * (n + 1)
    ps = []
    for p in range(2, n + 1):
        if tab[p]:
            ps.append(p)
            for i in range(p, n + 1, p):
                tab[i] = False
    return ps


def array_op(l1: list, l2: list, op: function):
    assert len(l1) == len(l2)
    return [op(a,b) for a,b in zip(l1,l2)]


def array_scal(x, l: list, op: function):
    return [op(x,a) for a in l]


def array_is_zero(l):
    return all([e==0 for e in l])

def array_is_inde(l1:list, l2:list, inv, add, mult):
    """Tests if l1 and l2 are independant vectors."""
    assert len(l1) == len(l2)
    if array_is_zero(l1):
        return not array_is_zero(l2)
    else:
        i=0
        while l1[i] == 0: i += 1
        if l2[i] == 0: return True
        else:
            return array_scal(inv(l1[i]), l1, mult) != array_scal(inv(l2[i]), l2, mult)


def row_echelon_form(m: list, inv, add, mult):
    """Returns the row echelon form of matrix `m`
    using Gaussian elimination"""
    p = len(m)
    q = len(m[0])
    for col in range(q):
        i = col
        while i < p and m[i][col] == 0:
            i += 1
        if i != p:
            m[i], m[col] = m[col], m[i]
            k = inv(m[col][col])
            m[col] = array_scal(k, m[col])
            for j in range(i+1,p):
                l = array_scal(-m[j][col], m[col], mult)
                m[j] = array_op(m[j], l, add)
    return m


def reduced_row_echelon_form(m: list, inv, add, mult):
    """Returns the row echelon form of matrix `m`
    using Gaussian elimination"""
    p = len(m)
    q = len(m[0])
    for col in range(q):
        i = col
        while i < p and m[i][col] == 0:
            i += 1
        if i != p:
            m[i], m[col] = m[col], m[i]
            k = inv(m[col][col])
            m[col] = array_scal(k, m[col])
            for j in range(p):
                if j != i:
                    l = array_scal(-m[j][col], m[col], mult)
                    m[j] = array_op(m[j], l, add)
    return m


if __name__ == "__main__":
    pass