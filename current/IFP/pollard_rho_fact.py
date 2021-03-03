#!/usr/bin/python3

from math import gcd,log

def pseudo_random(n:int, f):
    if f == None:
        def g(x): return (x**2 + 1) % n
    else:
        def g(x): return f(x) % n
    return g

def fact_PR(n:int, f=None, z=2):
    """Integer factorization using Pollard Rho algorithm."""
    g = pseudo_random(n,f)
    puiss = 2
    y = z

    while puiss < 4*n:
        x = y
        for k in range(puiss):
            y = g(y)
            factor = gcd(x-y,n)
            if factor > 1:
                if factor == n:
                    return fact_PR(n, f, z+1)
                else:
                    return factor
        puiss *= 2
    
    return 1 

def fact_PR_opti(n:int, f=None, y=2, step=100):
    """Integer factorization using Pollard Rho algorithm, with the multiplication improvement."""
    g = pseudo_random(n,f)
    last,cur = 2,0
    prod = 1
    puiss = 2

    while puiss < 4*n:
        x = y
        for k in range(puiss):
            cur += 1
            y = g(y)
            prod *= x-y
            if cur == step:
                factor = gcd(prod,n)
                if factor > 1:
                    if factor == n:
                        return fact_PR(n, f, last)
                    else:
                        return factor
                cur = 0
                last = y 
                    
        puiss *= 2
    
    return 1 


if __name__ == "__main__":
    print(fact_PR(144))
    print(fact_PR_opti(144))
    print(fact_PR_opti(2**32 + 1))
    print(fact_PR_opti(196))


