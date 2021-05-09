from primality.prime import is_prime
from IFP.factor import factor
from random import randint
from sage.all import EllipticCurve, GF
from point_counting.schoof import schoof
from groups.elliptic_curve import EllipticCurve as EC
from fields.finite import FiniteField
from time import time

def measure_time(func, *args):
    t = time()
    res = func(*args)
    print(f"Function took {time()-t} seconds to run")
    return res

#print(time("challenge", "challenges.exceptional_curves", ""))
a,b,p = 46, 74, 97
a,b,p = (1333, 1129, 3571)
""" E = EllipticCurve(GF(p), [a, b])
print(E.order()) """
F = FiniteField(p)
E2 = EC(F, a, b)
print(measure_time(schoof, E2))

""" E = EllipticCurve(GF(p), [a, b])
print(E.order())
F = FiniteField(p)
E2 = EC(F, a, b)
print(measure_time(schoof, E2))
"""
""" ec_test_values = [
(13, 215, 229),
(106, 166, 197),
(31, 16, 137),
(503, 367, 523),
(1333, 1129, 3571)
]

for a,b,p in ec_test_values:
    #E = EllipticCurve(GF(p), [a, b])
    #print(E.order())
    F = FiniteField(p)
    E2 = EC(F, a, b)
    print(measure_time(schoof, E2))
 """