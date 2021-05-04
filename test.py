from current.prime import is_prime
from current.factor import factor
from sage.all import factor as s_factor, is_prime as s_is_prime
import timeit
from random import randint


def time(function, source=None, *args):
    if source != None:
        t = timeit.Timer(
            "{}({})".format(function, *args),
            setup="from {} import {}".format(source, function),
        )
    else:
        t = timeit.Timer("{}({})".format(function, *args))
    return t.autorange()


# factorization correctness checks
N = 9876432123412789339 * 14484968830081347923412479832498112841
print(factor(7*5*13*53*31*47))
print(factor(1359561509*8169704801))
print(factor(N))

"""
# primality time checks
print(time('is_prime', 'current.main', '9876432123412789339'))
print(time('is_prime', 'sage.all', '9876432123412789339'))
print(time('is_prime', 'current.main', '14484968830081347923412479832498112841'))
print(time('is_prime', 'sage.all', '14484968830081347923412479832498112841'))

N = 100
MAX = 0x100000000000000000000000000000000
SUM = 0

for i in range(N):
    n = randint(2,MAX)
    prime = is_prime(n)
    assert _is_prime(n) == prime
    x,tx = time('is_prime', 'current.main', str(n))
    y,ty = time('is_prime', 'sage.all', str(n))
    rel_diff = (tx/x - ty/y)/(ty/y)
    SUM += rel_diff
    print(f'{i:02d} - {prime} - écart : {rel_diff}')

print(SUM/N)

# primality correctness checks
print(any(is_prime(x) != _is_prime(x) for x in range(100000)))
"""