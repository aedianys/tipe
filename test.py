#!/usr/bin/python3

from current.main import is_prime as prime
from sage.all import is_prime

print(any(is_prime(x) != prime(x) for x in range(100000)))