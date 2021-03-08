#!/usr/bin/sage

from sage.all import EllipticCurve, GF
from .decrypt import decrypt_flag
from ..ECDLP.pohlig_hellman import ECC_PH

# Define the curve
p = 310717010502520989590157367261876774703
a = 2
b = 3
E = EllipticCurve(GF(p), [a, b])

# generator's public key
G = E.point(
    (179210853392303317793440285562762725654, 105268671499942631758568591033409611165)
)

# alice's public key
A = E.point(
    (280810182131414898730378982766101210916, 291506490768054478159835604632710368904)
)
# bob's public key
B = E.point(
    (272640099140026426377756188075937988094, 51062462309521034358726608268084433317)
)


# Ciphertext
iv = "07e2628b590095a5e332d397b8a59aa7"
encrypted_flag = "8220b7c47b36777a737f5ef9caa2814cf20c1c1ef496ec21\
                  a9b4833da24a008d0870d3ac3a6ad80065c138a2ed6136af"

n = ECC_PH(G, A, E)
secret = (n * B).xy()[0]
print(decrypt_flag(secret, iv, encrypted_flag))
