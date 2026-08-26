#!/usr/bin/env python3
# pooltransfer.py <M> <out epool> [epool dir]
#
# Build an eta pool at level M for free, out of the pools already banked at
# every level that DIVIDES M.
#
# The handoff records that "the eta pool transfers by level M, not by base"
# (15_1 reused epool_15_2, etc.) -- that is the special case M_src = M.  The
# general statement is cheaper and much more useful: if M_src | M then every
# divisor of M_src is a divisor of M, so a level-M_src exponent vector embeds
# in the level-M one by padding with zeros, and the quotient is still modular
# on Gamma_0(M) since Gamma_0(M) < Gamma_0(M_src).
#
# Embedding is necessary but NOT sufficient, so each transferred vector is
# re-checked at the TARGET level against exactly what enum32m/eis32 demand:
#   * weight:     sum r = 3
#   * character:  s1 = sum d*r     = 0 mod 24
#                 s2 = sum (M/d)*r = 0 mod 24      <- this one changes with M
#                 parity bits sum vp(d,p)*r mod 2 = 1 at p = 2, 0 elsewhere
#   * holomorphy: ord_c = sum r_d gcd(c,d)^2/d >= 0 at EVERY cusp c | M
#                                                   <- M has new cusps
# s2 and the cusp list both depend on M, so this is a real filter, not a
# formality: at M = 4620 it rejects a third of the candidates.
#
# At M = 4620 this yields 294 distinct vectors in about a second, against the
# ~9 h that enum32m R=8 K=7 needs to enumerate the level from scratch.  It is
# not a replacement for enum32m in general -- the transferred vectors have NO
# support on divisors of M that divide no source level (231, 385, 1155, ...
# at 4620), so the pool lives in a proper coordinate subspace.  The residual
# is the arbiter, as always.
import ast, glob, os, sys
from fractions import Fraction
from math import gcd

M = int(sys.argv[1])
outp = sys.argv[2]
pdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(os.path.abspath(sys.argv[0]))

ds = [d for d in range(1, M + 1) if M % d == 0]
primes = sorted({p for p in range(2, M + 1) if M % p == 0
                 and all(p % q for q in range(2, p))})

def vp(d, p):
    v = 0
    while d % p == 0: v += 1; d //= p
    return v

def accept(r):
    if sum(r) != 3: return False
    if sum(d * x for d, x in zip(ds, r)) % 24: return False
    if sum((M // d) * x for d, x in zip(ds, r)) % 24: return False
    par = tuple(sum(vp(d, p) * x for d, x in zip(ds, r)) % 2 for p in primes)
    if par != tuple(1 if p == 2 else 0 for p in primes): return False
    for c in ds:
        if sum(Fraction(x * gcd(c, d) ** 2, d) for d, x in zip(ds, r)) < 0:
            return False
    return True

found = {}
for f in sorted(glob.glob(os.path.join(pdir, "epool_*.txt"))):
    base = os.path.basename(f)[6:-4]
    try:
        D, N = [int(x) for x in base.split('_')]
    except ValueError:
        continue
    Msrc = 4 * D * N if (D * N) % 2 else 2 * D * N
    if M % Msrc:
        print(f"  skip {base:8s} M={Msrc:5d}  does not divide {M}")
        continue
    dsrc = [d for d in range(1, Msrc + 1) if Msrc % d == 0]
    vecs = ast.literal_eval(open(f).read())
    new = 0
    for r in vecs:
        if len(r) != len(dsrc):
            print(f"  WARNING {base}: length {len(r)} != {len(dsrc)} divisors, skipped")
            break
        big = [0] * len(ds)
        for d, x in zip(dsrc, r):
            big[ds.index(d)] = x
        t = tuple(big)
        if accept(big):
            if t not in found: new += 1
            found[t] = 1
    print(f"  {base:8s} M={Msrc:5d}  {len(vecs):4d} vectors -> {new:4d} new accepted")

rows = sorted(found)
with open(outp, 'w') as fh:
    fh.write('[\n' + ',\n'.join('[' + ','.join(map(str, r)) + ']' for r in rows) + '\n]\n')
print(f"\nTRANSFER POOL at M = {M}: {len(rows)} distinct vectors -> {outp}")
