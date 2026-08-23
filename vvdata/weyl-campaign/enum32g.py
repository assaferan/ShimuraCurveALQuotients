#!/usr/bin/env python3
# General-base version of enum32.py: enumerate ALL holomorphic weight-3/2 eta
# quotients on Gamma_0(M) with the panel character of a base, reading the base
# data (M, ds, panel monomials) from its banked cusp6_{D}_{N}.out dump.
# Conditions on r over ds: sum r = 3; s1 = sum d r_d = 0 mod 24;
# s2 = sum (M/d) r_d = 0 mod 24; squarefree class of prod d^{r_d} = panel's;
# holomorphy sum_d r_d gcd(c,d)^2 (M/d) >= 0 for all c | M.
# Support <= 5 divisors, exponents in [-R, R].
#   enum32g.py <cusp6 dump> <out epool file> [R=8]
import numpy as np, re, sys, math
from itertools import combinations

dump, outp = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 8

txt = open(dump, errors='replace').read().replace('\\\n','')
mb = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[(.*?)\]', txt)
D, N, M = int(mb.group(1)), int(mb.group(2)), int(mb.group(3))
ds = [int(x) for x in mb.group(4).split(',')]
nd = len(ds)
mvecs = []
for m in re.finditer(r'MONO (\d+) \[(.*?)\]', txt):
    mvecs.append([int(x) for x in m.group(2).split(',')])
assert mvecs, "no MONO lines"
print(f"# base {D}_{N}: M={M}, {nd} divisors, {len(mvecs)} panel monomials", file=sys.stderr)

def vp(d, p):
    v = 0
    while d % p == 0: v += 1; d //= p
    return v
primes = sorted({p for d in ds for p in range(2, d+1) if d % p == 0 and all(p % q for q in range(2, p))})
def s1(r): return sum(d*x for d,x in zip(ds,r))
def s2(r): return sum((M//d)*x for d,x in zip(ds,r))
def parbits(r):
    return tuple(sum(vp(d,p)*x for d,x in zip(ds,r)) % 2 for p in primes)
chk = {(s1(r) % 24, s2(r) % 24, parbits(r)) for r in mvecs}
assert len(chk) == 1, f"panel does not share one character: {chk}"
S1, S2, T = chk.pop()
assert S1 == 0 and S2 == 0, (S1, S2)
print(f"# panel character: square class bits {dict(zip(primes,T))}", file=sys.stderr)

cs = [c for c in range(1, M+1) if M % c == 0]
H = np.array([[math.gcd(c,d)**2 * (M//d) for d in ds] for c in cs], dtype=np.int64)
w1 = np.array(ds, dtype=np.int64)
w2 = np.array([M//d for d in ds], dtype=np.int64)
pb = np.array([[vp(d,p) for p in primes] for d in ds], dtype=np.int64)
Tv = np.array(T)

rng = np.arange(-R, R+1)
found = set()
for sup in combinations(range(nd), 5):
    i1,i2,i3,i4,i5 = sup
    g = np.stack(np.meshgrid(rng, rng, rng, rng, indexing='ij'), -1).reshape(-1,4)
    r5 = 3 - g.sum(1)
    ok = np.abs(r5) <= 12
    g = g[ok]; r5 = r5[ok]
    full = np.zeros((len(g), nd), dtype=np.int64)
    full[:, [i1,i2,i3,i4]] = g
    full[:, i5] = r5
    keep = full[((full @ w1) % 24 == 0) & ((full @ w2) % 24 == 0)]
    if not len(keep): continue
    keep = keep[((keep @ pb) % 2 == Tv).all(1)]
    if not len(keep): continue
    keep = keep[(keep @ H.T >= 0).all(1)]
    for row in keep:
        found.add(tuple(int(x) for x in row))

found = sorted(found)
print(f"# holomorphic same-character weight-3/2 quotients: {len(found)}", file=sys.stderr)
with open(outp, 'w') as f:
    f.write('[\n' + ',\n'.join('[' + ','.join(str(x) for x in r) + ']' for r in found) + '\n]\n')
