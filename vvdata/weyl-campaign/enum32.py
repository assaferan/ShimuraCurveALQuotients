#!/usr/bin/env python3
# Enumerate ALL holomorphic weight-3/2 eta quotients on Gamma_0(60) with the
# 15_2 panel character, supported on <= 5 of the 12 divisors, entries in
# [-8, 8].  Conditions on the exponent vector r over ds:
#   sum r_d = 3;  s1 = sum d r_d = 0 mod 24;  s2 = sum (60/d) r_d = 0 mod 24;
#   squarefree class of prod d^{r_d} = the panel class t;
#   holomorphy: sum_d r_d gcd(c,d)^2 (60/d) >= 0 for every c | 60.
# Output: a Magma sequence literal for eis32.m's EF hook.
import numpy as np, re, sys, math
from itertools import combinations

CAMP = '/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/e91ada1b-fe8a-4a7b-a418-93559b9c5c11/scratchpad/campaign/vvdata/weyl-campaign'
ds = [1,2,3,4,5,6,10,12,15,20,30,60]
cs = ds
nd = len(ds)

# panel character from the banked cusp6 MONO lines
mono_txt = open(f'{CAMP}/cusp6_15_2.out', errors='replace').read().replace('\\\n','')
mvecs = []
for m in re.finditer(r'MONO (\d+) \[(.*?)\]', mono_txt):
    mvecs.append([int(x) for x in m.group(2).split(',')])
assert mvecs, "no MONO lines"
def s1(r): return sum(d*x for d,x in zip(ds,r))
def s2(r): return sum((60//d)*x for d,x in zip(ds,r))
def parbits(r):
    b = []
    for p in (2,3,5):
        e = 0
        for d,x in zip(ds,r):
            v = 0; dd = d
            while dd % p == 0: v += 1; dd //= p
            e += v*x
        b.append(e % 2)
    return tuple(b)
chk = {(s1(r) % 24, s2(r) % 24, parbits(r)) for r in mvecs}
assert len(chk) == 1, f"panel monomials do not share one character: {chk}"
(S1, S2, T) = chk.pop()
print(f"# panel character: s1={S1} s2={S2} mod 24, square class bits (2,3,5)={T}", file=sys.stderr)
assert S1 == 0 and S2 == 0

# holomorphy matrix (integers): H[c][d] = gcd(c,d)^2 * (60//d)
H = np.array([[math.gcd(c,d)**2 * (60//d) for d in ds] for c in cs], dtype=np.int64)
d_arr = np.array(ds, dtype=np.int64)
w1 = d_arr.copy()                      # s1 weights
w2 = np.array([60//d for d in ds], dtype=np.int64)
pb = np.array([[ (lambda dd,p: sum(1 for _ in iter(lambda: None, None)) )(d,p) for p in (2,3,5)] for d in ds])
# recompute p-adic valuations properly
def vp(d,p):
    v=0
    while d % p == 0: v+=1; d//=p
    return v
pb = np.array([[vp(d,p) for p in (2,3,5)] for d in ds], dtype=np.int64)

R = 8
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
    m1 = (full @ w1) % 24 == 0
    m2 = (full @ w2) % 24 == 0
    keep = full[m1 & m2]
    if not len(keep): continue
    pbv = (keep @ pb) % 2
    keep = keep[(pbv == np.array(T)).all(1)]
    if not len(keep): continue
    hol = (keep @ H.T >= 0).all(1)
    keep = keep[hol]
    for row in keep:
        found.add(tuple(int(x) for x in row))

found = sorted(found)
print(f"# holomorphic same-character weight-3/2 quotients found: {len(found)}", file=sys.stderr)
with open(sys.argv[1] if len(sys.argv) > 1 else 'epool_15_2.txt', 'w') as f:
    f.write('[\n' + ',\n'.join('[' + ','.join(str(x) for x in r) + ']' for r in found) + '\n]\n')
