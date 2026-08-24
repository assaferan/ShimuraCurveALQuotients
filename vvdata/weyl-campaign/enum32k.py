#!/usr/bin/env python3
# enum32g generalized to support size K (chunk-safe numpy sweep).
#   enum32k.py <cusp6/mono dump> <out epool file> [R=8] [K=6]
import numpy as np, re, sys, math
from itertools import combinations

dump, outp = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 8
K = int(sys.argv[4]) if len(sys.argv) > 4 else 6

txt = open(dump, errors='replace').read().replace('\\\n','')
mb = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[(.*?)\]', txt)
D, N, M = int(mb.group(1)), int(mb.group(2)), int(mb.group(3))
ds = [int(x) for x in mb.group(4).split(',')]
nd = len(ds)
mvecs = [[int(x) for x in m.group(2).split(',')] for m in re.finditer(r'MONO (\d+) \[(.*?)\]', txt)]
assert mvecs

def vp(d, p):
    v = 0
    while d % p == 0: v += 1; d //= p
    return v
primes = sorted({p for d in ds for p in range(2, d+1) if d % p == 0 and all(p % q for q in range(2, p))})
def s1(r): return sum(d*x for d,x in zip(ds,r))
def s2(r): return sum((M//d)*x for d,x in zip(ds,r))
def parbits(r): return tuple(sum(vp(d,p)*x for d,x in zip(ds,r)) % 2 for p in primes)
chk = {(s1(r) % 24, s2(r) % 24, parbits(r)) for r in mvecs}
assert len(chk) == 1, chk
S1, S2, T = chk.pop()
assert S1 == 0 and S2 == 0
print(f"# base {D}_{N}: M={M}, char bits {dict(zip(primes,T))}, R={R}, K={K}", file=sys.stderr)

cs = [c for c in range(1, M+1) if M % c == 0]
H = np.array([[math.gcd(c,d)**2 * (M//d) for d in ds] for c in cs], dtype=np.int64)
w1 = np.array(ds, dtype=np.int64)
w2 = np.array([M//d for d in ds], dtype=np.int64)
pb = np.array([[vp(d,p) for p in primes] for d in ds], dtype=np.int64)
Tv = np.array(T)

rng = np.arange(-R, R+1)
grids = np.stack(np.meshgrid(*([rng]*(K-1)), indexing='ij'), -1).reshape(-1, K-1)
rlast = 3 - grids.sum(1)
ok = np.abs(rlast) <= 12
grids = grids[ok]; rlast = rlast[ok]
found = set()
ncomb = 0
for sup in combinations(range(nd), K):
    ncomb += 1
    full = np.zeros((len(grids), nd), dtype=np.int64)
    full[:, list(sup[:-1])] = grids
    full[:, sup[-1]] = rlast
    keep = full[((full @ w1) % 24 == 0) & ((full @ w2) % 24 == 0)]
    if not len(keep): continue
    keep = keep[((keep @ pb) % 2 == Tv).all(1)]
    if not len(keep): continue
    keep = keep[(keep @ H.T >= 0).all(1)]
    for row in keep:
        found.add(tuple(int(x) for x in row))
print(f"# combos {ncomb}, holomorphic same-character quotients: {len(found)}", file=sys.stderr)
found = sorted(found)
with open(outp, 'w') as f:
    f.write('[\n' + ',\n'.join('[' + ','.join(str(x) for x in r) + ']' for r in found) + '\n]\n')
