#!/usr/bin/env python3
# enum32m: meet-in-the-middle replacement for enum32k's per-support grid sweep.
# Same interface + output format:  enum32m.py <mono dump> <out epool> [R=8] [K=7]
# Covers supports of size <= K (K = a + b split, default 3+4) with entries in [-R, R].
# Strategy: for each b-subset B (the K-a largest support indices), match congruence
# keys of its value-grid against cached keys of every a-subset A below min(B).
# A match satisfies ALL congruences (s1, s2 mod 24; parity bits; sum r = 3) by
# construction; holomorphy (H r >= 0) is then checked cusp-by-cusp with early exit.
import numpy as np, re, sys, math
from itertools import combinations

dump, outp = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 8
K = int(sys.argv[4]) if len(sys.argv) > 4 else 7
NA = 3 if K >= 6 else K // 2            # size of the low-index half
NB = K - NA

txt = open(dump, errors='replace').read().replace('\\\n', '')
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
np_ = len(primes)
def s1(r): return sum(d*x for d, x in zip(ds, r))
def s2(r): return sum((M//d)*x for d, x in zip(ds, r))
def parbits(r): return tuple(sum(vp(d, p)*x for d, x in zip(ds, r)) % 2 for p in primes)
chk = {(s1(r) % 24, s2(r) % 24, parbits(r)) for r in mvecs}
assert len(chk) == 1, chk
S1, S2, T = chk.pop()
assert S1 == 0 and S2 == 0
print(f"# base {D}_{N}: M={M}, char bits {dict(zip(primes, T))}, R={R}, K={K}={NA}+{NB}", file=sys.stderr)

cs = [c for c in range(1, M+1) if M % c == 0]
H = np.array([[math.gcd(c, d)**2 * (M//d) for d in ds] for c in cs], dtype=np.int64)
w1 = np.array(ds, dtype=np.int64) % 24
w2 = np.array([M//d for d in ds], dtype=np.int64) % 24
pb = np.array([[vp(d, p) for p in primes] for d in ds], dtype=np.int64) % 2
Tv = np.array(T, dtype=np.int64)

# key packing: (s1 mod 24, s2 mod 24, np_ parity bits, sum+SOFF) -> single int64
SOFF = R * K + 4                       # sums lie in [-R*K, R*K+3]
SRANGE = 2 * SOFF
PARM = 1 << np_
def pack(a1, a2, par, sm):
    return ((a1 * 24 + a2) * PARM + par) * SRANGE + sm

vals = np.arange(-R, R + 1, dtype=np.int64)   # 2R+1 values incl. 0
nv = len(vals)

def side_keys(idxs, n):
    """all value-assignments on n coords idxs: returns (rows [cnt,n], key components)."""
    grids = np.stack(np.meshgrid(*([vals] * n), indexing='ij'), -1).reshape(-1, n)
    a1 = grids @ w1[list(idxs)] % 24
    a2 = grids @ w2[list(idxs)] % 24
    par = np.zeros(len(grids), dtype=np.int64)
    for j, ii in enumerate(idxs):
        par ^= ((grids[:, j] % 2) * 0)  # placeholder, replaced below
    # parity bits: sum vp(d,p)*r mod 2 per prime, packed little-endian over primes
    parbits_ = (grids @ pb[list(idxs)]) % 2
    par = parbits_ @ (1 << np.arange(np_, dtype=np.int64))
    sm = grids.sum(1)
    return grids, a1, a2, par, sm

# cache all A-side (NA-subset) sorted key tables
Acache = {}
def get_A(idxs):
    if idxs not in Acache:
        grids, a1, a2, par, sm = side_keys(idxs, len(idxs))
        key = pack(a1, a2, par, sm + SOFF)
        order = np.argsort(key, kind='stable')
        Acache[idxs] = (grids[order], key[order])
    return Acache[idxs]

Hlist = H  # (ncusp, nd)
found = set()
ncand = 0
nB = 0
targ_par = int(Tv @ (1 << np.arange(np_, dtype=np.int64)))

for Bidx in combinations(range(nd), NB):
    lo = Bidx[0]
    na = min(NA, lo)
    Bgrids, b1, b2, bpar, bsm = side_keys(Bidx, NB)
    nB += 1
    # needed A-key so the pair satisfies all congruences exactly
    need_a1 = (-b1) % 24
    need_a2 = (-b2) % 24
    need_par = bpar ^ targ_par
    need_sm = 3 - bsm
    ok = np.abs(need_sm) <= R * NA               # A side cannot reach beyond this
    need = pack(need_a1[ok], need_a2[ok], need_par[ok], need_sm[ok] + SOFF)
    Bok = Bgrids[ok]
    for Aidx in combinations(range(lo), na) if na > 0 else [()]:
        if na == 0:
            # B alone must hit the target key
            hit = (b1 == 0) & (b2 == 0) & (bpar == targ_par) & (bsm == 3)
            cand_A = np.zeros((int(hit.sum()), 0), dtype=np.int64)
            cand_B = Bgrids[hit]
            pairsA, pairsB = cand_A, cand_B
        else:
            Agrids, Akey = get_A(Aidx)
            L = np.searchsorted(Akey, need, side='left')
            Rr = np.searchsorted(Akey, need, side='right')
            cnt = Rr - L
            hit = cnt > 0
            if not hit.any(): continue
            # expand multi-matches
            reps = cnt[hit]
            brow = np.repeat(np.nonzero(hit)[0], reps)
            starts = L[hit]
            arow = np.concatenate([np.arange(s, s + c) for s, c in zip(starts, reps)]) if len(reps) else np.zeros(0, dtype=np.int64)
            pairsA, pairsB = Agrids[arow], Bok[brow]
        if not len(pairsB): continue
        ncand += len(pairsB)
        # holomorphy, cusp by cusp with early exit
        cols = list(Aidx) + list(Bidx)
        rvals = np.concatenate([pairsA, pairsB], axis=1)
        alive = np.ones(len(rvals), dtype=bool)
        Hs = Hlist[:, cols]
        for c in range(Hs.shape[0]):
            if not alive.any(): break
            alive[alive] = (rvals[alive] @ Hs[c]) >= 0
        for row, colset in ((rvals[i], cols) for i in np.nonzero(alive)[0]):
            full = [0] * nd
            for j, cix in enumerate(colset): full[cix] += int(row[j])
            found.add(tuple(full))

print(f"# B-subsets {nB}, congruence candidates {ncand}, holomorphic same-character quotients: {len(found)}", file=sys.stderr)
found = sorted(found)
with open(outp, 'w') as f:
    f.write('[\n' + ',\n'.join('[' + ','.join(str(x) for x in r) + ']' for r in found) + '\n]\n')
