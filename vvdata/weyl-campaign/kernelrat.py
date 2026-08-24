# Stage 1 of the E* rationalization: the cusp subspace of the pool.
# Parse EMAT (constants of the 21 pool quotients at all cosets, 50 digits),
# compute the kernel of the constants map numerically at high precision,
# echelonize, snap to rationals, and VERIFY the snapped kernel exactly
# against the 50-digit data.
import re, sys
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 45

raw = open(sys.argv[1], errors='replace').read().replace('\\\n','')
# magma -b wraps long lines WITHOUT backslashes: join continuation lines
KEY = ('EMAT', 'RHOV', 'BASE', 'SOLVE', 'RHO ', '#', 'EPOOL', 'BETA', 'ECOEF', 'E ', 'oo-', 'ESTAR')
lines, cur = [], ''
for line in raw.splitlines():
    if any(line.startswith(k) for k in KEY):
        if cur: lines.append(cur)
        cur = line
    else:
        cur += line
lines.append(cur)
log = '\n'.join(lines)
E = {}
for m in re.finditer(r'EMAT (\d+) (\d+) (\d+) (\S+) (\S+)', log):
    a, wi = int(m.group(1)), int(m.group(2))
    E.setdefault(a, {})[wi] = mp.mpc(mp.mpf(m.group(4)), mp.mpf(m.group(5)))
nE = max(E)
cosets = sorted({wi for a in E for wi in E[a]})
print(f"# {nE} pool members, {len(cosets)} cosets with nonzero constants")
# matrix: rows = cosets (Re and Im separately), cols = pool members
rows = []
for wi in cosets:
    rows.append([E.get(a, {}).get(wi, mp.mpc(0)).real for a in range(1, nE+1)])
    rows.append([E.get(a, {}).get(wi, mp.mpc(0)).imag for a in range(1, nE+1)])
A = mp.matrix(rows)
# kernel via column-reduction: reduce A to column echelon; kernel = combos
# Use Gaussian elimination on A^T A? Simpler: SVD unavailable in mpmath for rect;
# do QR-less: row-reduce A (rows >> cols), track rank and back-solve nullspace.
M = [[A[i,j] for j in range(nE)] for i in range(A.rows)]
piv = []
r = 0
tol = mp.mpf(10)**(-30)
for c in range(nE):
    # find pivot
    pr = None; best = tol
    for i in range(r, len(M)):
        if abs(M[i][c]) > best: best = abs(M[i][c]); pr = i
    if pr is None: continue
    M[r], M[pr] = M[pr], M[r]
    pv = M[r][c]
    M[r] = [x/pv for x in M[r]]
    for i in range(len(M)):
        if i != r and abs(M[i][c]) > 0:
            f = M[i][c]
            M[i] = [x - f*y for x, y in zip(M[i], M[r])]
    piv.append(c); r += 1
free = [c for c in range(nE) if c not in piv]
print(f"# rank {len(piv)}, kernel dim {len(free)}")
kernel = []
for fc in free:
    v = [mp.mpf(0)]*nE; v[fc] = mp.mpf(1)
    for ri, c in enumerate(piv):
        v[c] = -M[ri][fc]
    kernel.append(v)
# snap to rationals and verify exactly against the data
for ki, v in enumerate(kernel):
    snapped = [Fraction(float(x)).limit_denominator(10**6) for x in v]
    err = max(abs(x - mp.mpf(s.numerator)/s.denominator) for x, s in zip(v, snapped))
    # verify: sum snapped_a * E[a][wi] ~ 0 for all cosets
    dev = mp.mpf(0)
    for wi in cosets:
        tot = mp.mpc(0)
        for a in range(1, nE+1):
            s = snapped[a-1]
            if s: tot += mp.mpf(s.numerator)/s.denominator * E.get(a, {}).get(wi, mp.mpc(0))
        dev = max(dev, abs(tot))
    den = 1
    for s in snapped: den = den*s.denominator // __import__('math').gcd(den, s.denominator)
    ints = [s*den for s in snapped]
    print(f"K{ki}: snap_err {mp.nstr(err,3)}  exact_dev {mp.nstr(dev,3)}  den {den}")
    print("   ", [int(x) for x in ints])
