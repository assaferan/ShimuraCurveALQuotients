# Stage 2: the pivot-supported particular solution of constants(beta) = rho.
# If it snaps to rationals, it is a canonical rational representative E_eis
# and W(m) = -a_{E_eis}(m) exactly.
import re, sys, math
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 45

raw = open(sys.argv[1], errors='replace').read().replace('\\\n','')
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
rho = {}
for m in re.finditer(r'RHOV (\d+) (\d+) (\S+) (\S+)', log):
    rho[int(m.group(1))] = mp.mpc(mp.mpf(m.group(3)), mp.mpf(m.group(4)))
pools = {}
for m in re.finditer(r'EPOOL (\d+) \[(.*?)\]', log):
    pools[int(m.group(1))] = [int(x) for x in m.group(2).split(',')]
nE = max(E)
cosets = sorted(set(rho) | {wi for a in E for wi in E[a]})
print(f"# {nE} pool, {len(cosets)} cosets in system")

# real-ified matrix rows over all cosets, columns = pool members
def colvec(a):
    return [E.get(a, {}).get(wi, mp.mpc(0)) for wi in cosets]
cols = [colvec(a) for a in range(1, nE+1)]
rhovec = [rho.get(wi, mp.mpc(0)) for wi in cosets]

# recompute pivots as in stage 1 (row echelon over the real-ification)
rows = []
for i, wi in enumerate(cosets):
    rows.append([cols[a][i].real for a in range(nE)])
    rows.append([cols[a][i].imag for a in range(nE)])
rhs = []
for i, wi in enumerate(cosets):
    rhs.append(rhovec[i].real); rhs.append(rhovec[i].imag)
M = [r[:] + [b] for r, b in zip(rows, rhs)]
piv = []; r = 0
tol = mp.mpf(10)**(-30)
for c in range(nE):
    pr = None; best = tol
    for i in range(r, len(M)):
        if abs(M[i][c]) > best: best = abs(M[i][c]); pr = i
    if pr is None: continue
    M[r], M[pr] = M[pr], M[r]
    pv = M[r][c]
    M[r] = [x/pv for x in M[r]]
    for i in range(len(M)):
        if i != r and abs(M[i][c]) != 0:
            f = M[i][c]
            M[i] = [x - f*y for x, y in zip(M[i], M[r])]
    piv.append(c); r += 1
# consistency: rows beyond rank must have ~0 rhs
maxbad = max((abs(M[i][nE]) for i in range(len(piv), len(M))), default=mp.mpf(0))
print(f"# pivots (0-based pool cols): {piv}; residual beyond rank: {mp.nstr(maxbad, 3)}")
beta = [mp.mpf(0)]*nE
for ri, c in enumerate(piv):
    beta[c] = M[ri][nE]
# snap
snapped = [Fraction(float(x)).limit_denominator(10**6) for x in beta]
err = max(abs(x - mp.mpf(s.numerator)/s.denominator) for x, s in zip(beta, snapped))
print(f"# pivot beta snap err: {mp.nstr(err, 3)}")
for c in piv:
    print(f"  beta[{c+1}] = {snapped[c]}   r = {pools.get(c+1)}")
if float(err) < 1e-25:
    print("RATIONAL E_eis FOUND")
