#!/usr/bin/env python3
# Analysis of the cusp5 per-monomial class sums: does the functional
# F(m_r) = sum_g T_g[r] factor through the principal part of m_r, and if so
# what are the channel weights -- solved on the FULL monomial span (rank >> 9),
# de-aliasing what the 9-form panels could not pin.
#
# Monomial m_r = prod_d eta(d tau)^{r_d}, sum r_d = 1 (weight 1/2).
# At oo: m_r = q^{S/24} prod_d prod_n (1 - q^{dn})^{r_d},  S = sum d r_d.
#   pp_oo(m_r) = coefficients at q-exponents <= 0 (grid 1/24).
# At 0 (slash S):  exponents on grid 1/(24 lcm(ds)); same product with d -> 1/d.
#
# usage: python3 mono5.py cusp5_15_2.out
import sys, re, math
from fractions import Fraction as F

path = sys.argv[1] if len(sys.argv) > 1 else 'cusp5_15_2.out'
txt = open(path, errors='replace').read().replace('\\\n', '')

mB = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[([^\]]+)\]', txt)
D, N, M = int(mB.group(1)), int(mB.group(2)), int(mB.group(3))
ds = [int(s) for s in mB.group(4).split(',')]
print(f"X0^{D}({N})  M={M}  ds={ds}")

monos = {}
for m in re.finditer(r'MONO (\d+) \[([^\]]+)\]', txt):
    monos[int(m.group(1))] = [int(s) for s in m.group(2).split(',')]
print(f"#monomials parsed: {len(monos)}")

TG = {}
for m in re.finditer(r'TG (\d+) g=(\d+) (\S+)\s+(\S+)', txt):
    ri, g = int(m.group(1)), int(m.group(2))
    TG.setdefault(ri, {})[g] = complex(float(m.group(3).replace('E','e')),
                                       float(m.group(4).replace('E','e')))

# ---- exact eta-product expansions ----
def eta_prod_pp(r, scale):
    """principal part + constant of prod_d prod_n (1 - q^{d n})^{r_d} * q^{S/24}
    with d -> d*scale... scale=1: at oo (grid 1/24). Returns dict {exp_in_grid: coeff}
    with grid = 1/(24*den), exponents integers on that grid, only exp <= 0.
    For the 0-cusp use r reversed trick: eta(tau/d) ~ d -> 1/d: use den = lcm."""
    # exponents in units of 1/(24*den): d_eff = d*24? We'll compute with den:
    # here scale is the common denominator: exponent of q for factor d is d*24/den...
    # We implement: exponents measured in 1/(24*den); factor (1-q^{dn}) has step dn*24*?
    raise NotImplementedError

def series_pp(steps, S24, depth_extra=0):
    """prod over (step, r) in steps of (1-x^step)^r  * x^{S24}, x = q^{1/grid-unit};
    S24 = leading exponent in grid units (can be negative). Return {e: c} for e <= 0,
    plus full series dict up to 0. Exact integers."""
    depth = max(0, -S24) + 1  # need coefficients of x^k for k in [0, -S24]
    ser = {0: 1}
    for step, r in steps:
        if r == 0 or step >= depth:
            continue
        # multiply by (1 - x^step)^r ; for negative r use expansion sum C(-r+j-1, j) x^{step j}
        new = {}
        if r > 0:
            base = {}
            for j in range(0, r + 1):
                e = step * j
                if e >= depth: break
                base[e] = ((-1) ** j) * math.comb(r, j)
        else:
            base = {}
            j = 0
            while step * j < depth:
                base[step * j] = math.comb(-r + j - 1, j)
                j += 1
        for e1, c1 in ser.items():
            for e2, c2 in base.items():
                e = e1 + e2
                if e < depth:
                    new[e] = new.get(e, 0) + c1 * c2
        ser = new
    return {S24 + e: c for e, c in ser.items() if S24 + e <= 0 and c != 0}

def pp_oo(r):
    S24 = sum(d * ri for d, ri in zip(ds, r))          # units of 1/24
    steps = [(24 * d, ri) for d, ri in zip(ds, r)]     # (1 - q^{dn}) = x^{24 d n}
    return series_pp(steps, S24)                        # exponents in 1/24 units

L = 1
for d in ds: L = L * d // math.gcd(L, d)
def pp_0(r):
    S24 = sum((L // d) * ri for d, ri in zip(ds, r))   # units of 1/(24 L)
    steps = [(24 * (L // d), ri) for d, ri in zip(ds, r)]
    return series_pp(steps, S24)                        # exponents in 1/(24 L) units

# ---- assemble the linear system: F(m_r) = sum_g T_g[r] ----
mult = {}
for ri in sorted(TG):
    tot = sum(TG[ri].values())
    mult[ri] = tot
# monomials with no TG lines at all have F = 0
for ri in monos:
    if ri not in mult: mult[ri] = 0j

# channel index: ('oo', e) for e in 1/24 units <= 0; ('0', e) in 1/(24L) units <= 0
rows, rhs, chans = [], [], []
for ri in sorted(monos):
    r = monos[ri]
    prof = {}
    for e, c in pp_oo(r).items(): prof[('oo', e)] = c
    for e, c in pp_0(r).items(): prof[('0', e)] = c
    for kk in prof:
        if kk not in chans: chans.append(kk)
    rows.append(prof); rhs.append(mult[ri])

chans.sort(key=str)
A = [[float(p.get(c, 0)) for c in chans] for p in rows]
print(f"system: {len(rows)} equations, {len(chans)} channels")

bmax = max(abs(z.imag) for z in rhs)
print("max |imag rhs|:", bmax)
use_complex = bmax > 1e-8

# Gaussian elimination on [A | b] (complex rhs allowed), partial pivoting.
nr, nc = len(A), len(chans)
Ab = [[complex(x) for x in row] + [rhs[i]] for i, row in enumerate(A)]
piv_of_col = {}
prow = 0
for col in range(nc):
    p = max(range(prow, nr), key=lambda i: abs(Ab[i][col]))
    if abs(Ab[p][col]) < 1e-9: continue
    Ab[prow], Ab[p] = Ab[p], Ab[prow]
    pl = Ab[prow]; s = pl[col]
    Ab[prow] = [x / s for x in pl]
    for i in range(nr):
        if i != prow and abs(Ab[i][col]) > 1e-14:
            f0 = Ab[i][col]
            Ab[i] = [x - f0 * y for x, y in zip(Ab[i], Ab[prow])]
    piv_of_col[col] = prow
    prow += 1
rank = prow
print(f"rank(P) = {rank} of {nc} channels, {nr} monomials")
incons = [i for i in range(rank, nr) if abs(Ab[i][nc]) > 1e-6]
print(f"factorization through pp: {'FAILS on ' + str(len(incons)) + ' rows, max resid ' + format(max(abs(Ab[i][nc]) for i in incons), '.3e') if incons else 'HOLDS (all quotient rows vanish)'}")
free = [c for c in range(nc) if c not in piv_of_col]
if free: print(f"free (aliased) channels: {[str(chans[c]) for c in free]}")
for col in sorted(piv_of_col):
    i = piv_of_col[col]
    # pinned iff the row expresses e_col alone: no free-column entries
    dep = [chans[c] for c in free if abs(Ab[i][c]) > 1e-8]
    w = Ab[i][nc]
    tag = 'PINNED' if not dep else f'dep on {[str(d) for d in dep]}'
    if abs(w) > 1e-9 or not dep:
        wq = F(round(w.real * 1680), 1680)
        wq = wq if abs(float(wq) - w.real) < 1e-7 else None
        print(f"  w[{chans[col]}] = {w.real:+.8f}{w.imag:+.2e}i ~ {wq} [{tag}]")
