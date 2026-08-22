#!/usr/bin/env python3
# Analysis of cusp6 dumps: derive the m=0 multiplier functional WITHOUT aliasing.
#
#   F(x) = sum_g T_g(x) is linear on the 86-dim monomial space.
#   V = {x : no poles at any cusp class except 1 (cusp 0) and M (cusp oo)},
#       optionally intersected with a fixed eta-multiplier character class.
#   On V, F should equal W . P(x) where P = (pp at oo, pp at 0, and the
#   constants) -- the channels of the two-channel conjecture. Solve W on V,
#   report residual, snap weights, translate oo-pole channels to the
#   (d, k) + square dictionary.
#
# usage: python3 mono6.py cusp6_15_2.out
import sys, re, math
from fractions import Fraction as F

path = sys.argv[1] if len(sys.argv) > 1 else 'cusp6_15_2.out'
txt = open(path, errors='replace').read().replace('\\\n', '')

mB = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[([^\]]+)\]', txt)
D, N, M = int(mB.group(1)), int(mB.group(2)), int(mB.group(3))
ds = [int(s) for s in mB.group(4).split(',')]
print(f"X0^{D}({N})  M={M}")

monos = {}
for m in re.finditer(r'MONO (\d+) \[([^\]]+)\]', txt):
    monos[int(m.group(1))] = [int(s) for s in m.group(2).split(',')]
nm = len(monos)
print(f"#monomials {nm}")

formc = {}
for m in re.finditer(r'FORMC (-?\d+) \[(.*?)\]\n', txt, re.S):
    k = int(m.group(1))
    body = re.sub(r'\s+', ' ', m.group(2))
    pairs = re.findall(r'<\s*(\d+)\s*,\s*([-0-9/]+)\s*>', body)
    formc[k] = {int(a): F(b) for a, b in pairs}
    ntup = body.count('<')
    if ntup != len(pairs):
        print(f"WARNING: form {k}: {ntup} tuples in dump, {len(pairs)} parsed")

TG = {}
for m in re.finditer(r'TG (\d+) g=(\d+) (\S+)\s+(\S+)', txt):
    ri, g = int(m.group(1)), int(m.group(2))
    TG.setdefault(ri, {})[g] = complex(float(m.group(3).replace('E', 'e')),
                                       float(m.group(4).replace('E', 'e')))
PP = {}
for m in re.finditer(r'PP (\d+) g=(\d+) W=(\d+) e=(-?\d+) (\S+)\s+(\S+)', txt):
    ri, g, W, e = (int(m.group(i)) for i in range(1, 5))
    PP.setdefault(g, {}).setdefault((e, W), {})[ri] = complex(
        float(m.group(5).replace('E', 'e')), float(m.group(6).replace('E', 'e')))

mult = {ri: sum(TG.get(ri, {}).values()) for ri in monos}

def snap(x, den=10**6):
    fr = F(x).limit_denominator(den)
    return fr if abs(float(fr) - x) < 1e-20 else None

# ---- consistency: panel totals ----
print("panel totals (c_eta*(0) per form, from per-monomial sums):")
for k in sorted(formc):
    tot = sum(float(c) * mult[ri].real for ri, c in formc[k].items())
    print(f"  form {k:3}: {tot:+.9f}")

# ---- eta-multiplier character class per monomial (Ligozat data) ----
# char determined by (sum d r_d mod 24, sum (M... use n = lcm ds) (n/d) r_d mod 24,
# prod d^{r_d} mod squares). Group monomials; the panel's forms live in one class.
L = 1
for d in ds: L = L * d // math.gcd(L, d)
def charclass(r):
    s1 = sum(d * rd for d, rd in zip(ds, r)) % 24
    s2 = sum((L // d) * rd for d, rd in zip(ds, r)) % 24
    pr = 1
    for d, rd in zip(ds, r):
        pr *= pow(d, rd % 2)
    for p in range(2, 100):
        while pr % (p * p) == 0: pr //= p * p
    return (s1, s2, pr)
cls_of = {ri: charclass(r) for ri, r in monos.items()}
panel_ris = set()
for k in formc: panel_ris |= set(formc[k])
panel_cls = {cls_of[ri] for ri in panel_ris}
print(f"panel character classes: {panel_cls}")

# ---- build V: no poles at classes other than 1 and M ----
rows_con = []
for g in PP:
    if g in (1, M): continue
    for (e, W), vec in PP[g].items():
        if e >= 0: continue
        rows_con.append([vec.get(ri, 0).real for ri in sorted(monos)])
        rows_con.append([vec.get(ri, 0).imag for ri in sorted(monos)])
use = sorted(monos)
if '--panelclass' in sys.argv:
    use = [ri for ri in use if cls_of[ri] in panel_cls]
    rows_con = [[row[sorted(monos).index(ri)] for ri in use] for row in rows_con]
print(f"pole constraints (rows): {len(rows_con)} on {len(use)} monomials")

# kernel via Gaussian elimination (floats)
def kernel(Arows, n, tol=1e-9):
    Ab = [row[:] for row in Arows]
    piv = {}
    pr = 0
    for col in range(n):
        if pr >= len(Ab): break
        p = max(range(pr, len(Ab)), key=lambda i: abs(Ab[i][col]))
        if abs(Ab[p][col]) < tol: continue
        Ab[pr], Ab[p] = Ab[p], Ab[pr]
        s = Ab[pr][col]
        Ab[pr] = [x / s for x in Ab[pr]]
        for i in range(len(Ab)):
            if i != pr and abs(Ab[i][col]) > 1e-14:
                f0 = Ab[i][col]
                Ab[i] = [x - f0 * y for x, y in zip(Ab[i], Ab[pr])]
        piv[col] = pr
        pr += 1
    freec = [c for c in range(n) if c not in piv]
    basis = []
    for fc in freec:
        v = [0.0] * n
        v[fc] = 1.0
        for col, r in piv.items():
            v[col] = -Ab[r][fc]
        basis.append(v)
    return basis

Vb = kernel(rows_con, len(use)) if rows_con else [[1.0 if j == i else 0.0 for j in range(len(use))] for i in range(len(use))]
print(f"dim V = {len(Vb)}")

# ---- channels on V: pp at oo (g=M), pp at 0 (g=1), incl constants e=0 ----
chan = []
chanvec = []
for g in (M, 1):
    tag = 'oo' if g == M else '0'
    for (e, W), vec in sorted(PP.get(g, {}).items(), key=lambda t: t[0][0]):
        if e > 0: continue
        chan.append((tag, e, W))
        chanvec.append([vec.get(ri, 0).real for ri in use])
        im = [vec.get(ri, 0).imag for ri in use]
        if max(abs(x) for x in im) > 1e-30:
            chan.append((tag + 'i', e, W))
            chanvec.append(im)

# system on V: for each basis vector v: sum_c W_c (chanvec_c . v) = F(v)
A = []
b = []
Fv = [mult[ri].real for ri in use]
for v in Vb:
    A.append([sum(cv[i] * v[i] for i in range(len(use))) for cv in chanvec])
    b.append(sum(Fv[i] * v[i] for i in range(len(use))))

# least squares via Gaussian elimination on [A|b] with pivoting
nr, nc = len(A), len(chan)
Ab = [row[:] + [b[i]] for i, row in enumerate(A)]
piv = {}
pr = 0
for col in range(nc):
    if pr >= nr: break
    p = max(range(pr, nr), key=lambda i: abs(Ab[i][col]))
    if abs(Ab[p][col]) < 1e-9: continue
    Ab[pr], Ab[p] = Ab[p], Ab[pr]
    s = Ab[pr][col]
    Ab[pr] = [x / s for x in Ab[pr]]
    for i in range(nr):
        if i != pr and abs(Ab[i][col]) > 1e-14:
            f0 = Ab[i][col]
            Ab[i] = [x - f0 * y for x, y in zip(Ab[i], Ab[pr])]
    piv[col] = pr
    pr += 1
res = [abs(Ab[i][nc]) for i in range(pr, nr)]
print(f"solve on V: rank {pr} / {nc} channels; leftover rows {nr - pr}, max resid {max(res) if res else 0:.3e}")
freec = [c for c in range(nc) if c not in piv]
print("weights (F = c_eta*(0) = 2*mult; divide by 2 for ledger convention):")
for col in sorted(piv):
    i = piv[col]
    dep = [chan[c] for c in freec if abs(Ab[i][c]) > 1e-8]
    w = Ab[i][nc]
    if abs(w) < 1e-9 and dep: continue
    wq = snap(w, 10**4)
    tag, e, W = chan[col]
    lbl = f"{tag} e={e}/{24*W}"
    if tag == 'oo' and e % 24 == 0 and e < 0:
        m0 = -e // 24
        lbl += f" [c_0(-{m0})"
        if m0 % N == 0: lbl += f" N|m d,k={m0 // N}"
        r0 = math.isqrt(m0)
        if r0 * r0 == m0: lbl += " SQUARE"
        lbl += "]"
    print(f"  w[{lbl}] = {w:+.8f} ~ {wq}" + (f"  (dep {dep})" if dep else ""))
