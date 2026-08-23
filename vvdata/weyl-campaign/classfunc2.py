#!/usr/bin/env python3
# Per-class functional, exact version: snap S_g to rationals (den <= 24), build
# the exact linear system S_g(f) = sum_ch W_g(ch) c_f(ch) over oo + 0-side
# channels, rref over Q, report rank / forced weights / inconsistency, and the
# inter-class linear relations among the S_g row vectors.
import re, collections, math, glob, os, sys
from fractions import Fraction

CAMP = sys.argv[1] if len(sys.argv) > 1 else '.'

def ff(s): return float(s.replace('\\','').replace('E','e'))

def snap(x, maxden=24, tol=None):
    f = Fraction(x).limit_denominator(maxden)
    return f

def load_cusp3(path):
    txt = open(path, errors='replace').read()
    m = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)', txt)
    D,N,M = int(m.group(1)), int(m.group(2)), int(m.group(3))
    rows = collections.defaultdict(list)
    joined = txt.replace('\\\n','')
    for chunk in joined.split('COSET ')[1:]:
        chunk = chunk.split('\nFORM')[0].replace('\n','')
        m2 = re.match(r'(-?\d+) c=(\d+) d=(\d+) gcd=(\d+) a0=\((\S+?),(\S+?)\) rho=\((\S+?),(\S+?)\) contrib=\((\S+?),(\S+?)\)', chunk)
        if not m2: continue
        k=int(m2.group(1)); g=int(m2.group(4))
        ct=complex(ff(m2.group(9)),ff(m2.group(10)))
        rows[(k,g)].append(ct)
    S = {}
    for (k,g),cts in rows.items():
        tot = sum(cts)
        S[(k,g)] = tot
    return D,N,M,S

def load_pp(D,N):
    for stem in (f'pp3_{D}_{N}.out', f'pp2_{D}_{N}.out'):
        p = os.path.join(CAMP, stem)
        if os.path.exists(p): break
    else:
        return None
    pp = collections.defaultdict(dict)
    txt = open(p, errors='replace').read().replace('\\\n','')
    for line in txt.splitlines():
        t = line.split()
        if not t or t[0] != 'COEF': continue
        _, d_, n_, k, cusp, m1, m2, c = t
        if int(d_) != D or int(n_) != N: continue
        pp[int(k)][(cusp, m1)] = Fraction(c)
    return pp

def rref(Ab, ncols):
    # rows of Fractions, length ncols+1 (augmented); returns (pivots, rows)
    rows = [r[:] for r in Ab]
    piv = []
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r,len(rows)) if rows[i][c] != 0), None)
        if pr is None: continue
        rows[r], rows[pr] = rows[pr], rows[r]
        rows[r] = [x/rows[r][c] for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f*y for x,y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
        if r == len(rows): break
    return piv, rows

for path in sorted(glob.glob(os.path.join(CAMP,'cusp3_*_full.out'))):
    D,N,M,S = load_cusp3(path)
    pp = load_pp(D,N)
    if pp is None: continue
    ks = sorted({k for k,_ in S}); gs = sorted({g for _,g in S})
    scale = max(abs(v) for v in S.values())
    tol = max(1e-6*scale, 1e-6)
    print(f"\n== X0^{D}({N}) M={M} scale {scale:.1e}")
    chans = sorted({ch for k in ks for ch in pp.get(k,{})},
                   key=lambda ch: (ch[0], Fraction(ch[1]) if '/' in ch[1] else int(ch[1])))
    # live classes with rational-snapped S rows
    Srows = {}
    for g in gs:
        ys = [S.get((k,g), 0) for k in ks]
        if max(abs(y) for y in ys) < tol: continue
        if max(abs(y.imag) for y in ys) > tol:
            print(f"  g={g}: IMAGINARY class sum, max imag {max(abs(y.imag) for y in ys):.2e}")
        sn = [snap(y.real) for y in ys]
        err = max(abs(float(a)-y.real) for a,y in zip(sn,ys))
        if err > max(2e-4*scale, 1e-6):
            print(f"  g={g}: snap err {err:.2e} at scale {scale:.1e} -- keeping (10-digit prints)")
        Srows[g] = sn
    lives = sorted(Srows)
    # inter-class relations: rref of the matrix whose rows are S_g vectors
    mat = [[*Srows[g]] for g in lives]
    piv, rows = rref([r+[Fraction(0)] for r in mat], len(ks))
    rank = len(piv)
    print(f"  live classes {lives}; rank of S-rows {rank}/{len(lives)}")
    if rank < len(lives):
        # nullspace of the transpose: coefficients lam with sum lam_g S_g = 0
        # solve via rref on transpose
        tmat = [[Srows[g][i] for g in lives] for i in range(len(ks))]
        pivT, rowsT = rref([r+[Fraction(0)] for r in tmat], len(lives))
        free = [j for j in range(len(lives)) if j not in pivT]
        for fj in free:
            lam = [Fraction(0)]*len(lives); lam[fj] = Fraction(1)
            for ri,c in enumerate(pivT):
                lam[c] = -rowsT[ri][fj]
            den = 1
            for l in lam: den = den*l.denominator//math.gcd(den, l.denominator)
            lamZ = [l*den for l in lam]
            terms = " ".join(f"{'+' if l>0 else '-'}{abs(l)}*S_{g}" for l,g in zip(lamZ,lives) if l)
            print(f"    RELATION: {terms} = 0")
    # per-class exact solve over channels
    for g in lives:
        Ab = []
        for i,k in enumerate(ks):
            row = [pp[k].get(ch, Fraction(0)) for ch in chans] + [Srows[g][i]]
            Ab.append(row)
        piv, rows = rref(Ab, len(chans))
        # inconsistency?
        bad = [r for r in rows if all(x==0 for x in r[:-1]) and r[-1] != 0]
        nf = len(chans) - len(piv)
        forced = []
        for ri, c in enumerate(piv):
            # pivot col c forced iff row ri has zeros in all free cols
            if all(rows[ri][j] == 0 for j in range(len(chans)) if j != c and j not in piv):
                forced.append((chans[c], rows[ri][-1]))
        msg = f"  g={g:>3}: rank {len(piv)}, {nf} free"
        if bad: msg += "  INCONSISTENT (oo+0 channels insufficient)"
        if forced:
            msg += "  forced: " + " ".join(f"W({c[0]},{c[1]})={v}" for c,v in forced)
        print(msg)
print("\ndone")
