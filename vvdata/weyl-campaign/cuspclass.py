#!/usr/bin/env python3
# Per-cusp-class analysis of the cusp3 dumps: for each base and each cusp class
# g = gcd(c, M), the class sum S_g(f) over the panel's forms is itself a linear
# functional of f's principal part. Extract the per-class values, snap to
# rationals, and solve each class's mini-ledger over the same (d,k)+square
# channels as ledger2 -- the closed form decomposes class by class.

import re, sys, glob, os, math, collections
from fractions import Fraction as F
from kappa_conjecture import parse, square_free_split
from weylfit import rref

HERE = os.path.dirname(__file__) or '.'

def ff(s): return float(s.replace('\\', '').replace('E', 'e'))

def read_dump(path):
    """-> {form: {gcdclass: complex sum}}, M"""
    txt = open(path, errors='replace').read()
    Mm = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)', txt)
    D, N, M = int(Mm.group(1)), int(Mm.group(2)), int(Mm.group(3))
    out = collections.defaultdict(lambda: collections.defaultdict(complex))
    for chunk in txt.split('COSET ')[1:]:
        chunk = chunk.split('\nFORM')[0].replace('\\\n', '').replace('\n', '')
        m = re.match(r'(-?\d+) c=(\d+) d=(\d+) gcd=(\d+) a0=\((\S+?),(\S+?)\) '
                     r'rho=\((\S+?),(\S+?)\) contrib=\((\S+?),(\S+?)\)', chunk)
        if not m: continue
        k = int(m.group(1)); g = int(m.group(4))
        out[k][g] += complex(ff(m.group(9)), ff(m.group(10)))
    return D, N, M, out

def snap(x, den=48):
    fr = F(round(x * den), den)
    return fr if abs(float(fr) - x) < 1e-15 else None

def is_square(m):
    r = math.isqrt(m); return r * r == m

def main():
    for path in sorted(glob.glob(os.path.join(HERE, 'cusp3_*_full.out'))):
        D, N, M, data = read_dump(path)
        ppf = os.path.join(HERE, f'pp3_{D}_{N}.out')
        if not os.path.exists(ppf): ppf = os.path.join(HERE, f'pp2_{D}_{N}.out')
        coefs = None
        if os.path.exists(ppf):
            _, _, coefs = parse(ppf)
        ks = sorted(data)
        classes = sorted({g for k in ks for g in data[k]})
        print(f"\n== X0^{D}({N})  M = {M}: per-class sums S_g(f)  [rows: forms] ==")
        hdr = "form  " + "".join(f"g={g:<8}" for g in classes) + "total"
        print(hdr)
        table = {}
        for k in ks:
            vals = []
            for g in classes:
                z = data[k].get(g, 0)
                v = snap(z.real, den=48)
                if v is None or abs(z.imag) > 1e-18:
                    v = snap(z.real, den=2*48*48)   # generous second try
                vals.append(v)
            tot = sum(v for v in vals if v is not None)
            table[k] = dict(zip(classes, vals))
            print(f"{k:4}  " + "".join(f"{str(v):<10}" for v in vals) + str(tot))
        if coefs is None:
            print("   (no pp dump; skipping per-class solves)")
            continue
        # per-class mini-ledger: S_g(f) = sum over channels of weights * coeffs
        for g in classes:
            rows, rhs, idx = [], [], []
            for k in ks:
                if k not in coefs or table[k][g] is None: continue
                prof = {}
                for (blk, m), c in coefs[k].items():
                    if blk != 'oo': continue
                    if m % N == 0:
                        d, kk = square_free_split(m // N)
                        prof[('dk', d, kk)] = prof.get(('dk', d, kk), F(0)) + c
                    if is_square(m):
                        prof[('sq',)] = prof.get(('sq',), F(0)) + c
                for i in prof:
                    if i not in idx: idx.append(i)
                rows.append(prof); rhs.append(table[k][g])
            idx.sort(key=str)
            A = [[r.get(i, F(0)) for i in idx] for r in rows]
            R, r2, piv = rref(A, rhs)
            incons = any(all(x == 0 for x in R[i]) and r2[i] != 0 for i in range(len(R)))
            sol = {str(idx[c]): str(r2[ri]) for ri, c in enumerate(piv)}
            free = [str(idx[c]) for c in range(len(idx)) if c not in piv]
            print(f"   class g={g}: {'INCONSISTENT (needs another channel)' if incons else 'ok'}"
                  f"  weights {sol}" + (f"  free {free}" if free else ""))

if __name__ == '__main__':
    main()
