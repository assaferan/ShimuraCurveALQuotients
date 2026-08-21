#!/usr/bin/env python3
# The constraints ledger: for every base with ground truth, print the exact linear
# constraints the measured multipliers place on the theta-class weights w(d, k)
# (index m = N*d*k^2, d squarefree). Solve where determined. This is the data any
# future closed form / Serre-Stark theta computation must reproduce.

import glob, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(__file__))
from kappa_conjecture import parse, gtsweep_truth, square_free_split  # reuse
from weylfit import rref

here = os.path.dirname(__file__)
gtt = gtsweep_truth(here)
seen = set()
for p in sorted(glob.glob(os.path.join(here, 'pp*_*.out'))):
    try:
        (D, N, M, det), forms, coefs = parse(p)
    except Exception:
        continue
    if (D, N) in seen: continue
    seen.add((D, N))
    if (D, N) in gtt:
        for row, val in gtt[(D, N)].items():
            if forms.get(row) is None: forms[row] = val
    known = [k for k in sorted(forms) if forms[k] is not None]
    if not known: continue
    # collect per-form equations over w[(d,k)]
    idxs = []
    rows, rhs = [], []
    for kf in known:
        r = {}
        for (blk, m), c in coefs[kf].items():
            if blk != 'oo' or m % N: continue
            d, k = square_free_split(m // N)
            r[(d, k)] = r.get((d, k), F(0)) + c
        for key in r:
            if key not in idxs: idxs.append(key)
        rows.append(r); rhs.append(forms[kf])
    idxs.sort()
    A = [[r.get(i, F(0)) for i in idxs] for r in rows]
    R, r2, piv = rref(A, rhs)
    incons = [i for i in range(len(R)) if all(x == 0 for x in R[i]) and r2[i] != 0]
    print(f"\nX0^{D}({N})  [{len(rows)} eqs, {len(idxs)} classes"
          f"{', INCONSISTENT!' if incons else ''}]  truths: "
          + ", ".join(f"{k}:{forms[k]}" for k in known))
    if not idxs:
        print("   (no supported indices; all multipliers must be 0)"
              + ("  OK" if all(v == 0 for v in rhs) else "  *** VIOLATED ***"))
        continue
    free = [c for c in range(len(idxs)) if c not in piv]
    for ri, ci in enumerate(piv):
        terms = [f"w{idxs[ci]}"]
        for cj in free:
            if R[ri][cj] != 0:
                terms.append(f"{'+' if -R[ri][cj] > 0 else '-'} {abs(R[ri][cj])}*w{idxs[cj]}")
        val = r2[ri]
        if len(terms) == 1:
            print(f"   w{idxs[ci]} = {val}   <-- FORCED")
        else:
            print(f"   w{idxs[ci]} = {val} {' '.join(terms[1:])}")
    for cj in free:
        print(f"   w{idxs[cj]} free")
