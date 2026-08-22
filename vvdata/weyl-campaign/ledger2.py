#!/usr/bin/env python3
# The TWO-CHANNEL ledger (post square-class discovery, 2026-08-22):
#     mult_f = sum over N|m classes  w(d, k) c_0(-m)   [level channel, m = N d k^2]
#            + w_sq * sum_k c_0(-k^2)                  [square-class / Zagier channel]
#
# Truth sources, in priority order:
#   1. cusp1_<D>_<N>.out  -- the cusp-expansion oracle's full panels (half: values)
#   2. gtsweep gt_*.out   -- TRUE multiplier lines / residual pairs (odd N)
#   3. FORM lines in pp*_ dumps
# Principal parts always from the pp*_ dumps (validated against the cusp scan on 22_3).
#
# Output per base: consistency, spares, w_sq (forced / free / absent), level-channel
# weights, and the cross-base tables for pattern hunting.

import glob, os, re, math
from fractions import Fraction as F
from kappa_conjecture import parse, gtsweep_truth, square_free_split
from weylfit import rref

HERE = os.path.dirname(__file__) or '.'

def cusp_truths(here):
    """Parse cusp1_<D>_<N>.out: FORM k: ... (half: v) lines -> {(D,N): {k: v}}.
    Uses the first nontrivial isotropic coset's half-value; requires all isotropic
    cosets to agree (they always have) and |Im| tiny. Values snapped to rationals
    with denominator <= 48."""
    out = {}
    for p in sorted(glob.glob(os.path.join(here, 'cusp1_*_*.out'))):
        mb = re.match(r'.*cusp1_(\d+)_(\d+)(?:_full)?\.out$', p)
        if not mb: continue
        D, N = int(mb.group(1)), int(mb.group(2))
        txt = open(p, errors='replace').read()
        # stitch wrapped lines
        txt = re.sub(r'\n(?!FORM|  c_eta|BASE|#)', '', txt)
        forms = {}
        for fm in re.split(r'FORM ', txt)[1:]:
            k = int(fm.split(':')[0])
            vals = []
            for cm in re.finditer(r'c_eta\(0\)\[(\d+)\][^=]*= (\S+) \+ (\S+) i\s+\(half:\s*(\S+)\)', fm):
                j = int(cm.group(1)); re_, im_, half = (float(cm.group(2)), float(cm.group(3)),
                                                        float(cm.group(4)))
                vals.append((j, re_, im_, half))
            triv = [v for v in vals if v[0] == 1]
            nontriv = [v for v in vals if v[0] != 1]
            if not nontriv: continue
            # gates: trivial ~ 0, all nontrivial agree, imaginary parts tiny
            if triv and abs(triv[0][1]) > 1e-25: continue
            hs = [v[3] for v in nontriv]
            if max(hs) - min(hs) > 1e-20: continue
            if max(abs(v[2]) for v in vals) > 1e-25: continue
            h = hs[0]
            fr = F(round(h * 48), 48)
            if abs(float(fr) - h) > 1e-18: continue
            forms[k] = fr
        if forms:
            key = (D, N)
            if key in out: out[key].update(forms)
            else: out[key] = forms
    return out

def is_square(m):
    r = math.isqrt(m); return r * r == m

def solve_base(D, N, forms, coefs, truths):
    known = [k for k in sorted(truths)]
    rows, rhs, idx = [], [], []
    for k in known:
        if k not in coefs: continue
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
        rows.append(prof); rhs.append(truths[k])
    idx.sort(key=str)
    A = [[r.get(i, F(0)) for i in idx] for r in rows]
    R, r2, piv = rref(A, rhs)
    incons = [i for i in range(len(R)) if all(x == 0 for x in R[i]) and r2[i] != 0]
    free = [c for c in range(len(idx)) if c not in piv]
    sol = {}
    dep = {}
    for ri, ci in enumerate(piv):
        sol[idx[ci]] = r2[ri]
        dd = {idx[cj]: -R[ri][cj] for cj in free if R[ri][cj] != 0}
        if dd: dep[idx[ci]] = dd
    return dict(eqs=len(rows), idx=idx, piv=piv, free=[idx[c] for c in free],
                sol=sol, dep=dep, incons=bool(incons),
                spare=len(rows) - len(piv))

def main():
    ct = cusp_truths(HERE)
    gt = gtsweep_truth(HERE)
    seen = set()
    sq_table, level_table = {}, {}
    for p in sorted(glob.glob(os.path.join(HERE, 'pp*_*.out'))):
        try:
            (D, N, M, det), forms, coefs = parse(p)
        except Exception:
            continue
        if (D, N) in seen: continue
        seen.add((D, N))
        truths = {}
        for k, v in forms.items():
            if v is not None: truths[k] = v
        for src in (gt.get((D, N), {}), ct.get((D, N), {})):
            for k, v in src.items():
                if k in truths and truths[k] != v:
                    print(f"  !! X0^{D}({N}) form {k}: truth conflict {truths[k]} vs {v}")
                truths[k] = v
        if not truths: continue
        r = solve_base(D, N, forms, coefs, truths)
        sq = ('absent' if ('sq',) not in r['idx'] else
              ('FREE' if ('sq',) in r['free'] else
               (r['sol'][('sq',)], r['dep'].get(('sq',)))))
        tag = 'INCONSISTENT' if r['incons'] else 'ok'
        nsrc = len(ct.get((D, N), {}))
        print(f"X0^{D}({N}): eqs {r['eqs']} (cusp-oracle {nsrc}), rank {len(r['piv'])}, "
              f"spare {r['spare']} -> {tag}; w_sq = {sq}")
        if r['incons']:
            print("    *** channels insufficient -- a third channel exists here ***")
        pretty = {}
        for i in r['idx']:
            if i == ('sq',): continue
            if i in r['sol']:
                v = r['sol'][i]
                d = r['dep'].get(i)
                pretty[str(i)] = f"{v}" + (f" (+dep on {list(map(str, d))})" if d else "")
            else:
                pretty[str(i)] = 'FREE'
        print("   ", pretty)
        sq_table[(D, N)] = sq
        level_table[(D, N)] = r
    print("\n== w_sq across bases ==")
    for k in sorted(sq_table): print(f"  X0^{k[0]}({k[1]}): {sq_table[k]}")

if __name__ == '__main__':
    main()
