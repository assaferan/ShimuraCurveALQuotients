#!/usr/bin/env python3
"""
Compare candidate SUPPORT RULES for the m=0 weights u_m, base by base.

A rule says which principal-part indices may carry a nonzero weight.  Imposing it turns the measured
multipliers into an overdetermined linear system; the evidence for a rule is
   (a) CONSISTENCY on every base, and
   (b) how many SPARE CONDITIONS it satisfies (equations beyond the rank -- checks that could fail),
   (c) how few residual free parameters it leaves (a rule that explains rather than absorbs).

Rules:
   diff   -- the space fails to represent m at exactly one finite place   (last session's rule)
   level  -- N | m                                                       (this session's rule)
   both   -- the intersection

usage: rule_test.py [base ...]
"""
import sys, os
from fractions import Fraction as F
import support_fit as S

RULES = ["diff", "level", "both"]


def load(base):
    for cand in (os.path.join(S.WT, f"allterms_{base}.txt"),
                 os.path.join(S.WT, "atdump", f"allterms_{base}.txt")):
        if os.path.exists(cand) and "DONE" in open(cand).read():
            return S.parse_allterms(cand)
    return None, None, None


def test(base, rule):
    D, N = (int(x) for x in base.split("_"))
    meas = dict(S.MEASURED.get(base, {}))
    meas.update({k: v for k, v in S.HAUPT.get(base, {}).items() if k not in meas})
    forms, _, Mlev = load(base)
    if forms is None or not meas:
        return None
    dpath = os.path.join(S.WT, f"diff_{base}.txt")
    diff = S.parse_diff(dpath) if os.path.exists(dpath) else None
    if rule in ("diff", "both") and diff is None:
        return None
    idxs = sorted({i for f in forms.values() for i in f}, key=lambda t: (t[0], t[1]))
    free = []
    for i in idxs:
        mm = S.index_m(i, Mlev)
        lvl = (mm % N == 0)
        dif = diff is not None and len(diff.get(mm, [])) == 1
        if {"diff": dif, "level": lvl, "both": lvl and dif}[rule]:
            free.append((i, mm))
    rows, tg = [], []
    for k, f in forms.items():
        if k in meas:
            rows.append([f.get(i, F(0)) for i, _ in free])
            tg.append(F(meas[k]))
    if not rows:
        return None
    piv, M, ok = S.rref(rows, tg)
    if not ok:
        return dict(n=len(free), eqs=len(rows), ok=False)
    sol = [F(0)] * len(free)
    for r, c in enumerate(piv):
        sol[c] = M[r][len(free)]
    resid = len(free) - len(piv)
    pinned = []
    if resid == 0:
        pinned = [(free[k][0], free[k][1], sol[k]) for k in range(len(free))]
    return dict(n=len(free), eqs=len(rows), ok=True, spare=len(rows) - len(piv),
                resid=resid, pinned=pinned)


def main():
    bases = sys.argv[1:] or ["6_5", "6_7", "6_11", "6_17", "6_19", "6_23", "10_3", "10_11", "15_2"]
    tot = {r: [0, 0, 0] for r in RULES}      # [spare, refuted, bases]
    print(f"{'base':>7} {'rule':>6} {'#free':>6} {'#eqs':>5} {'result':>12} {'spare':>6} {'resid':>6}   weights (when unique)")
    for b in bases:
        for r in RULES:
            res = test(b, r)
            if res is None:
                continue
            tot[r][2] += 1
            if not res["ok"]:
                tot[r][1] += 1
                print(f"{b:>7} {r:>6} {res['n']:>6} {res['eqs']:>5} {'REFUTED':>12} {'-':>6} {'-':>6}")
                continue
            tot[r][0] += res["spare"]
            w = "  ".join(f"{'z' if i[0]=='0' else 'u'}{m}={v}" for i, m, v in res["pinned"])
            print(f"{b:>7} {r:>6} {res['n']:>6} {res['eqs']:>5} {'consistent':>12} "
                  f"{res['spare']:>6} {res['resid']:>6}   {w}")
        print()
    print("TOTALS")
    for r in RULES:
        print(f"   {r:>6}: {tot[r][0]:>3} spare conditions satisfied, "
              f"{tot[r][1]} refutation(s) over {tot[r][2]} base(s)")


if __name__ == "__main__":
    main()
