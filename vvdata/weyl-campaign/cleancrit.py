#!/usr/bin/env python3
"""When is the infinity-only closed form exact?  When f is holomorphic away from infinity.

infonly.py showed the truncation

    mult(f) =? (1/2) sum_m c_oo(-m) W(m)

is not a general identity. This locates the criterion, over the three bases
with per-monomial ground truth (cusp6_*.out):

    if f has NO poles away from the infinity class, the truncation is EXACT.

Across 319 monomials there is not one counterexample. It is an iff at 21_2 and
22_3; at 15_2 exactly one monomial agrees anyway, by accidental cancellation,
so the implication is one-directional -- as it must be, since cancellation can
always produce coincidental agreement.

This explains the earlier counts: 22_3 has no monomial with poles away from
infinity (132/132); 21_2 has 80 such (80/101); 15_2 has 37, plus one accident
(38/86).

NOTE this does NOT yet explain the panel forms. Genuine Borcherds inputs DO
have 0-side principal parts (e.g. `COEF 15 2 -2 0 30 1/2 -2` in the pp dumps),
so they are not clean, yet infinity-only reproduces all 41 measured
multipliers over ten bases. Either those cancellations are systematic for
Borcherds inputs -- which would be the theorem to prove -- or 41 coincidences
have occurred. That is the open question gating production use.
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from closedcoef import W  # noqa: E402

BASES = [('cusp6_15_2.out', 15, 2, 60), ('cusp6_21_2.out', 21, 2, 84),
         ('cusp6_22_3.out', 22, 3, 132)]


def load(path):
    pp, tg = {}, {}
    for ln in open(path):
        t = ln.split()
        if not t:
            continue
        if t[0] == 'PP':
            r, g, Wg, e = int(t[1]), int(t[2][2:]), int(t[3][2:]), int(t[4][2:])
            if e < 0:
                pp.setdefault(r, {}).setdefault(g, {})[F(-e, 24 * Wg)] = float(t[5])
        elif t[0] == 'TG':
            r, g = int(t[1]), int(t[2][2:])
            tg.setdefault(r, {})[g] = float(t[3])
    return pp, tg


def main():
    tot = viol = 0
    for fn, D, N, M in BASES:
        pp, tg = load(os.path.join(HERE, fn))
        rs = sorted(tg)
        tt = ft = tf = 0
        for r in rs:
            truth = sum(tg[r].values())
            pred = sum(c * float(W(D, N, int(m)))
                       for m, c in pp.get(r, {}).get(M, {}).items() if m.denominator == 1)
            agree = abs(truth - pred) < 1e-9
            clean = not any(g != M for g in pp.get(r, {}))
            if clean and agree:
                tt += 1
            elif clean:
                tf += 1
            elif agree:
                ft += 1
        tot += len(rs)
        viol += tf
        print(f"X_0^{D}({N}):  {len(rs):3} monomials   clean&exact={tt:3}"
              f"   exact-anyway={ft}   clean&WRONG={tf}   iff={tf == 0 and ft == 0}")
    print(f"\n  {tot} monomials, {viol} counterexamples to"
          f' "holomorphic away from infinity => infinity-only exact"')


if __name__ == "__main__":
    main()
