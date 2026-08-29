#!/usr/bin/env python3
"""Empirical recovery of E's coefficients at the 0 cusp, on X_0^15(2).

infonly.py showed the infinity-only truncation is not a general identity
(38/86 on the eta monomials). Corollary 9.14 pairs over ALL cosets, so the
missing piece is a pairing against E's coefficients at the other cusps. At
15_2 only two cusp classes occur, g = 60 (infinity, width 1) and g = 1, so the
correction is a single extra block:

    mult(f) = (1/2) [ sum_m c_oo(-m) W(m) + sum_n c_0(-n) a_n ]

Fitting delta = truth - infinity-only against the 45 0-cusp slots over 86
monomials (overdetermined by 41) gives a well-conditioned HEAD whose
coefficients snap to small integers, and an ill-conditioned tail where few
monomials have poles that deep. Restricting to the 49 monomials whose 0-cusp
support lies inside the head, the corrected formula is exact: 49/49.

Recovered (exponent n at the 0 cusp -> a_n), all multiples of 4:

    60n :  1   2   3   4   5   6   7   8   9  10  11  14
    a_n :  4  -4   0   8   0  -4  -8   4   4   0  -8  -8

This does not yet give a closed form at the 0 cusp. It fixes the target that
any such formula must reproduce, and confirms Cor 9.14's shape numerically.
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from closedcoef import W  # noqa: E402

HEAD = {F(1, 60): 4, F(1, 30): -4, F(1, 20): 0, F(1, 15): 8,
        F(1, 12): 0, F(1, 10): -4, F(7, 60): -8, F(2, 15): 4,
        F(3, 20): 4, F(1, 6): 0, F(11, 60): -8, F(7, 30): -8}


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
    pp, tg = load(os.path.join(HERE, 'cusp6_15_2.out'))
    rs = sorted(tg)
    inside = [r for r in rs if all(s in HEAD for s in pp.get(r, {}).get(1, {}))]
    ok, bad = 0, []
    for r in inside:
        truth = sum(tg[r].values())
        pred = sum(c * float(W(15, 2, int(m)))
                   for m, c in pp.get(r, {}).get(60, {}).items() if m.denominator == 1)
        corr = sum(HEAD[s] * v for s, v in pp.get(r, {}).get(1, {}).items())
        if abs(truth - (pred + corr)) < 1e-9:
            ok += 1
        else:
            bad.append((r, truth, pred + corr))
    print("X_0^15(2): mult = (1/2)[ oo-pairing + 0-cusp pairing ]")
    print(f"  monomials with 0-cusp support inside the determined head: {len(inside)}")
    print(f"  reproduced exactly: {ok}/{len(inside)}")
    for r, t_, p in bad[:5]:
        print(f"    mono {r}: truth {t_:.6f} vs {p:.6f}")
    print("  (infinity-only alone: 38/86 overall -- see infonly.py)")


if __name__ == "__main__":
    main()
