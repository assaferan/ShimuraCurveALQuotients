#!/usr/bin/env python3
"""Is the infinity-only closed form a general identity?  NO.

    mult(f) =? (1/2) sum_m c_oo(-m) W(m)          [infinity block only]

holds on every genuine Borcherds input tested (41/41 over ten bases,
multcheck.py + 22_5), but it is NOT an identity for arbitrary weakly
holomorphic input of the right character.  Corollary 9.14 pairs over ALL
cosets, and this script shows the infinity-only truncation failing on the 86
eta monomials of 15_2, which are individually weakly holomorphic of weight 1/2
with the panel character.

Ground truth comes from cusp6_15_2.out:
  TG <r> g=<class> re im   -- class sums, ALREADY summed over the cosets in the
                              class, so c_eta*(0) for monomial r is sum_g TG.
  PP <r> g=60 W=1 e=<E> v  -- principal part at the INFINITY class (lower-left
                              c = 0 mod M, width 1); exponent is E/(24*W_g).
Validated by the file's own SELFC lines, which reproduce 2x the nine measured
panel multipliers exactly.
"""
import os
import sys
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from closedcoef import W  # noqa: E402


def load(path):
    pp, tg = {}, {}
    for ln in open(path):
        t = ln.split()
        if not t:
            continue
        if t[0] == 'PP':
            r, g, Wg, e = int(t[1]), int(t[2][2:]), int(t[3][2:]), int(t[4][2:])
            if g == 60 and e < 0:
                assert e % (24 * Wg) == 0
                pp.setdefault(r, {})[-e // (24 * Wg)] = Decimal(t[5])
        elif t[0] == 'TG':
            r = int(t[1])
            tg[r] = tg.get(r, Decimal(0)) + Decimal(t[3])
    return pp, tg


def main():
    pp, tg = load(os.path.join(HERE, 'cusp6_15_2.out'))
    rs = sorted(tg)
    agree, dis = 0, []
    for r in rs:
        truth = tg[r] / 2
        pred = Decimal(0)
        for m, c in pp.get(r, {}).items():
            w = W(15, 2, m)
            pred += Decimal(int(c)) * Decimal(w.numerator) / Decimal(w.denominator)
        pred /= 2
        if abs(truth - pred) < Decimal('1e-20'):
            agree += 1
        else:
            dis.append((r, truth, pred))
    print(f"X_0^15(2), {len(rs)} eta monomials")
    print(f"  infinity-only closed form matches ground truth: {agree}/{len(rs)}")
    print(f"  disagreements: {len(dis)}\n")
    for r, t_, p in dis[:6]:
        print(f"    mono {r:3}: truth {t_:+.6f}  pred {p:+.6f}  diff {t_ - p:+.6f}")
    print("\n  => the infinity-only truncation is NOT a general identity.")
    print("     It holds on genuine Borcherds inputs (41/41 over ten bases) but")
    print("     fails on individual monomials, whose 0-side principal parts the")
    print("     pairing sees and the truncation ignores.")


if __name__ == "__main__":
    main()
