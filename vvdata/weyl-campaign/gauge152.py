#!/usr/bin/env python3
"""Why Table 5.6 and the Prop 9.15 coefficients disagree index by index.

Both reproduce all nine measured multipliers of X_0^15(2), yet they differ at
five of the six indices the panel touches.  The reason is that the panel's
infinity-side principal-part matrix is rank deficient: 9 forms x 6 indices of
rank 4, so an infinity-only functional is pinned only modulo a 2-dimensional
space, and the two vectors sit in the same coset of it.

Consequence: the N | m support rule is a choice of representative inside that
ambiguity, not a property of the correction -- the Prop 9.15 coefficients are
nonzero at m = 3, 7, 13, 15, 27.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from closedcoef import W          # noqa: E402
from multcheck import load        # noqa: E402

TABLE_A = {2: F(-1), 10: F(1), 30: F(0)}   # Table 5.6, X_0^15(2)


def rank(rows, ncols):
    M = [r[:] for r in rows]
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r


def main():
    forms = load()[(15, 2)]
    idx = sorted({m for d in forms.values() for m in d['oo']})
    print("X_0^15(2): infinity-side principal parts, indices", idx, "\n")
    print(" form   measured   Table 5.6   (1/2)W  [Prop 9.15]")
    for fid, d in sorted(forms.items(), key=lambda x: int(x[0])):
        a = sum(c * TABLE_A.get(m, F(0)) for m, c in d['oo'].items())
        w = F(1, 2) * sum(c * W(15, 2, m) for m, c in d['oo'].items())
        mark = "agree" if a == w == d['mult'] else "DIFFER"
        print(f"  {fid:>3}   {str(d['mult']):>8}   {str(a):>9}   {str(w):>6}   {mark}")

    print("\nthe two representatives, index by index:")
    print("    m  :  A_m     (1/2)W(m)")
    for m in idx:
        print(f"  {m:3}  :  {str(TABLE_A.get(m, 0)):>4}  {str(F(1, 2) * W(15, 2, m)):>10}")

    rows = [[d['oo'].get(m, F(0)) for m in idx] for d in forms.values()]
    r = rank(rows, len(idx))
    print(f"\n  panel matrix {len(rows)} x {len(idx)}, rank {r}"
          f"  =>  {len(idx) - r} free directions")
    print("  both vectors lie in the same coset of that space, which is why they")
    print("  agree on every panel form and disagree at almost every index.")

    off = [m for m in range(1, 41) if m % 2 and W(15, 2, m) != 0]
    print(f"\n  Prop 9.15 nonzero at N nmid m, m < 40:  {off}")


if __name__ == "__main__":
    main()
