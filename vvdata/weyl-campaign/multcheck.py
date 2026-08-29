#!/usr/bin/env python3
"""End-to-end check of Prop 9.15 against independently measured ground truth.

The chain verified elsewhere runs closed form -> genus theta combination.  This
runs it all the way out to the CM data:

    mult(f)  ==  (1/2) * sum_m c_oo(-m) * W(m),        W(m) = -a_E(m), Prop 9.15

with mult(f) the multiplier measured by kernel-consistency or by the Hauptmodul
relation (5.1) -- neither of which evaluates an Eisenstein series -- and c_oo the
infinity-side principal part read from the pp*_*.out dumps.

This is NOT the same check as the "4731 of 4731" sweep, which compares the closed
form against the theta combination it was derived from.  Here the right-hand side
is computed from the closed form alone and the left-hand side is measured.

Scope: the infinity block alone suffices on the nine-form panels.  Corollary 9.14
pairs over all cosets, and the 0-side genuinely enters on the larger monomial
probe space of 8.10 -- that case is not covered here.
"""
import glob
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from closedcoef import W  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    """Parse the principal-part dumps: FORM D N id mult / COEF D N id blk m exp c."""
    bases = {}
    for path in sorted(glob.glob(os.path.join(HERE, 'pp*_*.out'))):
        D = N = cur = None
        for ln in open(path):
            t = ln.split()
            if len(t) >= 3 and t[0] == 'BASEINFO':
                D, N = int(t[1]), int(t[2])
            elif len(t) >= 5 and t[0] == 'FORM':
                cur = t[3]
                mu = None if t[4] == '?' else F(t[4])
                e = bases.setdefault((D, N), {}).setdefault(
                    cur, {'mult': mu, 'oo': {}, '0': {}})
                if e['mult'] is None:
                    e['mult'] = mu
            elif len(t) >= 8 and t[0] == 'COEF':
                blk, m, c = t[4], int(t[5]), F(t[7])
                d = bases[(D, N)][cur]
                if blk in d:
                    d[blk][m] = d[blk].get(m, F(0)) + c
    return bases


def main():
    bases = load()
    print("mult(f) == (1/2) sum_m c_oo(-m) W(m),  W from the Prop 9.15 closed form")
    print("(forms carrying an independently measured multiplier)\n")
    tot = ok = 0
    for (D, N) in sorted(bases):
        rows = [(fid, d['mult'],
                 F(1, 2) * sum(c * W(D, N, m) for m, c in d['oo'].items()))
                for fid, d in bases[(D, N)].items() if d['mult'] is not None]
        if not rows:
            continue
        hit = sum(1 for _, mu, pr in rows if mu == pr)
        tot += len(rows)
        ok += hit
        print(f"  X_0^{D}({N}):  {hit}/{len(rows)}"
              + ("" if hit == len(rows) else "   *** FAIL ***"))
        for fid, mu, pr in rows:
            if mu != pr:
                print(f"       form {fid}: measured {mu}, predicted {pr}")
    print(f"\n  TOTAL: {ok}/{tot}")
    return 0 if ok == tot else 1


if __name__ == "__main__":
    sys.exit(main())
