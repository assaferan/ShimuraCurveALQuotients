#!/usr/bin/env python3
"""Why do the PANEL forms obey the infinity-only formula?

cleancrit.py established: f holomorphic away from infinity => infinity-only
exact (319 monomials, 0 counterexamples). Genuine Borcherds inputs are not
obviously in that class, yet infinity-only reproduces all 41 measured
multipliers. This asks why.

CAUTION, and the reason this script exists rather than an ad-hoc command:
Magma wraps long FORMC lists mid-token ("<19,\n189180>"), so a naive regex
silently drops terms and yields garbage of order 1e7. Whitespace inside each
FORMC block must be normalised BEFORE matching. The parse is validated by
reproducing the paper's recorded multipliers: 15_2 -> 4,8,0,0,8,4,8,-4,4, and
twice the recorded values at 21_2 and 22_3.

FINDINGS.
 1. The infinity class contributes exactly 0 to the constant-term assembly at
    every form and base -- which is Thm 8.5(i) (rho-entry vanishes iff N | g,
    and g = M at infinity). Concerns the assembly, not the residue pairing.
 2. At 22_3 all nine panel forms are CLEAN, so infinity-only is exact there by
    the criterion, not by luck.
 3. At 15_2 and 21_2 each panel form carries exactly ONE off-infinity pole slot
    (n = 1/2 or 3/4 at 15_2; n = 1/4 at 21_2, two forms clean). So agreement
    for these reduces to a single coefficient per base vanishing.

PREDICTION, NOT VERIFIED: a_(1/2) = a_(3/4) = 0 at 15_2, a_(1/4) = 0 at 21_2.
Implied by 41/41 plus the corrected formula, but circular as a check, and the
monomial data cannot settle it -- no monomial has support in the determined
head together with one of those slots. Needs E at the 0 cusp, computed directly.
"""
import os
import re
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BASES = [('cusp6_15_2.out', 15, 2, 60), ('cusp6_21_2.out', 21, 2, 84),
         ('cusp6_22_3.out', 22, 3, 132)]
KNOWN = {15: [4, 8, 0, 0, 8, 4, 8, -4, 4]}


def parse(path):
    txt = open(path).read()
    pp, tg = {}, {}
    for ln in txt.split('\n'):
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
    forms = {}
    for k, body in re.findall(r'FORMC (-?\d+) \[(.*?)\]', txt, re.S):
        body = ' '.join(body.split())          # <-- the line-wrap fix
        forms[int(k)] = {int(a): F(b)
                         for a, b in re.findall(r'<(\d+), (-?\d+(?:/\d+)?)>', body)}
    return pp, tg, forms


def main():
    for fn, D, N, M in BASES:
        pp, tg, forms = parse(os.path.join(HERE, fn))
        print(f"=== X_0^{D}({N})  M={M} ===")
        for i, k in enumerate(sorted(forms)):
            cls = {}
            for r, c in forms[k].items():
                for g, v in tg.get(r, {}).items():
                    cls[g] = cls.get(g, 0.0) + float(c) * v
            inf = cls.get(M, 0.0)
            tot = sum(cls.values())
            off = {}
            for r, c in forms[k].items():
                for g, d in pp.get(r, {}).items():
                    if g == M:
                        continue
                    for n, v in d.items():
                        off[(g, n)] = off.get((g, n), 0.0) + float(c) * v
            big = {s: v for s, v in off.items() if abs(v) > 1e-6}
            chk = f" (known {KNOWN[D][i]})" if D in KNOWN else ""
            print(f"  form {k:3}: c_eta(0)={tot:+8.4f}{chk}  inf-class={inf:+.1e}"
                  f"  off-inf pole slots: {len(big)}"
                  + ("" if not big else "  " + ", ".join(f"n={n}:{v:+.3g}"
                                                         for (_, n), v in big.items())))
        print()


if __name__ == "__main__":
    main()
