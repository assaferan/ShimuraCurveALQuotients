#!/usr/bin/env python3
"""
Build the m=0 weights from Kudla-Yang Propositions 5.3/5.4 -- the sections that treat OUR lattice:
the ternary trace-zero lattice of an Eichler order, eq (5.3)  L_e = {(b a; c -b): a,b in Z_p, c in p^e Z_p}.

Prop 5.3 (e = 0):  W_p(s + 1/2, m, mu) = char(...) * L_p(s+1, chi_km) * b_p(km, s+1; D) *
                                          { zeta_p(2s+2)^-1  if B is split at p
                                          { 1                if B is ramified at p
Prop 5.4 (e > 0, B split, i.e. p | N), with c(l) = min(l/2, e - l/2), k = ord_p(conductor), X = p^-s:
      p | d      :  1 + (1-1/p) sum_{1<=l<=k} p^{c(2l)} X^{2l} - p^{c(2k+2)-1} X^{2k+2}
      p not | d  :  1 + (1-1/p) sum_{1<=l<=k} p^{c(2l)} X^{2l} + chi(p) p^{c(2k+1)-1/2} X^{2k+1}

EVALUATION POINT.  The special value for weight l = 3/2 is s = l - 1 = 1/2 in the Eisenstein variable,
which is the argument of W_p; since Prop 5.3/5.4 write the argument as s + 1/2, that is s = 0, so
W_p is taken at X = 1 -- but b_p, whose argument is s + 1, is then at X = 1/p.  (This is consistent with
Prop 5.4 at X = 1 reproducing the pipeline's level-prime factor exactly, which was verified.)

NORMALISATION.  The pipeline computes (h/w) * (en/ed) * prod_{Sc} W_p, where en/ed = prod_{Sc}
(1-chi(p)/p)/(1-1/p^2) strips the Euler factors at Sc from the class-number L-value.  KY's W_p already
CONTAINS L_p(s+1,chi) = (1-chi(p)/p)^-1, so those cancel; what survives per prime is coded below.
Variants are provided for the genuinely ambiguous choices rather than guessed silently.

usage: ky53.py [base ...]
"""
import sys, os, re, itertools
from fractions import Fraction as F
import support_fit as S
from rule_test import test
from ky_clean import pdivs, ordp, kron, parse_dump


def b_216(p, k, v, X):
    """KY (2.16) as a rational function, at the given X."""
    num = 1 - v * X + p ** k * v * X ** (1 + 2 * k) - p ** (k + 1) * X ** (2 * k + 2)
    return num / (1 - p * X ** 2)


def b_217(p, k, v, X):
    """KY (2.17) at the given X."""
    num = ((1 - v * X) * (1 - p ** 2 * X ** 2) - v * p ** (k + 1) * X ** (2 * k + 1)
           + p ** (k + 2) * X ** (2 * k + 2) + v * p ** (k + 1) * X ** (2 * k + 3)
           - p ** (2 * k + 2) * X ** (2 * k + 4))
    return num / (1 - p * X ** 2)


def w_54(p, k, v, e=1):
    """KY Prop 5.4 at X = 1 (the special value), e = 1.  v = chi_d(p); p | d  <=>  v == 0."""
    tot = F(1)
    for l in range(1, k + 1):
        c = min(F(l), F(e - l))          # c(2l) = min(l, e-l)
        tot += (1 - F(1, p)) * (F(p) ** c if c >= 0 else F(1, p ** int(-c)))
    if v == 0:                            # p | d
        c = min(F(k + 1), F(e - k - 1))
        tot -= (F(p) ** c if c >= 0 else F(1, p ** int(-c))) * F(1, p)
    else:                                 # p does not divide d
        # chi(p) * p^{c(2k+1) - 1/2} * X^{2k+1}; at X = 1 and e = 1 the half powers cancel for k = 0
        c2 = min(F(2 * k + 1, 2), F(e) - F(2 * k + 1, 2))
        ex = c2 - F(1, 2)
        assert ex.denominator == 1, (p, k, v, ex)
        tot += v * (F(p) ** ex if ex >= 0 else F(1, p ** int(-ex)))
    return tot


def predict(base, variant):
    D, N = (int(x) for x in base.split("_"))
    res = test(base, "level")
    if not res or not res.get("pinned"):
        return None
    terms = parse_dump(base)
    ram = set(pdivs(D))
    out = []
    for idx, m, u in res["pinned"]:
        t = terms.get(idx)
        if t is None:
            continue
        c = 2 * t["cond"]
        if c.denominator != 1:
            return None
        c, d = int(c), t["d"]
        val = F(t["h"], t["w"])
        if variant["cond"]:
            val *= t["cond"]
        for p in sorted(set(pdivs(2 * c * D * N))):
            k, v = ordp(c, p), kron(d, p)
            X = F(1, p)
            euler = F(1 - F(v, p), 1 - F(1, p * p))     # the pipeline's en/ed factor at p
            if p in ram:
                val *= euler * (1 - F(v, p)) ** -1 * b_217(p, k, v, X)
            elif p == N:
                w = w_54(p, k, v)
                val *= euler * ((1 - F(v, p)) ** -1 if variant["lp_at_N"] else F(1)) * w
            else:
                val *= euler * (1 - F(v, p)) ** -1 * b_216(p, k, v, X) * (1 - F(1, p * p)) ** -1 \
                    if variant["zeta_inv"] else euler * (1 - F(v, p)) ** -1 * b_216(p, k, v, X) * (1 - F(1, p * p))
        out.append((m, u, val))
    return out


def main():
    bases = sys.argv[1:] or ["6_5", "6_17", "6_19", "10_3", "22_7", "15_2"]
    best = None
    for variant in [dict(cond=a, lp_at_N=b, zeta_inv=c)
                    for a in (False, True) for b in (False, True) for c in (False, True)]:
        okb, detail = [], {}
        for base in bases:
            rows = predict(base, variant)
            if not rows:
                continue
            if any((u == 0) != (v == 0) for _, u, v in rows):
                continue
            nz = {u / v for _, u, v in rows if v != 0}
            detail[base] = rows
            if len(nz) <= 1:
                okb.append(base)
        if best is None or len(okb) > len(best[1]):
            best = (variant, okb, detail)
        if len(okb) == len(bases):
            print(f"FULL MATCH with variant {variant}")
            break
    variant, okb, detail = best
    print(f"best variant: {variant}")
    print(f"proportional on {len(okb)}/{len(bases)}: {okb}\n")
    for base, rows in detail.items():
        print(f"--- {base}")
        for m, u, v in rows:
            print(f"    m={m:<5} u_true={str(u):>6}  pred={str(v):>16}  ratio="
                  f"{str(u / v) if v != 0 else 'n/a':>14}")


if __name__ == "__main__":
    main()
