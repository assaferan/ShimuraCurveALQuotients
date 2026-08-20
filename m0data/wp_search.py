#!/usr/bin/env python3
"""
Settle the m=0 VALUE formula mechanically instead of by hand.

We have, for our actual lattice, the exact local Whittaker POLYNOMIAL W_p(m, x) at every relevant p
(dumped by wpoly_dump.m via the new WhittakerPolynomial intrinsic).  The pipeline is known to compute
these correctly -- an independent representation count reproduces them.  What is NOT settled is which
FUNCTIONAL of each polynomial enters the m = 0 term, and that is exactly a normalisation question:
KY evaluate b_p at X = 1/p (s = l-1 = 1/2 for weight 3/2), while the code's scaled variable puts its
critical point at x = 1, and the log p terms come from the s-derivative.

So: search the small space of natural functionals, independently for the three prime classes
(good / p | D / p = N), times the natural prefactor choices, and demand PROPORTIONALITY to the measured
weights across every base (one constant per base).  Exact rational arithmetic, so a hit is a hit.

usage: wp_search.py [base ...]
"""
import sys, os, re, itertools
from fractions import Fraction as F
import support_fit as S
from rule_test import test

WTE = "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/4f908ae6-6cc4-4e2c-8456-b6725f6198b0/scratchpad/e2e-6-5/m0data"


def pdivs(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def parse_wp(base):
    """(r, p) -> polynomial coefficient list (rationals), for the oo block (eta = 0)."""
    path = os.path.join(WTE, f"wp_{base}.txt")
    if not os.path.exists(path) or "DONE" not in open(path).read():
        return None
    out = {}
    for line in open(path):
        mm = re.match(r"\[WP\] form=(-?\d+) oo r=(\S+) c=(-?\d+) p=(\d+) deg=(-?\d+) coeffs=\[(.*?)\]", line)
        if mm:
            r = F(mm.group(2))
            p = int(mm.group(4))
            co = [F(x.strip()) for x in mm.group(6).split(",") if x.strip()]
            out[(r, p)] = co
    return out


def ev(co, x):
    v = F(0)
    for i, a in enumerate(co):
        v += a * x ** i
    return v


def dev(co, x):
    v = F(0)
    for i, a in enumerate(co):
        if i:
            v += i * a * x ** (i - 1)
    return v


FUNCS = {
    "W(1)":     lambda co, p: ev(co, F(1)),
    "W(1/p)":   lambda co, p: ev(co, F(1, p)),
    "W'(1)":    lambda co, p: dev(co, F(1)),
    "W'(1/p)":  lambda co, p: dev(co, F(1, p)),
    "W(1)|W'":  lambda co, p: ev(co, F(1)) if ev(co, F(1)) != 0 else dev(co, F(1)),
    "W'(1)/W(1)": lambda co, p: (dev(co, F(1)) / ev(co, F(1))) if ev(co, F(1)) != 0 else None,
}


def parse_terms(base):
    path = None
    for cand in (os.path.join(S.WT, f"allterms_{base}.txt"),
                 os.path.join(S.WT, "atdump", f"allterms_{base}.txt")):
        if os.path.exists(cand) and "DONE" in open(cand).read():
            path = cand
    out = {}
    for line in open(path):
        mm = re.match(r"\s*oo m=(\d+)\s+c=(-?\d+)\s+contrib=\S+\s+\| d=(-?\d+) h=(\d+) w=(\d+) "
                      r"cond/2=(\S+) en/ed=(\S+) ", line)
        if mm:
            out[("oo", int(mm.group(1)))] = dict(
                h=int(mm.group(4)), w=int(mm.group(5)),
                cond=F(mm.group(6)), ened=F(mm.group(7)))
    return out


def build(bases):
    data = []
    for b in bases:
        D, N = (int(x) for x in b.split("_"))
        res = test(b, "level")
        wp = parse_wp(b)
        if not res or not res.get("pinned") or wp is None:
            continue
        terms = parse_terms(b)
        rows = []
        for idx, m, u in res["pinned"]:
            t = terms.get(idx)
            if t is None:
                continue
            polys = {p: co for (r, p), co in wp.items() if r == F(m)}
            if not polys:
                continue
            rows.append(dict(m=m, u=u, t=t, polys=polys))
        if rows:
            data.append((b, D, N, rows))
    return data


PREFAC = {
    "h/w":              lambda t: F(t["h"], t["w"]),
    "cond*h/w":         lambda t: t["cond"] * F(t["h"], t["w"]),
    "h/w*en/ed":        lambda t: F(t["h"], t["w"]) * t["ened"],
    "cond*h/w*en/ed":   lambda t: t["cond"] * F(t["h"], t["w"]) * t["ened"],
}


def main():
    bases = sys.argv[1:] or ["6_5", "6_17", "6_19", "10_3", "22_7", "15_2"]
    data = build(bases)
    print(f"bases with both a unique level-rule fit and a Whittaker dump: "
          f"{[b for b, _, _, _ in data]}")
    if not data:
        return
    hits = []
    for pfname, pf in PREFAC.items():
        for gname, gf in FUNCS.items():
            for dname, df in FUNCS.items():
                for lname, lf in FUNCS.items():
                    okbases = 0
                    for b, D, N, rows in data:
                        ram = set(pdivs(D))
                        rats, bad = [], False
                        for row in rows:
                            val = pf(row["t"])
                            for p, co in sorted(row["polys"].items()):
                                f = df if p in ram else (lf if p == N else gf)
                                x = f(co, p)
                                if x is None:
                                    bad = True
                                    break
                                val *= x
                            if bad:
                                break
                            rats.append((row["u"], val))
                        if bad:
                            break
                        # zeros must correspond, and nonzero ratios must agree
                        if any((u == 0) != (v == 0) for u, v in rats):
                            continue
                        nz = {u / v for u, v in rats if v != 0}
                        if len(nz) <= 1:
                            okbases += 1
                    if okbases == len(data):
                        hits.append((pfname, gname, dname, lname))
    print(f"\ntried {len(PREFAC)*len(FUNCS)**3} combinations over {len(data)} base(s)")
    if hits:
        print("MATCHES (proportional on every base):")
        for h in hits:
            print(f"   prefactor={h[0]:<16} good={h[1]:<12} p|D={h[2]:<12} p=N={h[3]}")
    else:
        print("no combination reproduces every base")
        # report the best partial
        best = None
        for pfname, pf in PREFAC.items():
            for gname, gf in FUNCS.items():
                for dname, df in FUNCS.items():
                    for lname, lf in FUNCS.items():
                        okb = []
                        for b, D, N, rows in data:
                            ram = set(pdivs(D))
                            rats, bad = [], False
                            for row in rows:
                                val = pf(row["t"])
                                for p, co in sorted(row["polys"].items()):
                                    f = df if p in ram else (lf if p == N else gf)
                                    x = f(co, p)
                                    if x is None:
                                        bad = True
                                        break
                                    val *= x
                                if bad:
                                    break
                                rats.append((row["u"], val))
                            if bad:
                                break
                            if any((u == 0) != (v == 0) for u, v in rats):
                                continue
                            if len({u / v for u, v in rats if v != 0}) <= 1:
                                okb.append(b)
                        if best is None or len(okb) > len(best[1]):
                            best = ((pfname, gname, dname, lname), okb)
        print(f"best partial: {best[0]} works on {len(best[1])}/{len(data)}: {best[1]}")


if __name__ == "__main__":
    main()
