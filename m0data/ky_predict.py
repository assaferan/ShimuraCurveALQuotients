#!/usr/bin/env python3
"""
Predict the m=0 weights u_m from Kudla-Yang, "Eisenstein series for SL(2)" (Sci. China Math. 53 (2010)),
and test against the measured ground truth.

LOCAL FACTORS.  KY (9.1) -- B quaternion of discriminant D, R an Eichler order of level N, (N,D)=1 --
gives, with X = p^(-s) and k = ord_p(m),

    p not | DN :  (1-X)^-1 [1 - X^(k+1)]                        -> at X=1:  k + 1
    p | D      :  (1-X)^-1 [1 - pX + pX^(k+1) - X^(k+2)]        -> at X=1:  2 - k(p-1)
    p | N      :  (1-X)^-1 [1 - 2X + pX - pX^(k+1) + X^(k+2)]   -> at X=1:  k(p-1)

The p | N factor VANISHES iff p does not divide m -- this is the derivation of the level support rule.

KY (2.17) is the ternary/general version at p | D, in terms of the CONDUCTOR exponent k = ord_p(c)
where 4*kappa*m = d*c^2, and v_p = 1, 0, -1 as Q(sqrt(kappa m)) is split, ramified, inert at p:

    b_p(1) = [ (1-v)(1-p^2) + p^(k+2) - p^(2k+2) ] / (1-p)

TEST.  The pipeline's own machinery for the non-ramified, non-level places is Table-45 validated, so we
keep the dump's cond_half * (h/w) * (en/ed) * prod_{p not | DN} W_p and swap in the KY factors at p | D
and p = N.  Each base is allowed ONE overall constant (the global normalisation KY carries in C^+(l,D)
etc.), so the test is whether the predicted vector is PROPORTIONAL to the measured one.

usage: ky_predict.py [--pD=ky|code] [base ...]
"""
import sys, os, re
from fractions import Fraction as F
import support_fit as S
import joint_fit as J
from rule_test import test

PD_MODE = "ky"


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


def ordp(m, p):
    v = 0
    while m % p == 0:
        m //= p
        v += 1
    return v


def kronecker(d, p):
    """chi_d(p) for the fundamental discriminant d: +1 split, -1 inert, 0 ramified."""
    if d % p == 0:
        return 0
    if p == 2:
        return 1 if d % 8 == 1 else (-1 if d % 8 == 5 else 0)
    r = pow(d % p, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def b_ramified(p, k, v):
    """KY (2.17) evaluated at X = 1;  k = ord_p(conductor), v = v_p in {1,0,-1}."""
    num = (1 - v) * (1 - p ** 2) + p ** (k + 2) - p ** (2 * k + 2)
    return F(num, 1 - p)


def parse_dump(base):
    """per (block,index): c, cond_half, h/w, en/ed, {p: W_p}, d"""
    path = None
    for cand in (os.path.join(S.WT, f"allterms_{base}.txt"),
                 os.path.join(S.WT, "atdump", f"allterms_{base}.txt")):
        if os.path.exists(cand) and "DONE" in open(cand).read():
            path = cand
    out, Mlev = {}, None
    for line in open(path):
        mm = re.search(r"M = (\d+)", line)
        if mm and Mlev is None:
            Mlev = int(mm.group(1))
        mm = re.match(r"\s*oo m=(\d+)\s+c=(-?\d+)\s+contrib=\S+\s+\| d=(-?\d+) h=(\d+) w=(\d+) "
                      r"cond/2=(\S+) en/ed=(\S+) g=(\S+) locs\(p,W,v_p\(r\)\)=\[(.*)\]", line)
        if mm:
            idx = ("oo", int(mm.group(1)))
            out[idx] = dict(
                d=int(mm.group(3)), hw=F(int(mm.group(4)), int(mm.group(5))),
                cond=F(mm.group(6)), ened=F(mm.group(7)),
                W={int(p): F(w.strip()) for p, w, _ in
                   re.findall(r"<\s*(\d+),\s*([^,]+),\s*(-?\d+)\s*>", mm.group(9))})
    return out, Mlev


def predict(base):
    D, N = (int(x) for x in base.split("_"))
    res = test(base, "level")
    if not res or not res.get("pinned"):
        return None
    terms, Mlev = parse_dump(base)
    ram = set(pdivs(D))
    rows = []
    for idx, m, u_true in res["pinned"]:
        t = terms.get(idx)
        if t is None:
            continue
        # conductor c from 4m = |d| c^2
        c2 = F(4 * m, abs(t["d"]))
        c = F(c2.numerator, 1) ** 0  # placeholder; compute integer sqrt below
        num, den = c2.numerator, c2.denominator
        r = (num // den) if den == 1 else None
        c_int = None
        if r is not None:
            s0 = int(round(r ** 0.5))
            if s0 * s0 == r:
                c_int = s0
        pred = t["cond"] * t["hw"] * t["ened"]
        for p, w in t["W"].items():
            if p in ram:
                if PD_MODE == "ky" and c_int is not None:
                    pred *= b_ramified(p, ordp(c_int, p), kronecker(t["d"], p))
                else:
                    pred *= w
            elif p == N:
                pred *= F(ordp(m, N) * (N - 1))
            else:
                pred *= w
        rows.append((m, u_true, pred))
    return D, N, rows


def main():
    global PD_MODE
    args = sys.argv[1:]
    while args and args[0].startswith("--"):
        a = args.pop(0)
        if a.startswith("--pD="):
            PD_MODE = a.split("=", 1)[1]
    bases = args or ["6_5", "6_17", "6_19", "10_3", "22_7", "15_2"]
    print(f"[p|D factor from: {PD_MODE}]")
    nprop = ntot = 0
    for base in bases:
        r = predict(base)
        if not r:
            continue
        D, N, rows = r
        ratios = []
        print(f"--- {base}  (D={D}, N={N})")
        for m, u, pred in rows:
            rat = (u / pred) if pred != 0 else None
            ratios.append(rat)
            print(f"    m={m:<5} u_true={str(u):>6}  predicted={str(pred):>10}  ratio={str(rat):>10}")
        good = [x for x in ratios if x is not None]
        ok = len(set(good)) == 1 and good
        # zeros must line up too
        zeros_ok = all((u == 0) == (pred == 0) for _, u, pred in rows)
        ntot += 1
        nprop += 1 if (ok and zeros_ok) else 0
        print(f"    => proportional: {ok and zeros_ok}"
              + (f"   (constant = {good[0]})" if ok and zeros_ok else "")
              + ("" if zeros_ok else "   [ZERO PATTERN MISMATCH]"))
    print(f"\n{nprop}/{ntot} base(s) reproduced up to one overall constant")


if __name__ == "__main__":
    main()
