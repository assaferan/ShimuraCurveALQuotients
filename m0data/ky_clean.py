#!/usr/bin/env python3
"""
Assemble the Kudla-Yang ternary Eisenstein coefficient CLEANLY (no splicing into the code path -- that
double-counted, see ky_predict.py), then (a) test it and (b) if it fails, fit ONLY the p | D factor with
everything else pinned by theory.

KY, "Eisenstein series for SL(2)", Corollary 8.4 / (8.1): for the ternary lattice attached to a maximal
order in a quaternion algebra of discriminant D, the weight-l Eisenstein coefficients are

    H(l, m; D) = L(3/2 - l, chi_{kappa m}) * b(kappa m, 3/2 - l; D)

At l = 3/2 this is L(0, chi_d) * b(kappa m, 0; D), and for an imaginary quadratic discriminant d,
L(0, chi_d) = 2 h(d)/w(d) -- the class number, exactly the pipeline's (h/w).  The remaining factor is
b = prod_p b_p with, writing 4*kappa*m = d c^2, k = ord_p(c), X = p^-s at s = 0 (X = 1), and
v = chi_d(p) in {+1 split, 0 ramified, -1 inert}:

    p not | D   (2.16):   [ 1 - v + p^k v - p^(k+1) ] / (1 - p)          -- equals 1 when k = 0
    p | D       (2.17):   [ (1-v)(1-p^2) + p^(k+2) - p^(2k+2) ] / (1 - p)

NOTE the pipeline's cond_half and (en/ed) must NOT be reused here: in KY those are already inside b_p
at p | c.  That was the double-count.

Our lattice is ternary WITH Eichler level N, which KY section 8 does not cover (it is maximal-order
only).  The level factor is taken from (9.1)'s p | N case (quaternary): ord_p(m) * (p-1).  Whether that
is right in the ternary case is exactly what this script tests.

usage: ky_clean.py [--fit-pD] [base ...]
"""
import sys, os, re
from fractions import Fraction as F
import support_fit as S
from rule_test import test


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


def kron(d, p):
    """chi_d(p): +1 split, -1 inert, 0 ramified."""
    if d % p == 0:
        return 0
    if p == 2:
        return 1 if d % 8 == 1 else -1
    return 1 if pow(d % p, (p - 1) // 2, p) == 1 else -1


def b_good(p, k, v):
    """KY (2.16) at X = 1."""
    return F(1 - v + p ** k * v - p ** (k + 1), 1 - p)


def b_ram(p, k, v):
    """KY (2.17) at X = 1."""
    return F((1 - v) * (1 - p ** 2) + p ** (k + 2) - p ** (2 * k + 2), 1 - p)


def b_level(p, m):
    """KY (9.1), the p | N case, at X = 1: ord_p(m) * (p-1)."""
    return F(ordp(m, p) * (p - 1))


def parse_dump(base):
    path = None
    for cand in (os.path.join(S.WT, f"allterms_{base}.txt"),
                 os.path.join(S.WT, "atdump", f"allterms_{base}.txt")):
        if os.path.exists(cand) and "DONE" in open(cand).read():
            path = cand
    out = {}
    for line in open(path):
        mm = re.match(r"\s*oo m=(\d+)\s+c=(-?\d+)\s+contrib=\S+\s+\| d=(-?\d+) h=(\d+) w=(\d+) "
                      r"cond/2=(\S+) ", line)
        if mm:
            out[("oo", int(mm.group(1)))] = dict(
                d=int(mm.group(3)), h=int(mm.group(4)), w=int(mm.group(5)),
                cond=F(mm.group(6)))
    return out


def rows_for(base):
    D, N = (int(x) for x in base.split("_"))
    res = test(base, "level")
    if not res or not res.get("pinned"):
        return None
    terms = parse_dump(base)
    ram, out = set(pdivs(D)), []
    for idx, m, u_true in res["pinned"]:
        t = terms.get(idx)
        if t is None:
            continue
        c = 2 * t["cond"]
        assert c.denominator == 1, (base, m, c)
        c = int(c)
        d = t["d"]
        pieces = {}
        val = F(t["h"], t["w"])
        for p in sorted(set(pdivs(2 * c * D * N))):
            if p in ram:
                pieces[(p, ordp(c, p), kron(d, p))] = pieces.get((p, ordp(c, p), kron(d, p)), 0) + 1
            elif p == N:
                val *= b_level(p, m)
            else:
                val *= b_good(p, ordp(c, p), kron(d, p))
        out.append(dict(m=m, u=u_true, base=val, ram=pieces, c=c, d=d))
    return D, N, out


def main():
    args = sys.argv[1:]
    fit = "--fit-pD" in args
    args = [a for a in args if not a.startswith("--")]
    bases = args or ["6_5", "6_17", "6_19", "10_3", "22_7", "15_2", "14_11", "10_11"]
    print("=== direct test: KY (2.17) at p | D, (9.1) level factor, L(0,chi)=2h/w ===")
    allrows = []
    for b in bases:
        r = rows_for(b)
        if not r:
            continue
        D, N, rows = r
        print(f"--- {b} (D={D}, N={N})")
        rats = []
        for row in rows:
            pred = row["base"]
            for (p, k, v), mult in row["ram"].items():
                pred *= b_ram(p, k, v) ** mult
            rats.append(None if pred == 0 else row["u"] / pred)
            print(f"    m={row['m']:<5} c={row['c']:<3} d={row['d']:<6} u={str(row['u']):>6} "
                  f"pred={str(pred):>12} ratio={str(rats[-1]):>12}")
        good = [x for x in rats if x is not None]
        zok = all((row["u"] == 0) == (pr == 0) for row, pr in
                  zip(rows, [row["base"] * _prod(row["ram"]) for row in rows]))
        print(f"    => proportional: {len(set(good)) == 1 and bool(good) and zok}")
        allrows.append((b, D, N, rows))
    if fit:
        fit_pD(allrows)


def _prod(ram):
    v = F(1)
    for (p, k, s), mult in ram.items():
        v *= b_ram(p, k, s) ** mult
    return v


def fit_pD(allrows):
    """Everything pinned except the p|D factor: solve for it as an unknown per (p, k, v)."""
    print("\n=== fit: ONLY the p|D factor unknown, indexed by (p, ord_p(c), chi_d(p)) ===")
    print("    (multiplicative -> take one unknown per class; equations are u = C_base * base * prod A)")
    cols, colidx, rows_out = [], {}, []

    def col(k):
        if k not in colidx:
            colidx[k] = len(cols)
            cols.append(k)
        return colidx[k]

    # log-linear is wrong (signs/zeros); instead solve the BILINEAR system exactly as in joint_fit:
    # u * t_base - base * prod(A) = 0 is not linear in A when several ram primes occur, so restrict to
    # rows with exactly ONE ramified class occurrence (true for all our data: D is a product of 2 primes
    # and the conductor is small) and report which rows are excluded.
    eqs, excluded = [], []
    for b, D, N, rows in allrows:
        tb = col(("t", b))
        for row in rows:
            key = tuple(sorted((k, v) for k, v in row["ram"].items()))
            e = {tb: row["u"], col(key): -row["base"]}
            eqs.append(e)
    dense = [[F(e.get(j, 0)) for j in range(len(cols))] for e in eqs]
    piv, M, _ = S.rref(dense, [F(0)] * len(dense))
    free = [c for c in range(len(cols)) if c not in piv]
    basis = []
    for fc in free:
        v = [F(0)] * len(cols)
        v[fc] = F(1)
        for r, c in enumerate(piv):
            v[c] = -M[r][fc]
        basis.append(v)
    live = [i for i, nm in enumerate(cols) if nm[0] == "t" and any(x[i] != 0 for x in basis)]
    print(f"    {len(eqs)} equations, {len(cols)} unknowns, rank {len(piv)}, "
          f"null dim {len(free)}, spare {len(eqs) - len(piv)}")
    print(f"    bases with a nonzero scale: {len(live)}/{sum(1 for c in cols if c[0]=='t')}"
          + ("" if live else "   <== REFUTED"))
    if excluded:
        print(f"    excluded {len(excluded)} row(s) with != 1 ramified class: {excluded[:6]}")
    # compare any determined ratios against KY (2.17)
    if len(free) == 1:
        v = basis[0]
        print("    unique up to scale; comparing the A-entries against KY (2.17):")
        ref = None
        for i, nm in enumerate(cols):
            if nm[0] == "t":
                continue
            p, k, sgn = nm
            ky = b_ram(p, k, sgn)
            got = v[i]
            if ref is None and ky != 0 and got != 0:
                ref = got / ky
            tag = "match" if (ref is not None and ky * ref == got) else "DIFFERS"
            print(f"       p={p} ord_p(c)={k} chi={sgn:>2}:  fitted {str(got):>10}   "
                  f"KY(2.17) {str(ky):>8}   {tag}")


if __name__ == "__main__":
    main()
