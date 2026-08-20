#!/usr/bin/env python3
"""
Test the SUPPORT RULE and extract the VALUES of the m=0 weights u_m, on any base.

RULE (verified on 6_5, 15_2, 10_11): u_m != 0 only when the ternary quadratic SPACE fails to
represent m at EXACTLY ONE finite place.  Failing places lie over p | D.

Inputs per base:
  * principal parts   -- from an allterms_<D>_<N>.txt dump, or from the [PP] lines of an oracle log
  * Diff sets         -- from diffset.m output  (diff_<D>_<N>.txt)
  * measured deltas   -- the ground-truth multipliers (aggregate_oracle.py / fit_signs.py)

Output: whether the constrained system is consistent, how many SPARE CONDITIONS it satisfies (that is
the real evidence), and the solved weights (pinned ones flagged), tagged with the failing prime.

usage: support_fit.py <base>            e.g.  support_fit.py 6_5
       support_fit.py --all
"""
import sys, os, re, glob
from fractions import Fraction as F

WT = "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/fe0e5880-9d7a-4747-9e91-68091e8d7211/scratchpad/m0-multiterm"

# ground-truth multipliers (units of log N).  15_2/10_11 confirmed end-to-end; the rest from the sweep.
MEASURED = {
    "15_2":  {-2: 2, -1: 4, 9: 0, 10: 0, 11: 4, 12: 2, 13: 4, 14: -2, 15: 2},
    "10_11": {-2: 0, -1: 0, 9: 1, 10: 1, 11: 2, 12: 2, 13: 1, 14: 0, 15: 1},
    "6_5":   {-2: 0, -1: 0, 9: 6, 10: 3, 13: 3, 14: 0, 15: 3},
    "10_3":  {9: 6, 10: 3, 11: 3, 12: 3, 13: 3},
    "6_11":  {9: 3, 11: 6, 13: 3, 14: 0, 15: 3},
    "6_17":  {10: 3, 11: 6, 13: 3, 14: 0, 15: 3},
    "6_19":  {9: 0, 11: 1, 13: 1, 15: 2},
    "6_23":  {10: 2, 13: 2, 14: 0, 15: 3},
    "6_31":  {13: 0},
    "6_41":  {13: 3, 15: 4},
    "10_13": {9: 0, 10: 1, 13: 1, 14: 2},
    "10_17": {13: 1, 14: 0},
    "10_19": {9: 0, 11: 1, 13: 1, 14: 0, 15: 0},
    "14_11": {9: 1, 10: 0, 11: 1, 13: 1, 14: 1, 15: 2},
    "22_7":  {10: 0, 12: 1, 13: 1, 14: 1, 15: 1},
    "26_5":  {9: 1, 12: 0, 13: 1, 14: 1, 15: 0},
    "34_5":  {9: 2, 10: 0, 12: 1, 13: 1, 14: 1, 15: 1},
    "34_7":  {9: 0, 12: 1, 13: 1, 14: 0},
    "6_29":  {14: 0},
}
# hauptmodul multipliers measured by the cheap oracle (harvest_haupt); only where both are known
HAUPT = {"6_7": {-1: 1, -2: 2}, "10_11": {-1: 0, -2: 0}, "6_5": {-1: 0, -2: 0},
         "10_3": {-1: 0, -2: 0}, "6_11": {-1: 0, -2: 0}, "6_13": {-1: 3, -2: 0},
         "6_17": {-1: 0, -2: 0}, "6_19": {-1: 0, -2: 0}, "6_23": {-1: 0, -2: 0}}


def sqfree(n):
    s, d = 1, 2
    n = abs(n)
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            s *= d
        d += 1
    return s * n


def parse_allterms(path):
    """form -> {index: coeff}; index = ('oo', m) or ('0', j).  Also returns M and per-index local data."""
    forms, loc, M = {}, {}, None
    cur = None
    for line in open(path):
        m = re.search(r"M = (\d+)", line)
        if m and M is None:
            M = int(m.group(1))
        m = re.match(r"--------- form\[(-?\d+)\] ---------", line.strip())
        if m:
            cur = int(m.group(1))
            forms[cur] = {}
            continue
        if cur is None:
            continue
        m = re.match(r"\s*oo m=(\d+)\s+c=(-?\d+).*locs\(p,W,v_p\(r\)\)=\[(.*)\]", line)
        if m:
            idx = ("oo", int(m.group(1)))
            forms[cur][idx] = F(int(m.group(2)))
            loc[idx] = [(int(p), F(w.strip()), int(v))
                        for p, w, v in re.findall(r"<\s*(\d+),\s*([^,]+),\s*(-?\d+)\s*>", m.group(3))]
            continue
        m = re.match(r"\s*0\s+j=(\d+)\s+c=(-?\d+)\s+r=(\S+)\s", line)
        if m:
            idx = ("0", int(m.group(1)))
            forms[cur][idx] = F(int(m.group(2)))
            continue
    return forms, loc, M


def parse_pp(path):
    """Parse the (line-wrapped) [PP] dump of an oracle log."""
    txt = open(path).read()
    # rejoin magma's 80-col wrapping: a [PP] record runs until the next line starting a new record
    recs = re.findall(r"\[PP\] (\d+) (\d+) form=(-?\d+) oo=\[(.*?)\] zero=\[(.*?)\]", txt, re.S)
    forms = {}
    for D, N, k, oo, zero in recs:
        f = {}
        for a, b in re.findall(r"<\s*(\d+),\s*(-?\d+)\s*>", oo):
            f[("oo", int(a))] = F(int(b))
        for a, b in re.findall(r"<\s*(\d+),\s*(-?\d+)\s*>", zero):
            f[("0", int(a))] = F(int(b))
        forms[int(k)] = f
    return forms


def parse_diff(path, sgn=1):
    """m -> list of finite places where the SPACE fails to represent m."""
    out = {}
    for line in open(path):
        m = re.match(r"\[DIFF\] \d+ \d+ sgn=(-?\d+)\s+m=(\d+)\s+finite_fail=\[(.*?)\]", line)
        if m and int(m.group(1)) == sgn:
            out[int(m.group(2))] = [int(x) for x in re.findall(r"\d+", m.group(3))]
    return out


def rref(A, b):
    m, n = len(A), len(A[0]) if A else 0
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                fq = M[i][c]
                M[i] = [x - fq * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if all(M[i][c] == 0 for c in range(n)) and M[i][n] != 0:
            return piv, M, False
    return piv, M, True


def index_m(idx, Mlev):
    """the integer whose square class the Diff test should use for this principal-part index"""
    if idx[0] == "oo":
        return idx[1]
    r = F(idx[1], Mlev)          # cusp-0 index j <-> r = j/M
    return sqfree(r.numerator * r.denominator)


def run(base, verbose=True):
    D, N = (int(x) for x in base.split("_"))
    meas = dict(MEASURED.get(base, {}))
    meas.update({k: v for k, v in HAUPT.get(base, {}).items() if k not in meas})
    if not meas:
        return None
    at = os.path.join(WT, f"allterms_{base}.txt")
    if os.path.exists(at):
        forms, loc, Mlev = parse_allterms(at)
    else:
        forms, loc, Mlev = parse_pp(os.path.join(WT, "oracle_out", f"{base}.log")), {}, None
    dpath = os.path.join(WT, f"diff_{base}.txt")
    if not os.path.exists(dpath):
        print(f"{base}: no diff file ({dpath}) -- run diffset.m")
        return None
    diff = parse_diff(dpath)
    idxs = sorted({i for f in forms.values() for i in f}, key=lambda t: (t[0], t[1]))
    if Mlev is None and any(i[0] == "0" for i in idxs):
        print(f"{base}: cusp-0 indices present but M unknown (need an allterms dump)")
        return None
    free, forced = [], []
    for i in idxs:
        mm = index_m(i, Mlev)
        fset = diff.get(mm, [])
        (free if len(fset) == 1 else forced).append((i, fset))
    rows, tg, used = [], [], []
    for k, f in forms.items():
        if k not in meas:
            continue
        rows.append([f.get(i, F(0)) for i, _ in free])
        tg.append(F(meas[k]))
        used.append(k)
    if not rows:
        return None
    piv, M, ok = rref(rows, tg)
    spare = len(rows) - len(piv) if ok else None
    if verbose:
        print(f"===== {base}  (D={D}, N={N}) =====")
        print(f"  indices: {len(idxs)}   free (|F|=1): {len(free)}   forced 0: {len(forced)}")
        for i, fs in free:
            print(f"     FREE   {str(i):>10}  m={index_m(i, Mlev):<4} fails at {fs}")
        for i, fs in forced:
            print(f"     zero   {str(i):>10}  m={index_m(i, Mlev):<4} fails at {fs}")
        print(f"  {len(rows)} measured forms {used}")
        if not ok:
            print("  ==> INCONSISTENT -- support rule REFUTED on this base\n")
            return (base, False, None, None)
        nfree = len(free) - len(piv)
        print(f"  ==> CONSISTENT: rank {len(piv)}, {nfree} residual free param(s), "
              f"{spare} SPARE CONDITION(S) satisfied")
        sol = [F(0)] * len(free)
        for r, c in enumerate(piv):
            sol[c] = M[r][len(free)]
        basis = []
        for fc in [c for c in range(len(free)) if c not in piv]:
            v = [F(0)] * len(free)
            v[fc] = F(1)
            for r, c in enumerate(piv):
                v[c] = -M[r][fc]
            basis.append(v)
        for k, (i, fs) in enumerate(free):
            drift = [b[k] for b in basis]
            tag = "PINNED" if all(d == 0 for d in drift) else "free"
            print(f"     u{str(i):<12} = {str(sol[k]):>8}   [{tag}]  fail p={fs[0] if fs else '-'}")
        print()
        return (base, True, spare, {free[k][0]: sol[k] for k in range(len(free))})


if __name__ == "__main__":
    args = sys.argv[1:]
    bases = sorted(set(MEASURED) | set(HAUPT)) if (not args or args[0] == "--all") else args
    res = []
    for b in bases:
        r = run(b)
        if r:
            res.append(r)
    print("summary:")
    for b, ok, spare, _ in res:
        print(f"   {b:>7}  {'CONSISTENT' if ok else 'REFUTED':<11} spare={spare}")
