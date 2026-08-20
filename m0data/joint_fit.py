#!/usr/bin/env python3
"""
Fit the m=0 weights across MANY bases at once, by exploiting that the local factor at a ramified
prime q | D depends only on the local lattice at q -- hence is SHARED by every base with that q,
independently of N.  (The per-base caveat in memory applies to the N-part, not to the q-part.)

MODEL.  The code computes each term as
      contrib(m) = pre * c(m) * cond_half * (h/w) * (en/ed) * prod_p W_p(m)
and is believed correct except at the ramified (anisotropic) primes q | D, where the Eisenstein
coefficient is NOT the lattice density.  So write

      u_m  =  C_base * P(m) * prod_{q | D} A_q(cls_q(m))

      P(m)  = cond_half * (h/w) * (en/ed) * prod_{q not | D} W_q(m)      [taken from the dump]
      A_q   = UNKNOWN table indexed by the local square class of m at q
              cls_q(m) = (ord_q(m), unit class)   -- unit class = Legendre symbol (q odd), or mod 8 (q=2)

Each measured multiplier gives  M_k = sum_m c_{k,m} u_m, i.e.

      M_k * t_base  -  sum_m c_{k,m} P(m) prod_q A_q(cls_q(m))  =  0,      t_base := 1/C_base

which is LINEAR AND HOMOGENEOUS in the unknown vector (A-entries..., t_base...) as long as at most one
ramified prime has an unknown per term -- for |D| = 2 ramified primes the product A_2*A_3 is itself the
unknown, indexed by the PAIR of classes.  We therefore solve for the pair-indexed table T[c2,c3] (that
is exactly what the data can see) and only afterwards ask whether it factors as A_2 * A_3.

The null space should be 1-dimensional (overall scale).  Extra dimensions = genuine degeneracies, and
they are reported, never silently resolved.

usage: joint_fit.py [base ...]
"""
import sys, os, re, glob
from fractions import Fraction as F
from support_fit import (WT, MEASURED, HAUPT, parse_allterms, parse_diff, sqfree, index_m, rref)


ANSATZ = "full"
UNKNOWN_LEVEL = False
LEVEL_SUPPORT = False


def cls(q, m):
    """local square class of the integer m at q, as a hashable label.

    The ansatz controls how coarse the indexing is -- coarser = fewer unknowns = a real test
    (it can be REFUTED by inconsistency) rather than a fit."""
    v = 0
    while m % q == 0:
        m //= q
        v += 1
    unit = (m % 8) if q == 2 else (1 if pow(m % q, (q - 1) // 2, q) == 1 else -1)
    if ANSATZ == "full":
        return (q, v, unit)
    if ANSATZ == "val":                 # depends on the valuation only
        return (q, v)
    if ANSATZ == "val+unit2":           # unit class matters only at q = 2
        return (q, v, unit) if q == 2 else (q, v)
    if ANSATZ == "parity":              # depends on v mod 2 and the unit class
        return (q, v % 2, unit)
    if ANSATZ == "stab":                # valuation stabilises at 2
        return (q, min(v, 2), unit)
    raise SystemExit(f"unknown ansatz {ANSATZ}")


def parse_terms(path):
    """form -> list of (index, c, P, m_int); P = the code's non-ramified prefactor product.
       Needs D to know which primes are ramified."""
    D = N = Mlev = None
    forms, cur = {}, None
    for line in open(path):
        mm = re.search(r"X0\^(\d+)\((\d+)\)", line)
        if mm:
            D, N = int(mm.group(1)), int(mm.group(2))
        mm = re.search(r"M = (\d+)", line)
        if mm and Mlev is None:
            Mlev = int(mm.group(1))
        mm = re.match(r"--------- form\[(-?\d+)\] ---------", line.strip())
        if mm:
            cur = int(mm.group(1))
            forms[cur] = []
            continue
        if cur is None:
            continue
        mm = re.match(r"\s*(oo m|0  j)=(\d+)\s+c=(-?\d+)\s.*", line)
        if not mm:
            continue
        block = "oo" if mm.group(1).startswith("oo") else "0"
        idx = (block, int(mm.group(2)))
        c = F(int(mm.group(3)))
        det = re.search(r"cond/2=(\S+) en/ed=(\S+) g=(\S+) locs\(p,W,v_p\(r\)\)=\[(.*)\]", line)
        if det:
            cond = F(det.group(1))
            ened = F(det.group(2))
            locs = [(int(p), F(w.strip())) for p, w, _ in
                    re.findall(r"<\s*(\d+),\s*([^,]+),\s*(-?\d+)\s*>", det.group(4))]
            hw = re.search(r"h=(\d+) w=(\d+)", line)
            hwv = F(int(hw.group(1)), int(hw.group(2)))
            ram = set(_pdivs(D)) | ({N} if UNKNOWN_LEVEL else set())
            P = cond * hwv * ened
            for p, w in locs:
                if p not in ram:
                    P *= w
            forms[cur].append((idx, c, P, None))
        else:
            forms[cur].append((idx, c, None, None))   # cusp-0 header line: local data on eta lines
    return D, N, Mlev, forms


def _pdivs(n):
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


def build(bases):
    """returns (colnames, rows) of the homogeneous system"""
    rows, colnames, colidx = [], [], {}

    def col(name):
        if name not in colidx:
            colidx[name] = len(colnames)
            colnames.append(name)
        return colidx[name]

    info = []
    for base in bases:
        path = None
        for cand in (os.path.join(WT, f"allterms_{base}.txt"),
                     os.path.join(WT, "atdump", f"allterms_{base}.txt")):
            if os.path.exists(cand) and "DONE" in open(cand).read():
                path = cand
                break
        if path is None:
            continue
        D, N, Mlev, forms = parse_terms(path)
        meas = dict(MEASURED.get(base, {}))
        meas.update({k: v for k, v in HAUPT.get(base, {}).items() if k not in meas})
        meas = {k: v for k, v in meas.items() if k in forms}
        # skip terms we cannot model yet (cusp-0 blocks carry an extra eta dependence)
        usable = True
        for k in meas:
            if any(P is None for _, _, P, _ in forms[k]):
                usable = False
        if not meas or not usable:
            info.append((base, len(meas), "skipped (cusp-0 block or no measurement)"))
            continue
        ram = _pdivs(D) + ([N] if UNKNOWN_LEVEL else [])
        tcol = col(f"t[{base}]")
        n_before = len(colnames)
        for k, val in sorted(meas.items()):
            row = {}
            row[tcol] = F(val)
            for idx, c, P, _ in forms[k]:
                m_int = index_m(idx, Mlev)
                if LEVEL_SUPPORT and m_int % N != 0:
                    continue          # the level rule: this term carries no weight
                key = tuple(("lvl",) + cls(q, m_int)[1:] if q == N and UNKNOWN_LEVEL
                        else cls(q, m_int) for q in ram)
                row[col(f"A{key}")] = row.get(col(f"A{key}"), F(0)) - c * P
            rows.append(row)
        info.append((base, len(meas), f"{len(colnames) - n_before} new class column(s)"))
    dense = [[r.get(j, F(0)) for j in range(len(colnames))] for r in rows]
    return colnames, dense, info


def nullspace(rows, n):
    piv, M, _ = rref(rows, [F(0)] * len(rows))
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = [F(0)] * n
        v[fc] = F(1)
        for r, c in enumerate(piv):
            v[c] = -M[r][fc]
        basis.append(v)
    return piv, free, basis


def main():
    global ANSATZ
    args = sys.argv[1:]
    global UNKNOWN_LEVEL, LEVEL_SUPPORT
    while args and args[0].startswith("--"):
        a = args.pop(0)
        if a.startswith("--ansatz="):
            ANSATZ = a.split("=", 1)[1]
        elif a == "--unknown-level":
            UNKNOWN_LEVEL = True
        elif a == "--level-support":
            LEVEL_SUPPORT = True
    bases = args or sorted(set(MEASURED) | set(HAUPT))
    print(f"[ansatz={ANSATZ}, unknown-level={UNKNOWN_LEVEL}, level-support={LEVEL_SUPPORT}]")
    names, rows, info = build(bases)
    print("bases used:")
    for b, nm, note in info:
        print(f"   {b:>7}  {nm} equation(s)   {note}")
    if not rows:
        print("no equations")
        return
    print(f"\n{len(rows)} equations, {len(names)} unknowns")
    piv, free, basis = nullspace(rows, len(names))
    print(f"rank {len(piv)}, null space dimension {len(free)}")
    print(f"=> {len(rows) - len(piv)} spare condition(s) satisfied "
          f"(equations beyond what the unknowns can absorb)")
    # A null space in which every vector has t = 0 means NO model with a finite per-base scale
    # reproduces the measurements: the model is REFUTED, not merely degenerate.
    tcols = [i for i, nm in enumerate(names) if nm.startswith("t[")]
    live = [i for i in tcols if any(v[i] != 0 for v in basis)]
    print(f"bases with a nonzero scale available: {len(live)}/{len(tcols)}"
          + ("" if live else "   <== REFUTED (every solution forces the scale to 0)"))
    if not live:
        return
    # normalise: pin the first base's scalar to 1 and report what is then determined
    ref = live[0]
    e = [F(0)] * len(names)
    e[ref] = F(1)
    piv2, M2, ok = rref(rows + [e], [F(0)] * len(rows) + [F(1)])
    if not ok:
        print("\ninconsistent after normalisation (should not happen)")
        return
    part = [F(0)] * len(names)
    for r, c in enumerate(piv2):
        part[c] = M2[r][len(names)]
    free2 = [c for c in range(len(names)) if c not in piv2]
    b2 = []
    for fc in free2:
        v = [F(0)] * len(names)
        v[fc] = F(1)
        for r, c in enumerate(piv2):
            v[c] = -M2[r][fc]
        b2.append(v)
    print(f"\nnormalised at {names[ref]} = 1  ({len(free2)} residual free parameter(s)):")
    for i, nm in enumerate(names):
        drift = [v[i] for v in b2]
        if all(d == 0 for d in drift):
            print(f"   {nm:<26} = {str(part[i]):>10}   PINNED")
    print("   -- undetermined:", ", ".join(nm for i, nm in enumerate(names)
                                           if any(v[i] != 0 for v in b2)) or "(none)")


if __name__ == "__main__":
    main()
