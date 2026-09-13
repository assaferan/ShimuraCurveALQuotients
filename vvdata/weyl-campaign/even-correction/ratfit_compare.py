#!/usr/bin/env python3
"""Compare the RATIONAL linear fit (RationalConstraintsOnEquations) between a baseline
run and an even-perturbed run of the same base.

WHY THIS EXISTS.  The error raised at EquationsCovers.m:68,
    "Error in Schofer table values at rational points - no solution found!"
is raised INSIDE QuadraticConstraintsOnEquations but tests `kernels[j]`, which is built by
RationalConstraintsOnEquations.  So the failing object is the kernel of the rational linear
system y2 = f(s), deg f <= 2g+2 -- not anything the quadratic stage computes.  An empty kernel
means the rational CM values fit NO such polynomial.

Two readings, and this script separates them:
  (a) DEGREE.   The perturbed y2 is the old one times the square of a nonconstant function,
                so f_new has degree 2g+2+2d and the ansatz is simply too short.  Signature:
                the kernel reappears at a larger degree bound, and y2_pert/y2_base is a
                perfect square at every rational CM point.
  (b) INCONSISTENCY.  The values do not come from a curve of this shape at all.  Signature:
                M stays full rank at every degree bound tried.
"""
import re, sys
from fractions import Fraction

def parse(path):
    hdr, sweep, data = {}, {}, {}
    for line in open(path):
        m = re.match(r"RATFIT key (-?\d+) W=(\[[^\]]*\]) g=(\d+) #ds=(\d+) #svals=(\d+) "
                     r"#rat=(\d+) #quad=(\d+) dimB=(\d+)", line)
        if m:
            k = int(m.group(1))
            hdr[k] = dict(W=m.group(2), g=int(m.group(3)), nds=int(m.group(4)),
                          nsvals=int(m.group(5)), nrat=int(m.group(6)),
                          nquad=int(m.group(7)), dimB=int(m.group(8)))
            continue
        m = re.match(r"RATFIT   key (-?\d+) degbound=(\d+) nrows=(\d+) ncols=(\d+) "
                     r"rank=(\d+) dimker=(\d+)", line)
        if m:
            k = int(m.group(1))
            sweep.setdefault(k, []).append(tuple(int(m.group(i)) for i in range(2, 7)))
            continue
        m = re.match(r"RATFITDATA key (-?\d+) (s|y2)=\[(.*)\]\s*$", line)
        if m:
            k, which = int(m.group(1)), m.group(2)
            vals = [Fraction(t.strip().replace(" ", "")) for t in m.group(3).split(",")
                    if t.strip()]
            data.setdefault(k, {})[which] = vals
    return hdr, sweep, data

def is_square(q):
    if q < 0:
        return False
    n, d = q.numerator, q.denominator
    return round(n ** .5) ** 2 == n and round(d ** .5) ** 2 == d

def main(base_log, pert_log):
    bh, bs, bd = parse(base_log)
    ph, ps, pd = parse(pert_log)
    print(f"keys: baseline {sorted(bh)}\n      perturbed {sorted(ph)}\n")
    for k in sorted(set(bh) & set(ph)):
        b, p = bh[k], ph[k]
        print(f"--- key {k}  W={b['W']}  g={b['g']}")
        print(f"    #rat CM pts   base {b['nrat']:3d}   pert {p['nrat']:3d}"
              f"      #quad  base {b['nquad']:3d}  pert {p['nquad']:3d}")
        print(f"    dim kernel    base {b['dimB']:3d}   pert {p['dimB']:3d}"
              f"       <-- 0 on the perturbed side is the failure")
        print("    degree sweep (degbound, rank, dimker):")
        for tag, sw in (("base", bs.get(k, [])), ("pert", ps.get(k, []))):
            print(f"      {tag}: " + "  ".join(f"{d}:r{r}/k{kk}" for d, _, _, r, kk in sw))
        # the ratio test: same s-values?
        if k in bd and k in pd and "s" in bd[k] and "s" in pd[k]:
            sb, sp = bd[k]["s"], pd[k]["s"]
            if sb != sp:
                print(f"    ⚠ rational CM s-values DIFFER between runs "
                      f"({len(sb)} vs {len(sp)}) -- ratio test not applicable")
                continue
            yb, yp = bd[k]["y2"], pd[k]["y2"]
            ratios = [(s, (y2 / y1) if y1 != 0 else None)
                      for s, y1, y2 in zip(sb, yb, yp)]
            nsq = sum(1 for _, r in ratios if r is not None and is_square(r))
            ncon = len({r for _, r in ratios if r is not None})
            print(f"    y2 ratio pert/base: {nsq}/{len(ratios)} are perfect squares in Q; "
                  f"{ncon} distinct value(s)")
            for s, r in ratios:
                print(f"        s={s}  ratio={r}  square={is_square(r) if r is not None else '-'}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
