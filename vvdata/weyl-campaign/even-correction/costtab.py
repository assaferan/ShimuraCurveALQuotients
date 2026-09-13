#!/usr/bin/env python3
"""Cross-tabulate a correction discriminant's COST against its LEGALITY.

    cost      = amt * deg Z(disc)      -- the degree the perturbed f gains, hence the number of
                                         RATIONAL CM points RationalConstraintsOnEquations needs
    legality  = condition 2 (in image over Q) and condition 3 (integrally solvable over Z)

Inputs: the degz.m log (DEGZ lines) and a PROBE_INTSWEEP log (INTSWEEP lines).

The question it answers: at this base, are the CHEAP discriminants (deg Z = 1) ever legal?
If legality concentrates on large deg Z, the hatch's pricing is structural rather than a bad
choice of discriminant.
"""
import re, sys
from collections import defaultdict

def main(degz_log, sweep_log, amt):
    degz = {}
    for line in open(degz_log):
        m = re.match(r"DEGZ d -(\d+) degZ (\d+)", line)
        if m: degz[int(m.group(1))] = int(m.group(2))

    tot, ok, ram, notimg = defaultdict(int), defaultdict(int), {}, defaultdict(int)
    for line in open(sweep_log):
        m = re.search(r"INTSWEEP key (-?\d+) disc (\d+) isram (\w+) amt (-?\d+) "
                      r"inimage (\w+) intsol (\w+)", line)
        if not m: continue
        d, isram, a = int(m.group(2)), m.group(3) == "true", int(m.group(4))
        if a != amt: continue
        ram[d] = isram
        tot[d] += 1
        if m.group(6) == "true": ok[d] += 1
        if m.group(5) != "true": notimg[d] += 1

    print(f"amt = {amt}\n")
    print(f"{'disc':>6} {'degZ':>5} {'cost':>5} {'ram?':>5} {'keys':>5} {'intsol':>7} "
          f"{'not-in-img':>11}   legal & cheap?")
    print("-" * 78)
    cheap_legal = []
    for d in sorted(tot, key=lambda x: (degz.get(x, 99), x)):
        z = degz.get(d)
        cost = z * amt if z else None
        flag = ""
        if not ram[d] and z == 1 and ok[d] == tot[d] and tot[d] > 0:
            flag = "  <== CHEAP AND LEGAL"; cheap_legal.append(d)
        print(f"{-d:>6} {str(z):>5} {str(cost):>5} {str(ram[d]):>5} {tot[d]:>5} "
              f"{ok[d]:>7} {notimg[d]:>11}{flag}")

    print("\nlegality vs cost:")
    by_z = defaultdict(lambda: [0, 0])
    for d in tot:
        if ram[d] or degz.get(d) is None: continue
        by_z[degz[d]][0] += 1
        if ok[d] == tot[d] and tot[d] > 0: by_z[degz[d]][1] += 1
    for z in sorted(by_z):
        n, k = by_z[z]
        print(f"  deg Z = {z}: {k} of {n} non-ram discriminants integrally solvable at every key")
    print(f"\ncheap (deg Z = 1) AND legal at every key: "
          f"{cheap_legal if cheap_legal else 'NONE'}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 6)
