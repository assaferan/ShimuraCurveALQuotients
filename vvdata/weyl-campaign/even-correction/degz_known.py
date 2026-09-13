#!/usr/bin/env python3
"""DERIVE the deg Z(d) values that degz.m should validate against, at ANY base.

degz.m's self-validation is currently hardcoded to 34_3's five values.  This derives the
equivalent list at any base from two things the pipeline already emits:

  * the per-key Borcherds divisor  (PROBEEVEN "div_f" lines, minus the injected perturbation)
  * the degree of each key's committed f  (data/models/models_<D>_<N>.m)

using the two identities established at 34_3 (QUADCONSTRAINTS.md sec 7b):

  (a) f is a polynomial, so its only pole is s = infinity; the divisor's NEGATIVE entry is that
      pole, and     deg f = |mult at the pole disc| * deg Z(pole disc)
  (b) a function's divisor has degree 0:   sum over all entries of  mult * deg Z(d)  =  0

Together these are an overdetermined integer linear system in the unknowns deg Z(d).  At 34_3 it
is consistent and pins deg Z(3)=deg Z(24)=deg Z(51)=deg Z(408)=1, deg Z(68)=2, and reproduces all
7 polynomial degrees -- which is what makes degz.m's sweep trustworthy there.

Usage:  degz_known.py <probeeven.log> <models_D_N.m> [perturbed_disc]
Output: a KNOWN string to hand to degz.m, plus the consistency report.

STATUS: DRAFT -- not yet applied to the repo (jobs were running from the tree).  Wants a negative
control before it is trusted: perturb one input degree and check the system reports inconsistency.
"""
import re, sys
from fractions import Fraction

def parse_divs(path, pert_disc):
    divs = {}
    for line in open(path):
        m = re.match(r"PROBEEVEN key (-?\d+) div_f \[(.*)\]\s*$", line)
        if not m: continue
        k = int(m.group(1))
        ent = [(int(a), int(b)) for a, b in
               re.findall(r"<\s*(-?\d+),\s*(-?\d+)\s*>", m.group(2))]
        if pert_disc is not None:
            ent = [(d, c) for d, c in ent if d != pert_disc]
        divs[k] = ent
    return divs

def parse_degrees(path):
    degs = {}
    for line in open(path):
        m = re.match(r"models\[\[([^\]]*)\]\]\s*:=\s*\[\*\s*<(\d+),\s*P!\[([^\]]*)\]", line)
        if not m: continue
        W = tuple(int(x) for x in re.findall(r"-?\d+", m.group(1)))
        coeffs = [c.strip() for c in m.group(3).split(",") if c.strip()]
        deg = len(coeffs) - 1
        while deg >= 0 and coeffs[deg] in ("0",): deg -= 1
        degs[W] = deg
    return degs

def main(divlog, modelfile, pert=None):
    divs = parse_divs(divlog, pert)
    degs = parse_degrees(modelfile)
    print(f"parsed {len(divs)} divisors, {len(degs)} committed model degrees\n")
    print("per-key: pole disc, |pole mult|, and the implied deg Z(pole)")
    for k, ent in sorted(divs.items()):
        neg = [(d, c) for d, c in ent if c < 0]
        if len(neg) != 1:
            print(f"  key {k}: {len(neg)} negative entries -- identity (a) does not apply"); continue
        pd, pc = neg[0]
        print(f"  key {k}: pole at {pd}, mult {pc}   =>  deg f = {abs(pc)} * deg Z({pd})")
    print("\ndegree-0 relations, one per key:")
    for k, ent in sorted(divs.items()):
        print("  key %-6d " % k + "  +  ".join(f"{c}*z[{d}]" for d, c in ent) + "  =  0")
    print("\n(solve these together with deg f = |pole mult| * z[pole]; the system is")
    print(" overdetermined, and its consistency IS the validation.)")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else None)
