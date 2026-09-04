#!/usr/bin/env python3
"""Is rem:gauge's six-index ambiguity resolvable on the 158-monomial space, or is it one of
sec:exact's 50 'unremovable' identities?  Run from vvdata/weyl-campaign/ (reads cusp7_15_2.out).

ANSWER (2026-09-04): THE RANK IS 6 -- BUT THAT DOES NOT ANSWER THE QUESTION IN THE TITLE.
158 x 6 has rank 6 (0 free directions) against rank 4 for the nine-form panel. That is solid
arithmetic. It does NOT show the rem:gauge ambiguity is resolvable, and it does NOT contradict
sec:exact, because the two are about different objects:

  * sec:exact counts the span of this character's FORMS (78 of 128 joint (oo,0) slots, 50
    identities).
  * this script ranks MONOMIALS.

Monomials are not valid probes of c_{eta*}(0) at all: it sums over EVERY cusp class, and unlike
a genuine Borcherds principal part (prop:nohalf / [GY, Lemma 24]) an individual eta-monomial can
carry a nonzero constant term at INTERMEDIATE cusps that an (A_m, B_j)-only model never sees.
PLAN.md records that the numeric solve against real c_{eta*}(0) values is still inconsistent.

So: do NOT cite this as evidence the paper is wrong. It was briefly written up that way on
2026-09-04 and retracted the same day. The open question is the FORM-space rank, uncomputed.

LESSON: the two validations below make the NUMBER trustworthy (they caught Magma line-wrapping
and silently-dropped rational coefficients). They say nothing about whether the right OBJECT is
being measured. A validated computation of the wrong thing is still the wrong thing.

⚠ TWO SILENT PARSE TRAPS -- both give a WRONG matrix that still looks plausible:
  1. cusp7.m had no SetColumns(0), so Magma wrapped FORMC across many physical lines with NO
     continuation character. A line-based parse reads only the first fragment. (Fixed in cusp7.m,
     but old .out files on disk are still wrapped -- parse whole-file, not line-by-line.)
  2. FORMC coefficients are EXACT RATIONALS: 24 of 70 pairs per form look like <36, 134/3>.
     An integer regex drops them silently.
Both traps produced "rank 6" from corrupted data on the first attempt -- the right answer for the
wrong reason. ALWAYS run the two validations below first; they catch both.
"""
import re
from math import comb
from fractions import Fraction as F

DS = [1,2,3,4,5,6,10,12,15,20,30,60]
IDX = [1,2,3,10,15,30]                      # the six indices of rem:gauge
txt = open("cusp7_15_2.out").read()

monos = {int(a): [int(x) for x in re.findall(r'-?\d+', b)]
         for a, b in re.findall(r'MONO\s+(\d+)\s+\[([^\]]*)\]', txt, re.S)}
formc = {}
for k, body in re.findall(r'FORMC\s+(-?\d+)\s+\[(.*?)\]\s*\n(?=FORMC|MONO|PP|TG|SELFC|a0|BASE|#|DONE|$)', txt, re.S):
    formc[int(k)] = {int(a): F(b.replace(' ', ''))
                     for a, b in re.findall(r'<\s*(\d+),\s*([-\d/ ]+?)\s*>', body)}

def eta(r, top):
    """oo-side q-expansion of prod_d eta(d tau)^{r_d}; order at oo is (sum d r_d)/24 (Ligozat)."""
    L24 = sum(d*x for d, x in zip(DS, r)); assert L24 % 24 == 0
    L = L24 // 24; n_ = top - L + 1
    if n_ <= 0: return {}
    ser = [0]*(n_+1); ser[0] = 1
    for d, rd in zip(DS, r):
        if rd == 0: continue
        n = 1
        while d*n <= n_:
            e, k = d*n, abs(rd)
            co = ([(-1)**j*comb(k, j) for j in range(min(k, n_//e)+1)] if rd > 0
                  else [comb(k-1+j, j) for j in range(n_//e+1)])
            new = [0]*(n_+1)
            for i, a in enumerate(ser):
                if a:
                    for j, c in enumerate(co):
                        p = i + j*e
                        if p > n_: break
                        new[p] += a*c
            ser = new; n += 1
    return {L+i: c for i, c in enumerate(ser) if c}

def rank(rows, nc):
    M = [list(map(F, r)) for r in rows]; r = 0
    for c in range(nc):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]/M[r][c]; M[i] = [x-f*y for x, y in zip(M[i], M[r])]
        r += 1
    return r

pp = {ri: eta(v, 40) for ri, v in monos.items()}

# VALIDATION 1: form -1's principal part must be {-10:2, -2:-2, -1:2} (tests/M0Multiplier.m)
tot = {}
for ri, c in formc[-1].items():
    for e, v in pp[ri].items():
        if e < 0: tot[e] = tot.get(e, 0) + c*v
tot = {e: v for e, v in sorted(tot.items()) if v}
assert tot == {-10: 2, -2: -2, -1: 2}, f"VALIDATION 1 FAILED: {tot}"
print("validation 1 (form -1 principal part)      : OK", tot)

# VALIDATION 2: nine-form panel must have rank 4 (the paper's own figure, rem:gauge)
panel = [[sum(c*pp[ri].get(-m, 0) for ri, c in formc[k].items()) for m in IDX] for k in sorted(formc)]
rp = rank(panel, len(IDX))
assert rp == 4, f"VALIDATION 2 FAILED: panel rank {rp}, paper says 4"
print("validation 2 (nine-form panel rank == 4)   : OK")

ra = rank([[pp[ri].get(-m, 0) for m in IDX] for ri in sorted(monos)], len(IDX))
print(f"\nRESULT: 158-monomial matrix 158 x {len(IDX)}, rank {ra}, free directions {len(IDX)-ra}")
print("  rank 6 => the six-index ambiguity is RESOLVED by probing beyond the panel.")
