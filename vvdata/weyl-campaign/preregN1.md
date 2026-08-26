# PREREGISTRATION: the N=1 deconfounder (written BEFORE any measurement)

Date: 2026-08-26/27.  Bases: X0^10(1), X0^14(1), X0^22(1), X0^26(1) --
D = 2q with q PRIME (5, 7, 11, 13), N = 1, M = 20/28/44/52.

## The confound this resolves

Across the eleven banked bases, "composite q" and "N = 1" coincide
EXACTLY: the nine prime-q bases all have N > 1 (weights summing to 0,
third term -1/2 present), and the two composite-q bases 210_1 and 330_1
both have N = 1 (weights (-2,+3), summing to +1, no third term).  So
nothing measured so far can attribute the doubling to either cause.

Note the algebraic coincidence that produced the confound: q composite
forces omega(D) = 4 hence D >= 210, and those bases were run at N = 1;
while every base cheap enough to run before the closed form had q prime
and N > 1.  These four bases are the missing cell -- prime q AND N = 1 --
and they only became reachable when N=1 support (tracking the zero
coset, #iso = 2N-1 = 1) entered the pipeline for 210_1.

## The two readings, which separate here

Let s = D/q = 2 at all four bases.  The prime s=2 law is
w = (-1/(s-1), (s+1)/(2(s-1)), -1/2) = (-1, 3/2, -1/2).

* **Reading A -- the doubling tracks N=1 / the absent third term.**
  At N = 1 the third slot (DN,1) = (2q,1) has omega = 2, hence is
  INDEFINITE with formal mass 0 and cannot appear.  Prediction:
  **w = (-2, +3) on {(q,1),(q,2)}**, weight sum +1, exactly as at both
  composite bases.
* **Reading B -- the doubling tracks composite q.**  With q prime the
  scale is 1.  Prediction: **w = (-1, +3/2) on {(q,1),(q,2)}**, weight
  sum +1/2.

The readings are distinguishable by a factor of 2 on rank-certified
weights, and the weight sum (= a_0(E_eis), since genus averages are
theta-normalized to constant term 1) distinguishes them independently.

## Preregistered predictions

P1 (support): E_eis in span{ GROSS(q,1), GROSS(q,2) } mod cusp, with NO
third term, at all four bases.  Falsifier: any fit needing a third
Gross column, or support outside {(q,*)}.

P2 (the discriminator): w = (-2,+3)  [Reading A]  vs  (-1,+3/2)
[Reading B].  I record NO preference between them -- this is the point
of the experiment.  Whichever fires must fire at all four bases; a split
verdict across q = 5,7,11,13 would refute BOTH readings and mean the
scale depends on something else again.

P3 (mass-ratio law, bases 12-15): |w2/w1| = mass(q,2)/mass(q,1) = 3/2
under BOTH readings (the ratio is scale-invariant), extending the law
that has now held 11-for-11.  Masses COMPUTED IN ADVANCE (massN1.m/.log,
committed with this file before any fit was run):

| q | mass(q,1) | mass(q,2) | ratio |
|---|---|---|---|
| 5  | 1/12 | 1/8  | 3/2 |
| 7  | 1/8  | 3/16 | 3/2 |
| 11 | 5/24 | 5/16 | 3/2 |
| 13 | 1/4  | 3/8  | 3/2 |

All four ratios are 3/2 as required, so P3 is a genuine test of the fits
rather than of the masses.

P4 (structural): these are prime-q bases, so the prime-base constants
are expected -- Eisenstein dim 11 and constants-map rank 7 -- NOT the
composite values (rank 15, aliasing 25) seen at 210_1/330_1.  Soft
prediction; deviation is informative, not fatal.

## Protocol

eis32k.m (closed-form rhov, validation-gated at 210_1 in 46d6849), eta
pool from enum32m via a synthetic mono dump carrying the universal
character key (s1 = s2 = 0 mod 24, parity bit 1 at p=2 only), then
kernelrat -> eisrat -> genallgross -> allgross over all definite (D',R)
with D'R | D.  Gauge checks via genrestrict.py: the winning support must
carry E_eis alone, and competing supports must fail rank.  SOLVE resid
and "residual beyond rank" are the arbiters; eisrat's RATIONAL banner
alone is NOT trusted.

All four bases are small (36-84 cosets), so all four run to completion
before any verdict is read.
