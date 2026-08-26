# SCORECARD: the N=1 deconfounder -- measured 2026-08-26/27, same day

Preregistration: preregN1.md, committed f2d599e BEFORE any measurement.
Bases: 10_1, 14_1, 22_1, 26_1 (D = 2q, q = 5,7,11,13 PRIME, N = 1,
M = 20/28/44/52).  Residuals 2.9e-44, 2.6e-44, 1.7e-43, 7.8e-44.

## THE VERDICT: Reading A, unanimously -- the scale 2 tracks N=1

    E_eis = -2 GROSS(q,1) + 3 GROSS(q,2)   mod cusp,   q = 5, 7, 11, 13

**Identical to the composite bases 210_1/330_1, at PRIME q.**  So the
doubling is NOT caused by q being composite; it tracks N = 1 / the
absent third term.  Reading B (composite q drives the scale), which
predicted (-1,+3/2), is REFUTED at all four bases.

This retires the "composite-D scale 2" framing of the 210_1/330_1
scorecards: those two bases were confounded, and the composite half of
the description was the wrong half.  What survives is:

    N > 1 :  w = (-1/(s-1), (s+1)/(2(s-1)), -1/2),  weight sum 0
    N = 1 :  w = 2 x (first two terms only),        weight sum 1

with the third Gross slot (DN,1) absent at N=1 because it is indefinite
(formal mass 0).  Weight sum = a_0(E_eis) since genus averages are
theta-normalized to constant term 1.

## Scorecard against the preregistration

| # | prediction | outcome |
|---|---|---|
| P1 | support {(q,1),(q,2)}, no third term | **PARTIAL** -- that support does carry E_eis with no third term at all four, but it is NOT unique (see below) |
| P2 | discriminator: (-2,+3) [A] vs (-1,+3/2) [B]; no preference recorded | **READING A**, 4/4 unanimous, and robust to the support ambiguity |
| P3 | mass ratio 3/2 (masses pre-computed) | **HIT** 4/4 -- and it turns out to hold on the competing support too |
| P4 | prime-base constants (Eis dim 11, rank 7) | **REFUTED as stated** -- these are Eis dim 5, kernelrat rank 3, aliasing 1, uniformly at all four (levels 40-104 are far below the 120-408 family). New constants, not a failure of the fits |

## The support is genuinely degenerate at D = 2q, N = 1 (new)

Unlike 210_1/330_1 -- where the competing support was EXCLUDED by rank --
here the competitor {(2,1),(2,q)} also carries E_eis, with

    w = ( -2/(q-1), (q+1)/(q-1) ),  aliasing dim 0, weight sum 1

which is exactly 2x the s-law for the OTHER reading of the same base
(ramified prime 2, s = D/2 = q).  So at D = 2q both ramified primes give
a valid decomposition, each at scale 2.  The verdict above is therefore
support-INDEPENDENT: whichever prime you call q, the weights are twice
that support's s-law and sum to 1, never the unscaled version.  This is
what makes Reading A robust rather than an artifact of the convention.

## THE MASS LAW IS GAUGE-INVARIANT (new, and the strongest structural find)

|w2/w1| = mass2/mass1 holds on BOTH supports at all four bases:

| q | canonical {(q,1),(q,2)} | | competitor {(2,1),(2,q)} | |
|---|---|---|---|---|
|   | \|w2/w1\| | mass ratio | \|w2/w1\| | mass ratio |
| 5  | 3/2 | 3/2 | 3 | 3 |
| 7  | 3/2 | 3/2 | 4 | 4 |
| 11 | 3/2 | 3/2 | 6 | 6 |
| 13 | 3/2 | 3/2 | 7 | 7 |

The mass-ratio law does not depend on the "q = largest ramified prime"
convention at all -- it holds for either decomposition of the same
E_eis.  That makes the mass law the gauge-invariant statement and the
support choice a gauge, which is a much better footing for the
Siegel-Weil derivation than the convention-dependent weight formulas.
Canonical-support count is now 15-for-15, plus 4 alternate-support
confirmations.

## What remains open

The scale 2 is now correctly LOCALIZED (it is the N=1 / missing-third-
term phenomenon) but still underived: why does dropping an
indefinite-slot term double the surviving pair rather than renormalize
some other way?  The weight-sum jump 0 -> 1 is the sharp form of the
question.  The re-split rule (odd part vs largest definite divisor) is
untouched by these bases and still needs odd 4-prime D = 1155.
