# SCORECARD: the eleventh base X0^330(1) -- measured 2026-08-26, same day

Preregistration: prereg330_1.md, committed cb634d7 BEFORE the measurement.
Measurement: eis32k_330_1.out.gz (SOLVE resid 1.709e-42, |rho| 1.618),
allgross_330_1.log, restrict330A/B.log.  First base ever run at M = 660 --
only reachable because the rhov loop is now the Kubota closed form.

## 4-for-4, all preregistered predictions HIT

| # | prediction | outcome |
|---|---|---|
| P1 | support {(165,1),(165,2)}, q = D/2 = 165, s = 2, NO third term | **HIT** -- the full fit over all definite (D',R) with D'R \| 330 returns exactly -2 GROSS(165,1) + 3 GROSS(165,2) mod cusp |
| P2 | weights (-2,+3) (2x the s=2 prime law), NOT the (-10/3,+5) of the constant-\|w\|/mass reading | **HIT** -- (-2,+3), and the competing reading is REFUTED |
| P3 | \|w2/w1\| = mass(165,2)/mass(165,1) = 3/2 (masses 5/12, 5/8 computed in advance) | **HIT** -- 11th consecutive base for the mass-ratio law |
| P4 | kernelrat rank 15, aliasing dim 25 (the composite-D values) | **HIT** -- rank 15, aliasing 25 |

## The falsifier fired correctly (gauge checks, restrict330A/B)

* Competing "largest ramified prime" support {(11,1),(11,30)}:
  **EXCLUDED**, rank 135 vs 136 -- E_eis is NOT in that span.
* Preregistered support {(165,1),(165,2)} alone: rank 135 vs 135, E in
  span TRUE, weights (-2,+3) with **aliasing dim 0** -- the weights are
  UNIQUELY determined here, with no gauge freedom whatsoever.  (At 210_1
  the corresponding check left aliasing; this base pins them outright.)

## What the two composite bases jointly establish

|  | mass1 | mass2 | mass ratio | w | \|w\|/mass |
|---|---|---|---|---|---|
| 210_1, q=105 | 1/4 | 3/8 | 3/2 | (-2,+3) | 8 |
| 330_1, q=165 | 5/12 | 5/8 | 3/2 | (-2,+3) | 24/5 |

The WEIGHTS are identical across the two composite bases while
\|w\|/mass is not (8 vs 24/5) -- so the composite-D weight is the
constant (-2,+3) = 2x the s=2 prime law, and the constant-\|w\|/mass
reading dies.  The mass RATIO law survives at 11 bases.  Per-prime
multiplicativity of the masses also replicates: the 7 -> 11 swap
(105 -> 165) multiplies both masses by (11-1)/(7-1) = 5/3, predicted in
advance and confirmed.

## What this base deliberately does NOT decide (stated in advance)

"Odd part" vs "largest definite divisor" CANNOT separate here -- they
coincide identically on every even squarefree D.  This base is a
REPLICATION at fresh primes ({3,5,11} vs {3,5,7}), not a rule separator.
The genuine separator remains an odd D with four prime factors (minimal
D = 1155, M = 4620), blocked on (i) the non-elementary 2-part -- Z/4
Jordan components need Stromberg Def 2.15 extended in eis32k -- and (ii)
the O(\|D\|) YTab scans at \|D\| ~ 1e7.

Open, unchanged: the scale 2 itself (now measured at TWO composite bases,
still underived) and the 2^(omega-1) genus-mass factor.
