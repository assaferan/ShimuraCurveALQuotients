# PREREGISTRATION: X_0^1155(1), the four-odd-prime re-split test

Written before any part of the base is computed -- no pool, no rho, no
fit exists for `M = 4620`. Engineering scoping is in `scope-1155.md`;
this file is only the predictions.

`D = 1155 = 3 * 5 * 7 * 11`, `N = 1`, level `4DN = 4620`,
`|A| = 2,668,050`.

## Why this base

At every even `D` the three candidate re-split rules coincide:

    odd part of D  =  D with its smallest prime (2) removed
                   =  largest divisor that is a definite discriminant

`210_1` and `330_1` confirmed that common value (`q = 105`, `q = 165`)
but **could not separate the three** -- recorded at the time.

`1155` has no prime `2`. Its odd part is `D` itself, which has four
prime factors and is therefore **indefinite**: not a discriminant of a
definite quaternion algebra, so `theta(1155, R)` does not exist. The
"odd part" rule makes no prediction here at all, while the other two
still do. This is the first base where the readings come apart.

**Recorded in advance, honestly:** `1155` does NOT separate "remove the
smallest prime" from "largest definite divisor" -- those agree on every
squarefree `D` with an even number of primes, hence on every base this
campaign can reach. That distinction is not testable here.

## The weights are not free -- two banked laws fix them

At `N = 1` the fit is read at the zero coset (Remark "the tracked
coset"), where every banked base has **weight sum 1**; and the
**mass-ratio law** gives `|w2/w1| = rho`, the mass ratio of the support.
Together these determine the pair with no freedom:

    w1 = 1/(1 - rho),    w2 = rho/(rho - 1)

Checked against all ten banked `N=1` fits: `rho = 3/2 -> (-2,+3)`
(210_1, 330_1, and the four `D=2q` bases); `rho = 2 -> (-1,+2)`
(15_1, 21_1, 33_1); `rho = 3 -> (-1/2,+3/2)` (35_1). No exceptions.

So the measurement at `1155` is a **classification**: which support,
and the weights follow. Every candidate has a distinct answer.

## The candidates, with masses computed in advance

All ten masses measured (`mass1155.m`, this commit), not inferred:

| rule | `D'` | `s = D/D'` | masses | `rho` | predicted `(w1,w2)` |
|---|---|---|---|---|---|
| **largest definite divisor / drop smallest prime** | **385** | **3** | 5/4, 5/2 | **2** | **(-1, +2)** |
| drop the largest prime | 105 | 11 | 1/4, 3/2 | 6 | (-1/5, +6/5) |
| largest ramified prime | 11 | 105 | 5/24, 5 | 24 | (-1/23, +24/23) |
| drop 5 | 231 | 5 | 5/8, 15/8 | 3 | (-1/2, +3/2) |
| drop 7 | 165 | 7 | 5/12, 5/3 | 4 | (-1/3, +4/3) |
| odd part | 1155 | -- | -- | -- | **no prediction: indefinite** |

Every mass obeys
`mass(D',R) = prod_{p|D'}(p-1) * prod_{p|R}((p+1)/2) / (48 * 2^(w(D')-1))`,
the formula fitted to the sixteen banked masses, including the
`2^(omega-1)` composite factor.

## THE PREDICTIONS

**Q1 (primary).** Support `{(385,1),(385,3)}`, weights `(-1,+2)`,
`|w2/w1| = 2`. This is "largest definite divisor", the rule that
survived `210_1` and `330_1`.

**Q2.** No third term (`N=1`, zero coset -- as at all six banked `N=1`
bases).

**Q3.** Weight sum exactly `1`.

**Q4 (the falsifier, and the reason this base is worth its cost).** The
four losing candidate spans **FAIL** the rank test: `E_eis` does not lie
in `span{theta(D',1), theta(D',s)} + cusp` for `D' = 105, 11, 231, 165`.
Basis for expecting exclusivity rather than the two-prime degeneracy:
at four primes `210_1` excluded its competitor at rank 87 vs 88 and
`330_1` at 135 vs 136, whereas every two-prime base is degenerate
(Observation "the odd-level quartet"). If instead several supports
carry `E_eis` here, then four-prime selection was an artifact of the
even-`D` cases and the whole re-split rule is a convention -- which
would be the most informative outcome of the campaign so far.

**Q5.** The residue pairing holds: `rho` lies in the holomorphic
weight-3/2 eta span, residual `~1e-42` at working precision. This has
held on all nineteen bases; a failure would indicate a driver defect at
odd level rather than arithmetic, and the `CTL` gate should be run
first to distinguish them.

**Q6.** The mass-ratio law holds on whatever support wins (test 20).
Note Q1 and Q6 are independent here: Q6 is satisfied by ANY of the five
candidates, since the predicted weight pair was derived from that
candidate's own mass ratio. Q6 fails only if the measured `|w2/w1|`
matches no candidate's mass ratio.

**No prediction** is offered for the Eisenstein dimension, constants
rank, or aliasing at level 4620: those are per-level-family constants
(Observation 9.2), and 4620 opens a new family.

## What each outcome means

* **Q1 + Q4 hit**: the re-split rule is "largest definite divisor",
  confirmed where it is not equivalent to "odd part", and the campaign's
  central empirical law is complete up to the smallest-prime ambiguity
  recorded above.
* **A different support wins with Q4 intact**: the rule is one of the
  other readings, and `210`/`330` were consistent with it by accident of
  having `2` as the removed prime.
* **Q4 fails (several supports carry it)**: four-prime selection was an
  even-`D` artifact; the support is a convention at every `D`, and only
  the mass-ratio law survives as content.
* **Q6 fails**: first failure of the mass law in twenty bases.
