# PREREGISTRATION: the odd-level N=1 bases 15_1, 21_1, 33_1, 35_1

Written before running the Gross-span recognition on any of them.
Supersedes nothing; opens the odd-`DN` regime (level `4DN`).

## Why these bases

Every one of the six banked `N=1` bases has **s = 2**: the four small
ones are `D = 2q` outright, and `210_1` / `330_1` re-split as `q = D/2`,
`s = 2`. So in the banked data

> "the N=1 weights are the constant (-2,+3)"

and

> "the N=1 weights are 2x the s-law, i.e. (-2/(s-1), (s+1)/(s-1))"

are **perfectly confounded** -- they agree at s = 2 and nowhere else.
This is the same confound structure as composite-q vs N=1 (killed by the
D=2q quartet) and it needs the same treatment: find the missing cell.

The missing cell is `N = 1` with `s > 2`, which forces **odd D**, hence
**odd DN**, hence level `4DN` instead of `2DN` -- a regime the pipeline
had never entered. Four cheap bases live there:

| base | q | s | M = 4DN | \|A\| |
|---|---|---|---|---|
| 15_1 | 5 | 3 | 60 | 450 |
| 21_1 | 7 | 3 | 84 | 882 |
| 33_1 | 11 | 3 | 132 | 2178 |
| 35_1 | 7 | 5 | 140 | 2450 |

15_1 shares its level (60) with 15_2, so the banked `epool_15_2.txt`
serves it unchanged.

## Driver status (settled BEFORE any fit, so it is not tuned to the answer)

Odd `DN` gives a 2-part of the discriminant form that is `Z/2` of ODD
type (`4*Qbar = 1 mod 4`), uniformly on all four bases -- not the
`(Z/2)^3` of the even-`DN` family, and NOT the `Z/4` component that was
feared. `eis32k.m` gains an `r2 = 1` branch (`x_c` = the single
generator); the `r2 = 3` search is left byte-for-byte as banked.

New `CTL:=1` mode recomputes every `rhov` entry by FFT and reports the
worst deviation. Both gates run and passed already:

* 15_2 (r2 = 3, regression): worst `|closed - FFT| = 7.25e-80`, 144/144 words.
* 15_1 (r2 = 1, the new path): worst `1.32e-79`, 144/144 words, 120 nonzero.

15_1 also already solves: `SOLVE resid = 8.00e-42` against
`|rho| = 2.1909`, so the weight-3/2 residue pairing holds at odd level.
**The Gross-average decomposition has NOT been run on any of the four.**

## Masses, computed in advance

| base | canonical support | masses | ratio | competitor | masses | ratio |
|---|---|---|---|---|---|---|
| 15_1 | (5,1),(5,3) | 1/12, 1/6 | **2** | (3,1),(3,5) | 1/24, 1/8 | 3 |
| 21_1 | (7,1),(7,3) | 1/8, 1/4 | **2** | (3,1),(3,7) | 1/24, 1/6 | 4 |
| 33_1 | (11,1),(11,3) | 5/24, 5/12 | **2** | (3,1),(3,11) | 1/24, 1/4 | 6 |
| 35_1 | (7,1),(7,5) | 1/8, 3/8 | **3** | (5,1),(5,7) | 1/12, 1/3 | 4 |

All match `mass(q,R) = (q-1)(R+1)/96` (omega = 1, so no `2^(omega-1)`
factor). Each ratio equals `(R+1)/2` for the relevant second level.

## THE PREDICTIONS

**P1 (primary -- 2x the s-law).** Support `{(q,1),(q,s)}` with `q` the
larger ramified prime, no third term, weight sum 1, and

    w = ( -2/(s-1),  (s+1)/(s-1) )

so `15_1, 21_1, 33_1 -> (-1, +2)` and `35_1 -> (-1/2, +3/2)`.

**P2 (the rival -- constant N=1 weight).** All four give `(-2, +3)`,
the weights of all six banked `N=1` bases.

P1 and P2 differ on every one of these four bases and agree nowhere in
the banked data. **Discriminator: `|w2/w1|` = 2, 2, 2, 3 under P1;
= 3/2 under P2.**

**P3 (the mass-ratio law, 16th through 19th tests).** `|w2/w1|` equals
the mass ratio of whatever support carries `E_eis`. Note P3 implies P1
and REFUTES P2 -- so P2 winning would be the first failure of the mass
law in 15 bases. This is the sharpest thing at stake here.

**P4 (support degeneracy).** At `D = 2q, N=1` the support was degenerate
(the competitor `{(2,1),(2,q)}` also carried `E_eis`, at its own mass
ratio). Prediction: at odd `D` the support is NOT degenerate and the
larger prime wins -- i.e. the competitor span in the table above FAILS
the rank test. Recorded as the one I am least sure of; if it is
degenerate instead, P3 must still hold on both supports.

**P5 (the level constants).** 15_1 has the same level (60) as 15_2, so
Eisenstein dimension 11 and constants-rank 7 should recur there;
aliasing may differ, since the Gross pool is built from `DN = 15` rather
than 30. No prediction for the other three levels (84, 132, 140).

## What each outcome means

* **P1 hits**: the `N=1` doubling is uniformly `2 x (s-law)`, the coset
  account of Remark "the tracked coset" is complete, and the s-law is
  one law across both normalisations and both parities of `DN`.
* **P2 hits**: the mass-ratio law fails for the first time, and the
  `N=1` regime is not a rescaling of the `N>1` one -- which would
  re-open the scale question that the `EST:=0` experiment closed.
* **Support elsewhere** (neither pair): the re-split rule is more
  general than "largest ramified prime" even at two primes, which would
  bear directly on the `1155_1` question these bases are a warm-up for.
