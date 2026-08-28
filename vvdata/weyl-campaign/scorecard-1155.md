# SCORECARD: X_0^1155(1), measured against `prereg-1155.md`

Preregistration banked `c4e99df`, before any part of the base was computed.
Nothing below was read against the prereg until the numbers existed.

## The measurements

The fit (`eis32s`, PR=120 with `etaRed`, 354-vector merged pool):

    #E's after external pool = 354
    E constant table built (coset-outer, 276.87 s)
    SOLVE resid = 2.217329371E-40  (|rho| = 2.446306519)
    RHO IS IN THE SPAN: E* found          285 nonzero BETA rows of 354

`eisrat`: **RATIONAL E_eis FOUND**, rank 31, residual beyond rank 1.75e-46,
pivot snap error 1.08e-44.  (The Gross driver independently reports
`Eis 31` at level 4620 -- an unplanned cross-check that agrees.)

Q4, by the constant-term route (`ctq4_1155.m`; see below for why not the
q-expansion route):

| support | weights | sum | \|w2/w1\| | E_eis in span |
|---|---|---|---|---|
| **(385,1)+(385,3)** | (-1, +2) | 1 | 2 | **yes** |
| (105,1)+(105,11) | (-1/5, +6/5) | 1 | 6 | **yes** |
| (11,1)+(11,105) | -- | -- | -- | **no** (rel. resid 0.961) |
| (231,1)+(231,5) | (-1/2, +3/2) | 1 | 3 | **yes** |
| (165,1)+(165,7) | (-1/3, +4/3) | 1 | 4 | **yes** |

all to 115 digits, relative residuals 1.7e-114 on the four that pass.

## The verdicts

**Q1 (primary) -- NOT DECIDED AS INTENDED.**  Support {(385,1),(385,3)} does
carry E_eis with exactly the predicted weights (-1,+2).  But so do three of
the four supports the prereg expected to fail.  Q1 asked which support is
selected; the answer is that the test does not select one.

**Q2 -- MET.**  No third term: every fit is exhausted by its two theta terms.

**Q3 -- MET.**  Weight sum exactly 1 on all four passing supports, to 115
digits.

**Q4 (the falsifier) -- FAILS.**  The four losing spans do NOT all fail.
Only `D' = 11` fails; 105, 231 and 165 pass with their own predicted weights.

**Q5 -- MET.**  The residue pairing holds: rho lies in the holomorphic
weight-3/2 eta span, relative residual 9.1e-41 against a 1e-20 criterion.

**Q6 -- MET.**  `|w2/w1| = rho` on every support that passes, exactly.  As
the prereg noted in advance, Q6 is independent of Q1 and is satisfied by any
candidate, since each predicted pair was derived from that candidate's own
mass ratio.

## What Q4's failure actually means -- and it is NOT what the prereg guessed

The prereg's contingency reading was: *"four-prime selection was an artifact
of the even-D cases, and the support is a convention at every D."*  The first
half is wrong, and in an instructive way.

Running the same test on the earlier bases shows the SAME non-uniqueness
there.  At 330_1 (`ctcontrol330.m`), of the eight supports
{(D',1),(D',330/D')}:

    (2,1)+(2,165)    rel resid 0.955   fail        <- D' one prime
    (3,1)+(3,110)    rel resid 0.980   fail        <- D' one prime
    (5,1)+(5,66)     rel resid 0.988   fail        <- D' one prime
    (11,1)+(11,30)   rel resid 0.990   fail        <- D' one prime
    (30,1)+(30,11)   3.4e-115  PASS  w=(-1/5,+6/5)  ratio 6
    (66,1)+(66,5)    3.4e-115  PASS  w=(-1/2,+3/2)  ratio 3
    (110,1)+(110,3)  3.4e-115  PASS  w=(-1,+2)      ratio 2
    (165,1)+(165,2)  3.4e-115  PASS  w=(-2,+3)      ratio 3/2

The same at 210_1 (`ctcontrol210.m`), of {(D',1),(D',210/D')}:

    (2,1)+(2,105)    rel resid 0.961   fail        <- D' one prime
    (3,1)+(3,70)     rel resid 0.982   fail        <- D' one prime
    (5,1)+(5,42)     rel resid 0.989   fail        <- D' one prime
    (7,1)+(7,30)     rel resid 0.990   fail        <- D' one prime
    (30,1)+(30,7)    1.0e-115  PASS  w=(-1/3,+4/3)  ratio 4
    (42,1)+(42,5)    1.0e-115  PASS  w=(-1/2,+3/2)  ratio 3
    (70,1)+(70,3)    1.0e-115  PASS  w=(-1,+2)      ratio 2
    (105,1)+(105,2)  1.0e-115  PASS  w=(-2,+3)      ratio 3/2

So the pattern is not about four primes at all, and it holds at THREE bases
-- 210_1, 330_1, 1155_1 -- identically:

    D' with ONE prime      -> FAILS
    D' with THREE primes   -> PASSES, weights fixed by its own mass ratio

**Why this was not seen before.**  The banked campaign tested exactly ONE
competitor per base, and in both cases it was a single-prime D':

    restrict330A:  (11,1)+(11,30)   rank 135 vs 136  ->  excluded
    restrict210A:  (7,1)+(7,30)     rank  87 vs  88  ->  excluded

Those are precisely the supports that fail here too, so the banked
exclusions are CORRECT -- they were simply not the discriminating test they
were taken for.  The composite competitors (30, 66, 110 at 330_1) were never
tested.  The claim "at four primes 210_1 and 330_1 genuinely exclude the
competitor" is true as far as it goes and does not support "four primes
select a support".

## The result is verified against the original method, not just the new one

The constant-term route is new, so it does not get to overturn a banked claim
on its own authority.  Two independent confirmations:

1. **It reproduces the banked answer.**  At 330_1 it returns
   ct(E_eis) = -2 ct(theta(165,1)) + 3 ct(theta(165,2)) exactly: residual
   7.87e-115, 0 of 3456 cosets missing.
2. **Head-to-head on a support the banked campaign never tested.**  The
   ORIGINAL q-expansion machinery, run on (110,1)+(110,3) at 330_1 with
   DEPTH=400:

       Gross-average+cusp span: rank 135 vs 135 -> E in span: true
       E_eis = -1 * GROSS(110,1) + 2 * GROSS(110,3)
       aliasing dim 0

   Identical weights, and `aliasing dim 0` says the decomposition is unique
   for that support.

   Extended to ALL FOUR passing supports at 330_1, every one by the original
   q-expansion machinery, every one on a support the campaign never tested:

       (30,1)+(30,11)   rank 135 vs 135  true  -1/5 * G(30,1) + 6/5 * G(30,11)
       (66,1)+(66,5)    rank 135 vs 135  true  -1/2 * G(66,1) + 3/2 * G(66,5)
       (110,1)+(110,3)  rank 135 vs 135  true   -1  * G(110,1) +  2  * G(110,3)
       (165,1)+(165,2)  rank 135 vs 135  true   -2  * G(165,1) +  3  * G(165,2)   [banked]

   all with `aliasing dim 0` and all matching the constant-term weights
   exactly.  So the non-uniqueness of the SUPPORT is confirmed four times over
   by the method the campaign has used all along, and does not rest on the new
   route.  Note what `aliasing dim 0` means here: given a support the weights
   are UNIQUE.  The non-uniqueness is in WHICH SUPPORT, not in the fit.

Degeneracy was also ruled out before drawing any conclusion: the ten
ct(theta) vectors at 1155 have **rank 7 of 10**, so the four passing supports
are four genuinely distinct 2-planes, not one plane counted four times.  The
four (D',1) vectors are mutually orthogonal to ~1e-121 while each overlaps
its own partner strongly (0.87, 0.55, 0.75, 0.66).

## What survives

* the **mass-ratio law** -- `|w2/w1| = rho` and weight sum 1, now on many
  more supports than before, at 115 digits.  This is the content.
* the **one-prime / three-prime distinction**: E_eis is carried by every
  three-prime D' and by no single-prime D'.  That is a real, reproducible
  fact, and it is what the banked exclusions were actually measuring.

## What does not survive

* **"the re-split rule"** as a selection principle.  There is no evidence the
  Gross-span test picks 385 over 105, 231 or 165 at 1155, or 165 over 30, 66
  or 110 at 330_1.  "q = the larger ramified prime", "odd part", "drop the
  smallest prime", "largest definite divisor" -- these agree on a value that
  the test does not single out.  The scale-2 note already recorded one
  convention masquerading as arithmetic ([[scale-2-is-forced]]); this is a
  second, and a larger one.

## Caveat, stated plainly

`aliasing dim 0` at (110,1)+(110,3) says the WEIGHTS are unique given that
support.  Nothing here bears on whether some finer invariant -- outside the
Gross-average-plus-cusp span -- distinguishes the supports.  The claim is
only that THIS test does not.
