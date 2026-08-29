# Non-integral principal parts: a predictor of failure, and where they come from

2026-08-29. Tools: `ppint.m` (integrality of every form's principal part, both cusps),
`whichclause.m` (replays a logged Hauptmodul table through `HauptmodulM0Residuals`).
Logs under `vvdata/weyl-campaign/ppint/`.

Follows [[base-39-2-malformed-form]], which found non-integral principal parts at `39_2`
and `46_3` and asked whether that predicts the backlog's failures.

## The measurement

`ppint.m` dumps each Borcherds form's principal part and reports whether the coefficients
are integral -- Schofer / Guo-Yang Thm B assume `c_eta(m) in Z` for `m <= 0`. It skips
`ValuesAtCMPoints`, so it runs on bases where the pipeline dies.

| base | integrality | pipeline outcome |
|------|-------------|------------------|
| 15_2 | INTEGRAL (both cusps) | works, 49/49 ground truth |
| 21_2 | INTEGRAL (both cusps) | **works: residual 0, mult (2,2)** -- measured here |
| 22_3 | INTEGRAL | **works: residual 0, mult (3,3)** -- measured here |
| 22_5 | INTEGRAL | works, 41/41 |
| 26_3 | INTEGRAL | not swept |
| 34_3 | INTEGRAL | not swept |
| 35_2 | INTEGRAL | works, residual 0, mult (0,0) |
| 38_3 | INTEGRAL | **works: residual 0, mult (-1,-1)** -- measured here |
| 6_29 | INTEGRAL | known-working model (Tier 2 control) |
| **33_2** | **NONINTEGRAL** (keys -1,-2 only) | **fails: table inconsistent** |
| **39_2** | **NONINTEGRAL** | **fails: M0MultiplierExact crash** |
| **46_3** | **NONINTEGRAL** | closed form gives 117/8, 3/8 -- not half-integers |

No counterexample in either direction. The claim is a **correlation with a mechanism**, not
a theorem: nine integral bases, three non-integral ones.

## 21_2 was never actually retested

[[handoff-2026-08-29-evening]] recorded "Tier 0: six for six no longer crash in
RationalNumber". The `21_2` log backing that is dated **Aug 21**, one day BEFORE `ca4d9fd`
shipped `M0MultiplierExact`; only five bases were re-run on Aug 29. Re-run here: `21_2`
sweeps clean in 117 s, residual 0. The claim is now true, but it was not evidenced before.

## 33_2's "NOT MEASURABLE" is a DEFECT, not a sweep limitation

`gtsweep` prints that message whenever `HauptmodulM0Residuals` returns empty, and the
message names two possible causes without saying which fired:

* (a) fewer than two informative points -- an early return at `SchoferFormula.m:1272`;
* (b) no `(r1,r2)` in `[-12..12]^2` fits the table -- a genuine inconsistency.

`whichclause.m` replays the logged table: `33_2` has **4 informative rows**, so (a) did not
fire. The row at `d = -15` is unsatisfiable for **every** `r1`:

    a = (9/4)*2^r1,  b = 5/16,  need  g1*a + g2*b = 1  =>  2^r1 in {11/36, 7/12}

Neither is a power of 2. And `d = -15` is an **anchor** (a divisor point of the Hauptmodul
forms), not one of the three discriminants admitted through the `Keep` hook -- so this is
not an artifact of the sweep's own widening. Contrast `46_3`, whose 25 surviving pairs are
the genuinely underdetermined case, and `35_2`/`21_2`, which have only 2 informative rows
each and still resolve to a single pair.

## Where the non-integrality comes from -- and why solving over Z will not fix it

Probe at `BorcherdsForms.m:817` (`Solution(coeffs_trunc, target_v)`), reporting the kernel
dimension, the denominator of the returned solution, and whether an **integral** solution
exists (`IsConsistent` over `Z` after clearing denominators -- verified to solve over `Z`,
not `Q`, on a system with a rational-only solution):

    33_2:  keys 9..15   kerdim 19-42   solden 1   intsol TRUE
           keys -1,-2   kerdim 18      solden 2   intsol FALSE
    15_2:  all keys     kerdim 8-22    solden 1   intsol TRUE

Two things follow.

1. **The solution space is large and `Solution()` returns one arbitrary point of it.** At
   `33_2` keys `-1,-2` that is an 18-dimensional affine space.
2. **But the family contains no integral point at all.** So this is *not* a matter of
   picking a better representative, and replacing the solve with an integral one would
   simply fail. *(I initially guessed the opposite -- that a positive-dimensional kernel
   meant the non-integrality was an artifact of the choice of point. `intsol false`
   refuted it.)*

The defect is therefore in the **divisor**, and it lands exactly where the divisor is a
free choice. Keys `9..15` have divisors forced by the covers' ramification and are integral
everywhere measured; keys `-1,-2` get their divisors from two CM points picked arbitrarily
at `BorcherdsForms.m:701-702`, inside a `CartesianPower` search that **breaks at the first
choice admitting a rational solution** and never checks integrality.

## The fix, PROTOTYPED AND WORKING

Add integrality as an acceptance criterion in that search: reject an `(infty, other_pts)`
triple whose solve has no integral solution and try the next one, and take the solution
from `IsConsistent` over `Z` rather than `Solution` over `Q`.

Done in the campaign copy of `BorcherdsForms.m` (~8 lines at the `Solution` call).
**`33_2` now reports `oo:INTEGRAL all:INTEGRAL`.** The trace shows the search walking past
the bad divisors instead of committing to the first one (`fix_33_2.log`, 105 PROBE lines):

    key -1:  ram [-12,-4]  intsol false -> rejected
             ram [-88,-4]  intsol false -> rejected
             ram [-132,-4] intsol TRUE  -> accepted
    key -2:  ram [-12,-4]  intsol false -> rejected
             ram [-88,-4]  intsol false -> rejected
             ram [-148,-4] intsol TRUE  -> accepted

So a valid choice existed all along. The old search accepted the first divisor admitting a
merely RATIONAL solution and never looked further. Keys `9..15` are `intsol true` on every
triple, confirming again that only the free-choice divisors were ever at fault.

Two caveats on the prototype. It recomputes keys `9..15` identically for every triple
(visible as the repetition in the log) -- correct but wasteful, worth hoisting before this
ships. And the probe `printf` is unconditional; it needs a `vprintf` before merging.

## BUT INTEGRALITY IS NOT SUFFICIENT -- 33_2 still fails

Sweeping `33_2` against the fixed forms (`gtsweep_33_2_fixed.log`):

    Runtime error in 'RationalNumber': s does not represent a rational number!

at `gtsweep.m:52`, i.e. BEFORE `HauptmodulM0Residuals`. The failure MOVED, it did not
clear. Unfixed `33_2` produced a rational but inconsistent table; fixed `33_2` produces a
table whose entries are not rational at all.

The prototype is not the problem: the divisor assertion at `BorcherdsForms.m:828` passes,
so the forms have the correct divisors AND integral principal parts. What this refutes is
the stronger reading -- **non-integrality predicts failure, but integrality does not
guarantee success.** `33_2` carries a second, independent defect.

Note the table is not comparable row-for-row with the old one: the new Hauptmodul divisors
are `-132` and `-148`, so those became anchors and the `Keep`-admitted extras changed from
`[-132,-148,-168]` to `[-88,-168,-232]`. The meaningful question was only whether the new
table is self-consistent, and it never got that far.

**Next diagnostic (cheap, not yet done):** print the raw `abs_tab` values before
`RationalNumber` and check whether they are quadratic irrationalities -- an `s` carrying a
`sqrt(2)` would point at a half-integral m=0 multiplier for the new forms, which is a
different bug from the one fixed here. Do NOT assume that without measuring it.

Also not established: that a sweep passing means a model can be BUILT. `gtsweep` drives the
pipeline only as far as the Hauptmodul table and the multiplier; it validates the m=0 term,
not model generation. `22_3`'s `(3,3)` says the multiplier is right there, nothing more.

Related: [[base-39-2-malformed-form]], [[post-m0-backlog-sweep-plan]],
[[handoff-2026-08-29-evening]].
