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

## 39_2 REFUTES THE FIX, AND NARROWS THE ROOT-CAUSE STORY

Running the fix on `39_2` (`ppint/fix_39_2.log`):

    PROBE key -1  ram [-52,-24]   intsol true  solden 1
    PROBE key -2  ram [-312,-24]  intsol true  solden 16384
    PROBE key 9   ...             intsol true  solden 8192
    ... intsol TRUE on every key ...

    VERD -2 oo 78 13 denom 4  1e48  NONINTEGRAL
    VERD -1 oo 14 13 denom 32 1e53  NONINTEGRAL
    ... every key oo:NONINTEGRAL, every key 0:INTEGRAL ...
    BASEVERD 39 2 oo:NONINTEGRAL all:NONINTEGRAL

Two things follow, both correcting what is written above.

1. **`intsol` is the WRONG CRITERION.** An integral solution vector in the eta basis does
   not give an integral principal part -- the echelon basis `ech_etas` has non-integral
   q-expansions, so the map `sol -> principal part` is not integral. At `33_2` the two
   denominators happened to coincide, which is the only reason the fix looked right there.
   The correct condition is integrality of the **q-expansion principal part**: find `x`
   with `x*coeffs_trunc = target_v` AND `x*P` integral, where `P` sends a basis vector to
   its principal-part coefficients. Still a lattice problem, but a different one.

2. **The free-choice-divisor story is 33_2-specific.** At `39_2` EVERY key is
   oo-non-integral, including `9..15` whose divisors are forced by ramification. So there
   is no arbitrary choice to blame; the whole weakly-holomorphic basis is bad, at magnitude
   1e48. `39_2` and `33_2` are two different failure modes that both surface as
   "non-integral principal part".

**But 46_3 IS rescued.** Running the fix on `46_3` gives `oo:INTEGRAL all:INTEGRAL` --
`solden 8` at keys `-1,-2`, replaced by the integral solution, and the principal parts come
out integral. So the fix repairs **2 of the 3** non-integral bases:

    33_2  NONINTEGRAL -> INTEGRAL   (base still fails downstream, separate defect)
    46_3  NONINTEGRAL -> INTEGRAL
    39_2  NONINTEGRAL -> NONINTEGRAL

and the two failure modes separate cleanly:

* **33_2, 46_3** -- an integral coefficient vector DOES give an integral principal part, so
  `intsol` is the right criterion and the fix works.
* **39_2** -- it does not, because the weakly-holomorphic BASIS has non-integral
  q-expansions (denominators 4 and 32, magnitude 1e48). No choice of coefficient vector can
  repair that; the target is `basis_of_weakly_holomorphic_forms` / `WeaklyHolomorphicBasis`.

The predictor itself is untouched: non-integrality still tracks failure. What is refuted is
that `intsol` is the criterion IN GENERAL, and that one cause explains every base.

## WARNING: THE CAMPAIGN BRANCH HAS NO M0MultiplierExact

`grep -c "intrinsic M0MultiplierExact" VectorValuedForm.m` is **0 on this branch** and 1 on
`tier1-models`. So a `gtsweep` run launched from THIS worktree applies a different m=0 term,
and dies outright at `gtsweep.m:114` if it ever reaches the multiplier stage. That is why
the previous session ran gtsweep from the main worktree and merely wrote its logs here.

This **confounds every fixed-form comparison below that was run from the campaign worktree**
(`gtsweep_33_2_fixed.log`, `gtsweep_34_3_mn14.log`, `gtsweep_46_3_fixed.log`,
`gtsweep_22_3_mn10.log`), because the fix lives here but the multiplier does not:

* "33_2 still fails downstream in RationalNumber, a second independent defect" -- **suspect.**
  A `RationalNumber` failure is exactly the signature of a MISSING m=0 term, which is what
  this branch has. Re-running with the patch applied to the main worktree instead.
* "46_3's r(-1) moved 0 -> 3 with the fixed forms" -- **not established**; the applied
  multiplier differed too.
* the uniform `(-3,-3)` re-basing at `34_3`/`22_3`, attributed above to `MaxNum` -- more
  simply explained by the absent `M0MultiplierExact`.

Runs launched from the MAIN worktree (`21_2`, `22_3`, `38_3`, `26_3`, `34_3`) are unaffected.
**Lesson: the two worktrees are not interchangeable packages. Check which intrinsics a
branch actually defines before attributing a result to the thing you changed.**

## RE-RUN WITH THE CORRECT PACKAGE: THE FIX DOES REPAIR 33_2

Applying the prototype to the MAIN worktree (patch via `git diff`/`git apply`, reverted
after) and sweeping there:

* **`33_2` no longer crashes in `RationalNumber`** -- that crash was entirely the campaign
  branch's missing `M0MultiplierExact`. The "second independent defect" was my artifact.
* **`33_2` moved from clause (b) to clause (a)** -- i.e. from INCONSISTENT to merely
  under-sampled. Replayed through `whichclause.m`:

      d = -88   exps [0,0]  row [11/36, 25/36]  satisfiable   (sums to exactly 1)
      d = -168  exps [0,0]  row [1/2, 1/2]      satisfiable
      d = -232  exps [0,0]  row [1/18, 17/18]   satisfiable
      d = -15   exps [1,1]  row [1/4, 5/4]      satisfiable   (-1/4 + 5/4 = 1)

  Before the fix, `d = -15` was unsatisfiable for EVERY `r1`. Now **all four rows satisfy
  the find_signs relation exactly at (0,0)**, three of them as consistency checks that
  pass. The sweep still declines only because it needs TWO informative rows and there is
  one. So the multiplier is not formally determined, but the table is now consistent and
  everything points at `(0,0)`.
* **`46_3`** gives 25 pairs all with `r(-1) = 0` -- identical to the ORIGINAL unfixed run,
  so the fix leaves its residual unchanged. The `0 -> 3` shift reported earlier was the
  package confound, not the forms.

So the corrected scorecard for the fix: **33_2 repaired** (inconsistent -> consistent),
**46_3 well-formed but unchanged**, **39_2 untouched** (different failure mode).

## WHAT REMAINS: NOT ENOUGH DEGREE-1 FIRING DISCRIMINANTS

`33_2`, `46_3` and `34_3` now share one limitation:

`46_3`'s underdetermination is STRUCTURAL, not a sampling shortfall. With `i_s0 = -3`
(non-firing, `kB = 0`) and `i_st0 = -8` (firing, `kA = 1`), the three usable degree-1 rows
are `-24, -312, -372`, all non-firing, all with `exps = [k-kA, k-kB] = [-1, 0]`. **`r2`
appears in no row at all**, so it is free by construction -- hence exactly 25 surviving
values, one per allowed exponent. The only firing degree-1 discriminant is the normaliser
`-8` itself; `-35` fires but is degree 2, and `HauptmodulM0Residuals` skips `degs ne 1`.

`34_3` has the same shape (2 informative rows, both `[1,1]`), and repaired `33_2` has ONE
informative row (`-15`), its only other firing discriminant `-55` being degree 2. So the
lever for all three is **degree-1 FIRING discriminants specifically** -- and `MN`/`MaxNum`
cannot ask for that: at `34_3` it spent the entire extra budget on degree-2 points and
changed nothing. Selecting for degree-1 firing discriminants is the concrete fix to the
sweep, and it would potentially determine three bases at once.

Note also that the 0-side is INTEGRAL throughout at `39_2` while the oo-side is not --
which is why splitting the verdict into `oo:` and `all:` was worth doing.

Also not established: that a sweep passing means a model can be BUILT. `gtsweep` drives the
pipeline only as far as the Hauptmodul table and the multiplier; it validates the m=0 term,
not model generation. `22_3`'s `(3,3)` says the multiplier is right there, nothing more.

Related: [[base-39-2-malformed-form]], [[post-m0-backlog-sweep-plan]],
[[handoff-2026-08-29-evening]].
