# Handoff — 2026-09-13

**The newest section is this one; everything after it is older and kept for provenance.** Earlier
material still says things like "34 of 43" or "23 of 34 tests check involutions" — those counts are
STALE.

✅ **2026-09-13: everything is COMMITTED AND PUSHED on both branches**, and the branch-divergence
invariant prints nothing against `origin`. ⚠ lava's clone is still stale at `8dac84c` — `git fetch
&& git reset --hard origin/main` there before any new run, but NOT while a job is alive.

**➡ For what to do next, see `PLAN.md`.** This file records *what happened*; when the two disagree
about state, this file wins.

## Handoff — 2026-09-13 (later) — the post-condition-4 blocker, examined

**Supersedes the "blocker then MOVES to `QuadraticConstraintsOnEquations`" claim in the section
below.** Full account with every count: `vvdata/weyl-campaign/even-correction/QUADCONSTRAINTS.md`
(campaign, `9a949b4`), with `ratfit_compare.py` and `ratfit-cmextra.patch` beside it.

### ⚠ THE STAGE ATTRIBUTION WAS WRONG — read the `require`, not the traceback

    Runtime error in 'QuadraticConstraintsOnEquations':
    Error in Schofer table values at rational points - no solution found!

That `require` (`EquationsCovers.m:68`) tests `kernels[j]`, an ARGUMENT, built one intrinsic earlier
in `RationalConstraintsOnEquations` (`EquationsCovers.m:31`). When it fires the quadratic stage has
done no arithmetic. It also **cannot** be the blocker in principle: its relations are built from
`B[1]` and solved inside `P(B)`, so they select within an existing solution space and can never
create one. At `34_3` it is moot anyway — **`#quad = 0`**, so the stage is a no-op in both runs.

⇒ The real failure is the **rational linear fit**: no `f` of degree `<= 2g+2` matches the rational
CM values. Empty at all 7 cover keys, at FULL COLUMN RANK (baseline: `dimB = 1` at all 7).

### ✅ THE HATCH'S FOUNDING PREMISE HOLDS — and more strongly than "branch locus"

* branch locus preserved: **0 of 42** (rational CM point, key) cells change zero-ness;
* `h = y2_pert / y2_base` is a perfect **6th power at 58 of 58 cells** once the bad primes of
  `D*N = 102` are divided out, and its 6th root depends **only on `s`, not on the cover key**.

⇒ `h = c_key * G(s)^6`, one shared `G`, `c_key` a bad-prime constant (the re-chosen
`find_y2_scales` row scale). `amt = 6` is the exponent — exactly what adding `6*Z(164)` to the
divisor predicts. Predicted first, then measured.

⇒⇒ **`G^6 = (G^3)^2` is a perfect SQUARE, so `y^2 = c*f_old*(G^3)^2` is the SAME CURVE as
`y^2 = c*f_old`.** The even correction preserves the whole cover, not just its branch divisor, and
**the pipeline is discarding a correct answer** — its ansatz is simply too short.

### ⇒ THE BINDING QUANTITY IS A COST, NOT A FIFTH CONDITION

    degree gain    = amt * deg Z(disc)
    the fit needs    2g+5 + amt*deg Z(disc)   RATIONAL CM points

`deg Z(164) >= 2` at `34_3` (`= 1` refuted on an exactly-determined system at the `g=0` keys), so
the demand is `>= 19` against a supply of **10** — and 10 is genuine: `#ds = 19` still yields
`#rat = 10`, the six extra points all quadratic. **The hatch is not obstructed at `34_3`; it is
PRICED OUT.**

⚠ No selection rule so far considers `deg Z(disc)`, and `PROBE_EVEN` actively prefers the
**largest** `|disc|` — hence large class number, hence large `deg Z(disc)`. **That heuristic works
directly against the cost.** The cheapest legal correction is `amt = M` at a discriminant with
`deg Z(disc) = 1`.

### ⚠⚠ A THIRD NULL-RUN TRAP, same family as the two below — a knob at a dead call site

The first `CMEXTRA` knob went into `EquationsOfCovers` (the `[4/6]` path). **`genmodels.m` never
calls it** — `AllEquationsAboveCovers` has its own copy of the `num_vals` computation. Both runs
came back byte-identical in the same 207 s, which reads as "more CM points don't help": a clean
refutation of the degree story **from a run that added no points**. Caught because `[4/6]` appeared
in no log, baseline included. ⇒ The knob now prints `CMEXTRA num_vals` unconditionally.

⚠ Also: **the degree sweep is VACUOUS whenever `ncols > nrows`**, which it is at the default
`MaxNum = 7` (`#rat = 6` = `2g+4` for `g=1`, a square system). A kernel appearing above the true
bound there means nothing.

### ⚠ `genmodels.m` output is NOT byte-comparable to `data/models/`

The unperturbed run gives 15 cover keys against the committed 10, and shared entries differ — but
they are **the same curves in a different normalisation** (`W=[1,2,17,34]` fresh = committed/9;
`W=[1,17]` fresh(x) = committed(17x/3); `W=[1,102]` both). `genmodels.m` does not pin `base_label`.
A byte-diff is not a valid reproduction check for this driver.

## Handoff — 2026-09-13 (newest)

Continues 09-12. Everything below is committed AND PUSHED on both branches; the branch-divergence
invariant prints nothing against `origin`.

### ✅ CONDITION 4 — a fourth requirement on any even correction, and it is SATISFIABLE

The instrumentation asked for on 09-12 was built and run. `RationalNumber` fails **iff some prime
carries a NON-INTEGRAL exponent** (`LogSum.m:137`) -- a `LogSm` is a formal `sum_p coeff_p log p`,
so the failure is a fractional EXPONENT, not an irrationality, and naming the prime is the whole
diagnostic. Measured chain:

* the prime is always a RAMIFIED prime of `D` (17 at `D=34`, 5 at `D=35`), never the level prime;
* it is present BEFORE the final rescaling (`scale = -1/4`, denominator 4, wrong prime);
* the principal-part coefficients are INTEGERS, **362 of 362**;
* `Kappa0`'s OWN log-`p` coefficients are fractional -- 235 at denominator 3, 113 at 9, natural at
  `p=17` where `p+1=18`.

⇒ Those fractions are INTRINSIC and are supposed to cancel: a legitimate divisor makes
`sum_m c(-m) kappa_p(m)` an INTEGER, and the perturbation breaks that. So a usable perturbation
needs FOUR conditions, not three:

    1. EVEN                        -- cover unchanged                 (parity survey 28/28)
    2. phi(target) = 0             -- Borcherds' criterion            (always solvable, gcd(phi)=1)
    3. integral solution           -- a form exists at all            (170/834 candidates at 34_3)
    4. sum_m c(-m) kappa_p(m) in Z at every ramified p | D, every CM d   <-- NEW, and not implied by 3

**Condition 4 IS satisfiable** -- at `34_3`, `35_1` and `21_2` there are amounts giving ZERO
non-rational cells, clearing `ValuesAtCMPoints` for the first time since 2026-08-30. ⚠ The blocker
then MOVES to `QuadraticConstraintsOnEquations` ("Schofer table values at rational points -- no
solution found"), a NEW and UNEXAMINED stage that may be CM supply, a known rescue axis.

### ⚠ THE MODULUS LAW IS UNRESOLVED — SIX FORMULAS FITTED AND REFUTED

    D=34 (2*17)   M = 6  EXACTLY    (N=3 and N=7 AGREE -- so it is NOT N-dependent)
    D=35 (5*7)    M = 12 EXACTLY
    D=21 (3*7)    M | 6             (2,4 untestable there -- they fail condition 3)

Refuted: `2N` (at `34_7`), `2*oddpart(p+1)` for `p=17` (at `34_3`), `2(pmin+1)` (at `21_2`),
`2(pmax+1)` (at `34_3`), `2*gcd` (at `35_1`), `2*lcm` (at `34_3`). **Do not propose a seventh.**

⚠ **UNCONTROLLED CONFOUND, possibly making the question mis-posed**: the perturbation DISCRIMINANT
differs across the three bases (164, 32, 16), and nothing shows the modulus is a property of the
BASE rather than of the DISCRIMINANT. ⇒ The next experiment is the disc-dependence control -- vary
the disc at FIXED base -- because it decides whether "the modulus at base X" is even well defined.

⚠ **THE 28-BASE DIVISIBILITY SCREEN IS VOID** (it used `2N`). **The hatch's REACH IS UNKNOWN** --
neither "3-4 of 28" nor the original "49" is supported. Do not quote either.

### ⚠⚠ TWO NULL-RUN TRAPS, SAME CAUSE — always print the perturbed-key count

`PROBE_EVEN_COPRIME` silently excludes the chosen discriminant and perturbs NOTHING while returning
0 cells. Hit twice: `180/+4` at `34_3`, and THREE runs at `21_2` (`N=2`, disc 16 even) that read as
"all amounts clear" -- which would have refuted the `p_min` formula **on no computation at all**.
That filter rests on the non-coprimality hypothesis, which the even-correction README itself records
as REFUTED. ⇒ **Run with `PROBE_EVEN_COPRIME=0`, and count `perturb disc` lines before believing any
verdict.**

### ✅ tests/ConicClasses.m — the first check on a genus-0 twist class (`4f4667a`)

`ModelChecks`' four tests are STRUCTURALLY BLIND to a quadratic twist at genus 0: a conic and its
non-square twist share genus, genus formula, the trivial Weil polynomial AND the point count over
every `F_p` (every smooth conic over a finite field is isotropic, so both are `P^1` with `p+1`
points). **282 of 822 committed entries are genus 0 across 75 files** -- the largest exempted class
in the repo, previously unvalidated, and why `10_3`'s `[1,2]` drift was invisible to CI.

The new test needs no theory: entries under one `(D,N,W)` key are the same quotient over different
bases, so they must share a class in `Br(Q)[2]`, computed as the quaternion algebra `(a, disc)`.
**215 conics, 38 multi-entry keys, all consistent, 0.05 s.** NEGATIVE-CONTROLLED in situ: twisting
one `10_3` entry by `-2` makes it fail and name the key. Suite **76/76**.
⚠ It catches INTERNAL inconsistency only -- it does NOT resolve `10_3`, whose committed entries
agree with each other. Predicting WHICH conic is right is the open arbiter question; `W=[1]` entries
come out RAMIFIED, consistent with `X^D(R) = empty`, which points at Ogg's real-points criterion.

### Runtime on lovelace: where the time goes

Measured scaling (at `51_1`): pool grows LINEARLY in pole order, but `qexps ~ PO^3.3`,
`EchelonForm ~ PO^4.4`, `ech_etas ~ PO^2.5`. Extrapolated to the pole orders the live jobs reach
(`119_1` hit 1309, `111_1` 1665): ~0.6 h per pool build at `PO=1300`, ~1.7 h at `PO=1700`, and
`EchelonForm` OVERTAKES `qexps` around `PO~1200-1600` -- so these runs sit exactly at the crossover.
⚠ The recorded "EchelonForm is only 6 s of 222 s" is from a small pole order and does NOT
extrapolate.

Levers, ranked: (1) build the DEFICIT PREDICTOR -- the deficit is a rank comparison independent of
any divisor choice, so obstruction is computable from the Weil representation WITHOUT the pipeline;
nobody has built it, and it is the only lever that changes the campaign's complexity rather than one
run's. (2) Attack COEFFICIENT GROWTH in the elimination (multi-modular + CRT): the pool is already
100% triangular when sorted by valuation, so `EchelonForm` never searches for pivots -- but "skip
the RREF" is REFUTED, it only relocates 33-digit coefficients downstream. (3) Reuse q-expansions
across the t-ladder (`qexp(t^j f) = qexp(t)^j qexp(f)`), untested, needs higher absolute precision.
(4) Audit `Prec` per base -- precision is `M^2`.
⚠ The scaling table is ONE base; verify the exponents elsewhere before investing.

⚠ Checked: `bfp` (`95_1`) and `vxfix` (`159_1`) DO contain the vx fix `d9b52d0`, so those multi-day
runs are on valid code. All five lovelace jobs still alive; `34_11` at **81.6 GB** with 1450 GB free.

## Handoff — 2026-09-12

### ⇒ LIVE JOBS: `X0_111_1` IS DONE AND PASSED; the five on lovelace are still running

**`X0_111_1` SUCCEEDED on lava in 58595 s (16.3 h).** It survived the 1678-vector pool that the
last handoff flagged as near the ~2000-vector / ~11 GB wall -- no OOM. Result line:

    X0^111(1): 1 curve comparison(s), 0 involution comparison(s), 1/1 expected covers matched;
               4 committed model cover(s) re-derived (0 CRV skipped)

It anchors on the FULL genus-7 curve, which is a stronger anchor than `93_1`'s quotient-only one.
⇒ **Both previously un-re-derived Guo-Yang bases now have passing re-derivation tests.**
⚠ Nothing has been copied off lava; the log is `$HOME/x0_111_1.log` there and its clone is still at
`8dac84c`, i.e. behind `main`.

**lovelace: all five still alive** (`34_11` INTSOL at 8 d 2 h, `95_1`, `159_1`, `69_1`, `119_1`), all
elapsed ~ CPU so all progressing. ⚠ **lovelace is USABLE AGAIN** -- load 40 on 256 cores, down from
324. The last handoff's "do not launch there" no longer holds.

### ⚠⚠ THE `A_m` MAIN LINE IS WITHDRAWN — the hatch is blocked on INTEGRALITY, not on a theorem

Full argument and data: `vvdata/weyl-campaign/even-correction/AM-REASSESSMENT.md` on campaign
(`bb3700e`), which also carries a superseding header on that directory's README.

`PLAN.md` named the `A_m` theorem as the main line because it "unblocks 49 obstructed bases and is
the only item that does". **That justification does not survive its own evidence.**

* `A_m` is DEFINED by `sum_m c(-m) A_m = mult(f)` -- that identity is how the values were solved,
  not a property proved of them. So inserting `A_m` into `Kappa0` adds exactly `mult(f) log N` per
  firing CM point.
* `SchoferFormula.m:1024` ALREADY adds precisely that, and it is `prop:kappa0`'s conclusion verbatim.
* That code was LIVE at `619051a`, the commit both hatch branches were cut from -- so the recorded
  17 non-rational cells were measured WITH the correction applied.
* Re-run on current code: **12 of 18 bad cells are at NON-FIRING discriminants** (`-24 -51 -228
  -408`), where no level term is owed and the `A_m` defect is outside its own stated scope; at the
  6 firing cells the outer term DID fire with nonzero `m0mult` and they are non-rational anyway.
* `rem:gauge` already says the `N | m` support rule is a GAUGE, which explains without any new
  theorem both why `N | m` had to be imposed by hand and why `prop:closedcoef`'s `-a_E` "reproduces
  1 of 13".

**⚠ A TRAP THIS CREATES:** implementing `A_m` inside `Kappa0` WITHOUT removing the outer m=0 term
would DOUBLE-COUNT. `SchoferFormula.m:609` says the correction "actually belongs" there -- true as
bookkeeping, but it is a MOVE, not an ADDITION.

**What actually blocks it:** the perturbed form is NON-INTEGRAL. `m0mult = (1/2)c_eta(0)` goes from
integers (baseline) to quarter-integers (perturbed), so `c_eta(0)` is half-integral; a fractional
multiplier puts a fractional exponent on a prime, which is exactly the `RationalNumber` failure --
and unlike the log-`N` story it predicts failures at firing AND non-firing discriminants, which is
what both runs measure. `IntegralSolution := true` does not rescue it (the perturbed divisor admits
no integral form), **and that verdict is safe only because its control was run: the unperturbed
baseline passes cleanly under the same flag** (0 cells, 12 keys).

⚠⚠ **CORRECTED LATER THE SAME DAY — "blocked on integrality" is an OVER-CLAIM.** The sweep meant
to confirm it refuted it. Integrality is real but PARTIAL: baseline 0 cells, the original heuristic
run 18, and **any integral perturbation exactly 11** (164/+2, 164/+4, 56/+4 all give 11), with
`m0mult` integral again. The residual 11 is INVARIANT in discriminant and amount, so it does not
depend on which even divisor is added. The remaining cause is NOT identified, and the next step is
INSTRUMENTATION (print what `RationalNumber` is handed at one stable bad cell), not a fifth
single-cause story — four have now been refuted by controls.
⚠ A genuine defect WAS found in the old instrumentation: its "prefer the largest |disc|" heuristic
sent 4 of 7 keys to disc 296, which is integrally solvable at NO key, while 164 is solvable at all.
⚠ Useful accident: the `180/+4` run applied no perturbation at all (180 is divisible by N=3 and the
probe requires coprimality), giving a NULL CONTROL that returns 0 cells — so the 11s are caused by
the perturbation, not by the harness.

⇒ The hatch is a SEARCH for an even perturbation that is ALSO integral -- two conditions, not
one -- rather than a wait for an open theorem; but both conditions together are still not enough. ⚠ `A_m`/`b` remains a genuine open question in its
own right (no product of local densities reproduces `b`); it is simply not what blocks the 49.

⚠ Not a byte-reproduction of 2026-08-31: 18 cells vs 17, because the CM evaluation set differs
(`-56/-68` then, `-228/-408` now). Same phenomenon, different sample -- do not read the two counts
as a change in the effect.

⚠ Recorded, not chased: **`mult(f)` is NOT determined by `div(f)`.** The INTSOL and default
baselines pick forms differing by a trivial-divisor kernel element and report different `m0mult`
vectors, while both give 0 bad cells and the same 12 keys.

### WHAT THE OBSTRUCTION ACTUALLY IS — and it is NOT the eta quotients failing to span

Asked directly this session, so it is written down here. **The eta-quotient basis is not the
problem, and "the space is one pole order too small" is REFUTED BY MEASUREMENT**
([[borcherds-obstruction-is-real]], probe at `38_5`):

    bump 0:  poleord 190  rows 164  cols 36  rank 35      deficit 1
    bump 8:  poleord 198  rows 172  cols 38  rank 37      deficit 1

Enlarging the weakly holomorphic space adds forms AND divisor-columns at the same rate, so the
deficit is invariant. Decisively, the annihilator `phi` is **stable under enlargement** -- equal
entry-for-entry on shared discriminants and merely extending to the new ones -- which is the
signature of reading coefficients off a FIXED modular form, not of a truncated basis.

**The actual reason is Borcherds' criterion.** A divisor is the divisor of a Borcherds product iff
it pairs to zero against every form of the obstruction space -- the weight-3/2 cusp forms of the
lattice's dual Weil representation. At `38_5` that space is 1-dimensional with generator `phi`, and
`phi(target) = -22 != 0`, so the requested ramification divisor **is not the divisor of any
Borcherds product**. That is why all 96 triples fail identically: the search is futile by
construction, not unlucky.

Two corollaries worth keeping:
* **The working bases have NO obstruction space at all** (`34_3`, `38_7`: rank = cols, deficit 0 at
  every key, target found on triple 1). The obstruction is absent there, not dodged.
* **It is not a size threshold.** `38_7` is strictly LARGER than `38_5` in every dimension and is
  fully surjective. The cokernel dimension is an ARITHMETIC INVARIANT of the discriminant form, not
  a monotone function of `DN` -- which is why `14_19`/`14_29` work while `14_17`/`14_23`/`14_31`
  fail.

⚠ **DO NOT CONFUSE THIS WITH THE ETA-QUOTIENT EXPLOSION**, which is a different failure mode
(TIMEOUT, not form-failure) and was root-caused to a `D0` bug and FIXED
([[odd-d-etaquotient-explosion]]).

⚠ And note the two axes are independent: the deficit persists with a fully integral basis and under
the integral solve. So the obstruction (is the target in the image at all?) and integrality (does
the preimage contain an integral point?) are DIFFERENT questions -- which is exactly why the even
perturbation can fix the first and break the second.

### ✅ `EquationsByRebase` now runs under a pinned `base_label`; both `model_drift_ok` flags are off

`main` `47ea828`. PLAN item B. The `base_label eq 0` gate is gone and the pin is threaded into the
stage's inner `EquationsAbovePointlessConics`, which had been silently reverting to the default base
(propagation adds NEW bases to `re_eqns`, so that was a real hole).

⚠ **`STAR` IS NOT FORCED TO THE PINNED BASE, AND THE FIRST ATTEMPT THAT DID SO WAS WRONG.** Forcing
`STAR := base_label` is what PLAN item B literally specifies; at `26_3` it runs the stage and fills
NOTHING. Measured `base_count` there: `<8092,1> <8098,1> <8103,3> <8104,3> <8105,7>` -- the pinned
base carries 3 first-level equations, the heuristic picks `8105` with 7, and only `8105` admits a
usable Hauptmodul root.

⚠⚠ **AND THE COMMITTED DATA IS A MIXTURE.** `models_26_3.m`'s `[1,2]` and `[1,13]` were filled by
`7a923ae` on a DEFAULT run -- the old gate was `base_label eq 0`, so that run cannot have been
pinned. So the file carries 13 keys in the `base_label := 8103` presentation plus 2 from an unpinned
rebase, and reproducing it REQUIRES the unpinned `STAR`.

**The check that caught it is the one worth keeping**: `model_drift_ok` tolerates MISSING keys, so
"Success!" with the flag on proved nothing. The flagless run is what failed, with exactly
`[1,13], [1,2] NOT PRODUCED AT ALL`. Both flags are now off and both tests pass on their own merits
(`X0_26_3` 189 s, `X0_10_13` 773 s), with both Guo-Yang oracles green -- and that oracle holds
Guo-Yang's own curves for exactly the two keys the rebase fills. Full suite **75/75, 0 failures**.

### ✅ ModelRegen retargeted from a MEASURED sweep — and it immediately found a real drift

All 38 bases with no `X0_*.m` test and no recorded cost were measured, one magma process each
(so an OOM costs one base, not the batch), capped at 600 s. **Only 12 of 38 finished.**
`CHEAP_BASES` 8 -> 14: `+191 s` for `+40` comparisons, against the old list's 802 s for 114.

    added:      34_1 5.8s->4   46_1 13.3s->4   6_5 22.3s->23
                106_1 38.6s->3  122_1 49.6s->3  118_1 61.4s->3
    left out:   34_3 204s->10   178_1 269s->3   202_1 374s->2        (poor value, not failure)
    capped:     the ENTIRE D=6 large-prime-N family (6_23 .. 6_83, 11 bases), plus
                10_17 10_29 10_31 10_37 10_41 10_53 10_61 14_13 14_19 14_29 22_13 34_5 34_7 38_7 58_5

⚠ **`10_3` DRIFTS, and it is REAL AND PRE-EXISTING** -- 3 non-isomorphic entries at `W=[1,2]`,
0 missing. Confirmed NOT caused by this session's change: it drifts identically with
`EquationsCovers.m` reverted to HEAD. The entries are genus-0 conics and two differ from fresh
output by exactly **-1/2, a NON-SQUARE** -- a quadratic twist, i.e. the unpinned-y2-scale class
([[committed-models-can-be-unreproducible]]), not a lost cover. **Which side is correct is NOT
settled: there is no Guo-Yang oracle at `10_3`.** It is deliberately left OUT of the default list
rather than given an undiagnosed `MR_KNOWN_DRIFT` row -- an undiagnosed row is how `26_3`'s went
inert.

⚠ **`15_4` cannot run here at all**: `N = 4` is not squarefree and `BorcherdsForms.m:55` asserts.
A METHOD BOUNDARY, not drift, and it fails in 0.7 s.

### A near-miss worth recording: I nearly reported a clean suite as truncated

`run_tests.m` prints `Tests failed:` **only when `#failed gt 0`**, so its ABSENCE means zero
failures -- it is NOT the truncation signature that `CLAUDE.md` warns about. My file count also used
a pattern that only matches when `Success!` lands on the same line as the filename, which is false
for every test that prints output first, giving 40 instead of 75. Counting the right thing:
**75 expected, 75 started, 75 `Success!`, 0 `Fail!`.** ⇒ When checking for truncation, count
`Success!` occurrences against the suite's OWN file filter, and read the summary's print condition
before treating its absence as evidence.

## Handoff — 2026-09-10

### ⇒ LIVE JOBS AT HANDOFF TIME — collect these before starting anything

**`X0_111_1` is RUNNING ON LAVA** and will outlive this session. To collect it:

    ssh -J lovelace lava
    tail -50 $HOME/x0_111_1.log            # look for "Success!" / "Fail!"
    pgrep -u $USER magma                   # empty => finished (or died)

It was at `m_idx=3 of 7` after 5.6 h. ⚠ Its pool reached **1678 vectors**; the recorded wall is
Magma dying around **~2000 vectors / ~11 GB**, so if the process is gone with no verdict in the log,
suspect OOM rather than a code fault, and re-run with a smaller `Prec` or on a bigger box.
⚠ Its clone is `$HOME/ShimuraCurveALQuotients` on lava at `8dac84c`; **it is now behind `main`** —
`git fetch && git reset --hard origin/main` BEFORE any new run there, but **NOT while that job is
alive** (`AttachSpec` compiles on demand; see [[never-update-a-clone-with-jobs-running]]).

**FIVE jobs on lovelace** (`ps -u $USER -o pid,etime,time,cmd | grep magma`), all ~100% CPU:
`34_11` (INTSOL=1, 5 d 16 h — PLAN's old item 1), `95_1`, `159_1`, `69_1`, `119_1`.
⚠ Do not `git pull` those checkouts while they run. ⚠ lovelace itself is SATURATED by other users
(load 324/256) — launch new work on **lava**, not there.

Everything else from this session is committed and pushed on both branches; the evidence for the
49/49 refresh is at `vvdata/weyl-campaign/obstructed-rerun-2026-09-10/` on campaign.


### ✅ THE OBSTRUCTED CLASS RE-RUN AGAINST CURRENT CODE: 49 of 49, ZERO FLIPS

Every OBSTRUCTED verdict on record was taken **2026-09-01/02**, and `BorcherdsForms.m` has had six
commits since — including **`d9b52d0` (09-05), "shift the oo-side basis by its own valuation, not
the 0-side n0"**, the vx fix, which is a CORRECTNESS fix to the very stage that raises "Failed to
find all Borcherds forms". So the 49-base figure justifying `A_m`'s priority rested on pre-fix
verdicts. Re-run 2026-09-10 with `spanprobe.m` at `PROBE_BUMP=0`:

    49 bases re-run    49 still obstructed    0 flipped    0 failing for another reason
    runtimes 18 s (38_5) to 1349 s (34_19)

⇒ **The obstruction is not an artifact of the pre-vx-fix code**, and `A_m`'s justification is now
refreshed evidence rather than a stale tally. A prediction recorded before the first run ("still
obstructed, ~60/40") held.
⚠ **The 49 was recovered, not assumed**: harvesting every obstructed verdict across `sweep122`, the
triage waves and the span probes yields EXACTLY 49 distinct bases, independently confirming the
"known 28 + 21 new" figure as the union of recorded verdicts.
⚠ **`38_5` returned in 18 s against 901 s recorded** (~6x, from the q-expansion bootstrap), and its
`pole_order=190 pool=164` reproduces the recorded `poleord 190 rows 164` — the same computation, not
merely another failure. ⚠ **Level does NOT predict cost** here either (18 s to 1349 s, uncorrelated
with M) — the third time that lesson recurred in one day.

### ✅ `X0_93_1` PASSES — 13389 s (3.7 h)

`tests/_offline/X0_93_1.m` (new): 1 external comparison against Guo-Yang's typo-corrected `[1,93]`
plus **3 committed model covers re-derived** (1 CRV skipped by design). This base mattered most
because `models_93_1.m` regenerates ONLY since the vx fix, so a silent regression there would have
left every committed artifact looking fine. Pre-flighted before the run, not after: their
`(3s^3-7s^2-3s-1)(3s^3+s^2-3s-9)` is isomorphic to the committed entry and all three refuted typo
repairs still fail.
`tests/_offline/X0_111_1.m` (new) is running on **lava** — it anchors on the FULL CURVE (genus 7,
hyperelliptic, published), which is stronger than 93_1's quotient-only anchor.
⚠ At m 3 of 7 its pool is **1678 vectors**, near the recorded ~2000-vector / ~11 GB wall where
Magma dies. If it disappears, that is the likely cause, not a code fault.

### ⚠ REMOTE MACHINES: lovelace is saturated, and PLAN's "four blockers" is FIVE

`lovelace` load **324 on 256 cores**, dominated by other users (`xw132`'s `k3rank` since Sep 06) —
the memory entry's warning that "idle is not a durable fact" holds. **Do not launch there.**
`lava` (`ssh -J lovelace lava`) was load 0.04 on 32 cores and is where `111_1` runs; it needed its
own clone, and the committed `polymake/` cache came with it.
⚠ **`PLAN.md` says four blockers; there are FIVE Magma jobs**, and the fifth is
**`34_11` with `INTSOL=1`, 5 d 16 h elapsed at ~100% CPU** — PLAN's old item 1, "the best-value
thing here". All five show elapsed ~ CPU, so they are progressing, not wedged.

### The math: two hypotheses formed, two retracted

Both concerned `A_m`; neither survived contact with the sources, and the record is worth more than
the hypotheses were.

1. **RETRACTED: "I derived the level-prime factor at general m."** Both "results" are already in
   `paper/level-prime-kappa.tex` — Result 1 IS `thm:closed` (`W_{m,N}(1) = (N-1) ord_N(m)`, with
   `cor:support` for the `N | m` vanishing, verified there over 180 checks against my 18), and
   Result 2 is in `sec:open`, which carries the same `alpha_k`/`G(X)` recipe AND the counts. Cause:
   I read the memory's "the next theorem is general `m` at a nonzero isotropic coset" as meaning the
   LEVEL PRIME was open at general `m`; it is not — the sentence means the intersection with the
   `D`-part. **Every number was right; I was wrong about which object was already known.**
   ⇒ **READ THE PAPER BEFORE DERIVING.** Memory entries and code are not a substitute for the
   30-page document in the repo.
2. **CHECKED AND DROPPED BEFORE REPORTING: "the `prop:closedcoef` refutation is a wrong-object
   comparison."** `rem:gauge` does say `-a_E` and `A_m` are two representatives disagreeing
   pointwise while both reproducing the multipliers — but (i) the memory POSTDATES `rem:gauge` by
   two days, (ii) its literal claim "`A_m` does not follow from `prop:closedcoef`" is TRUE, and
   (iii) decisively, `SchoferFormula.m:589` specifies the code needs the log-`N` coefficient of
   `Kappa0`, "nonzero exactly when `N | m`" — the LEVEL-supported object, whose support `cor:support`
   governs, not `-a_E`'s embedding support. **The memory is correct; the hatch is genuinely blocked.**

**What survives of the math:** `prop:closedcoef`, transcribed and evaluated against the repo's own
`Hurwitz`, reproduces `rem:gauge`'s stated values EXACTLY (`0,0,1,2,1,2` at `X_0^15(2)`) — a small
reusable confirmation that the closed form and its implementation agree.


    X0_*.m cover comparisons:   126 hand-written + 337 model-derived over ALL 34 bases
                                (was 126, and NOTHING else); 34 of 34 tests pass
    committed cover keys:       863 across 88 model files
      with a re-derivation test:  343 on 37 bases
      with none:                  520 on 51 bases   <- 476 of them validated ONLY by ModelChecks
    X0_* census:                32 of 34 pass; the 2 failures are MISSING-KEYS-ONLY and diagnosed

### The `X0_*` tests re-derived 41% of the covers. They now re-derive all of them, for free.

`test_AllEquationsAboveCoversSingleCurve` compared ONLY the keys hand-written into `cover_data`,
and `if not is_def then continue` dropped the rest IN SILENCE: 128 hand-written cover_data KEYS against
309 populated model keys, nine bases checking 1 of 15. (⚠ KEYS, not comparisons -- a key holds one
entry per base, so the comparison counts above are the larger multiset figures. Different objects;
do not quote one for the other.)

⇒ **The fix was not transcription.** `AllEquationsAboveCovers` is ALREADY PAID FOR by each test,
and `tests/_offline/ModelRegen.m` already had the right comparison -- it was offline only because
it paid for a SECOND pipeline run per base. So that comparison now runs as a second pass inside
the helper, reusing the run it already did:

* `tests/_modelfile.m` (NEW) -- `ReadModelSet(D,N)`. Isolated in its own file because an `eval`
  inside a procedure that closes over an outer variable segfaults Magma 2.29 (the trap that forces
  ModelChecks.m and ModelRegen.m into top-level form). PROBED in isolation before being built on.
* `tests/BorcherdsProducts.m` -- ModelRegen's MULTISET matching (so a 3 -> 2 loss cannot hide
  behind two committed entries matching one survivor), `<genus, f, h>` handling, CRV entries
  skipped AND COUNTED, a zero-comparison guard, and `model_drift_ok`.

**MEASURED at 6_11: 1 comparison -> 1 + 17, in 119.8 s against a 121.5 s baseline.** The check is
free; the pipeline run was the cost all along.
**NEGATIVE-CONTROLLED:** perturbing one entry to a same-genus DIFFERENT curve and adding a key the
AL group cannot produce makes it fail, naming both causes separately. It could have failed.

⚠ **IT IS A DRIFT CHECK, NOT A VALIDATION.** It says "current code still produces this", not "this
is correct" -- the committed file is what the pipeline itself wrote. Correctness still comes from
ModelChecks (Eichler-Selberg point counts) and the Guo-Yang oracles. The hand-written `cover_data`
entries must NOT be deleted in favour of it: those are Guo-Yang's PUBLISHED equations, and they
are the only entries carrying labelled involutions.

### ⚠ A pinned `base_label` loses EXACTLY the keys `EquationsByRebase` filled

Two tests fail, both MISSING-KEYS-ONLY, zero non-isomorphic anywhere in 34 bases:
`10_13` (`[1,2] [1,5] [1,26]`) and `26_3` (`[1,2] [1,13]`).

`AllEquationsAboveCovers` gates `EquationsByRebase` on `base_label eq 0` (`EquationsCovers.m:1061`),
so a test pinning a non-zero `base_label` cannot reproduce a key the rebase FILLED on a default run.
Both model files say so in their own headers -- `models_26_3.m` even names `[1,2]` and `[1,13]` as
the two that were empty and were "filled, unlocked by EquationsByRebase".

⚠ **THE CONTROL GROUP is what makes this a diagnosis and not an excuse.** `14_3`, `21_2` and `6_17`
also pin a `base_label` and ALL THREE PASS: `14_3`'s empties were fixed by the COPRIME FILTER FLIP,
not the rebase, and the other two never had any. The gate costs the rebase-filled keys and nothing
else.
⚠ **A PREDICTION WRITTEN DOWN BEFORE THE RUN WAS HALF WRONG, AND THAT IS WHY THE RULE IS NOW EXACT.**
It predicted drift at `10_13 14_3 21_2 6_17` and a PASS at `26_3`; the opposite happened for four of
the five. Had the flag been set from the prediction, three tests would have been needlessly
weakened and `26_3`'s real cause never found.

⇒ `model_drift_ok` therefore tolerates **MISSING keys only**. A key the pipeline DOES produce must
still be the committed curve, whatever `base_label` was pinned -- silencing both with one flag would
hide the failure that actually matters.

⇒ **OPEN, and well-evidenced: relax the `base_label eq 0` gate.** `EquationsByRebase` only ever
fills keys that are ALREADY EMPTY, so running it under a pinned `base_label` should not disturb the
pinned presentation's other covers -- and it would make both these tests reproduce their full model.
A pipeline change, so it needs oracle validation, not just a green test.

### The re-derivation gap is bigger than "128 of 309" -- that counted only the tested bases

    88  model files, 863 cover keys
    37  bases have a re-derivation test (CI or offline)  ->  343 keys
    51  bases have NONE                                 ->  520 keys (60%)
    44  bases are validated ONLY by ModelChecks          ->  476 keys (55%)

"Only ModelChecks" is not nothing -- genus, Weil divisibility and Eichler-Selberg point counts,
none of which touch the Borcherds machinery. But it NEVER RUNS THE PIPELINE, so drift there was
invisible to everything in the repo.

⚠ **And ModelRegen's default `CHEAP_BASES` had become PURE DUPLICATION: all nine had an `X0_*`
test.** Retargeted at bases with none. MEASURED PER BASE, because a batch total cannot tell a
3-minute base from a 26-minute one:

    6_1 10_1 14_1 22_1 6_7 6_13   72 comparisons, ~5 min for all six
    10_7            15 keys       26 comparisons, 185 s
    26_5  804 s | 14_11 1475 s | 22_7 1591 s | 65_1 813 s     <- measured, LEFT OUT for cost

⚠ **KEY COUNT DOES NOT PREDICT COST**: `10_7` has 15 keys and costs 185 s; `65_1` has 4 and costs
813 s. `14_43` was killed at 7 h 44 m unfinished.

⚠ **TWO CLAIMS I FIRST WROTE HERE WERE WRONG, both caught by being challenged rather than by a test.**
1. *"The new list is all even `D`, which is a hole."* **It is not a hole**, and "both D parities" from
   the old comment is itself the stale part. **10 of the 14 odd-`D` model bases have an `X0_*` test**
   (`15_1 15_2 21_2 35_1 39_1 51_1 55_1 57_1` in CI, `39_2 87_1` offline) and every such test now
   re-derives EVERY cover key, so odd-`D` model building is well exercised without ModelRegen. And
   **there is no D-parity branch in the code ModelRegen drives**: the only live `IsEven(D)` uses are
   in AL fixed-point code (`ShimuraQuotients.m:842`, `GeneralizedComplicatedFixedPoints.m:125,186`)
   reached from the FILTER/triage pipeline, never from `AllEquationsAboveCovers`;
   `BorcherdsForms.m:9`'s `assert IsEven(D)` is commented out. Parity mattered when ModelRegen was
   the only re-derivation for those bases; it is not any more.
2. *"`65_1` is the only odd `D` among the 51."* It is the only odd `D` among the **44** with neither
   a test nor an oracle mention. Among the **51** without a re-derivation test there are **four**:
   `111_1`, `15_4`, `65_1`, `93_1`. I quoted a figure for one set while naming the other.

⇒ Both were inherited framing rather than measured claims — the first copied from the comment being
replaced, the second a set I had computed earlier for a different purpose. Spend the ModelRegen
budget on COST, not parity.

**RUN END TO END with the new default: 8 of 8 reproduce, 0 drifted, 114 comparisons, 802 s** — inside
the ~20 min the old list cost, and none of the 7 new bases is re-derived anywhere else.

⚠ **AND THAT RUN FALSIFIED MY OWN REASON FOR ONE ENTRY.** `26_3` was included "to keep the
known-drift path exercised"; it reports `OK (16 compared, 1 CRV skipped)` — it REPRODUCES. The drift
`MR_KNOWN_DRIFT` records for `26_3` is entirely in its `W={1}` entry, that entry is a `"CRV"` entry,
and ModelRegen SKIPS every CRV entry. So **the `26_3` row of `MR_KNOWN_DRIFT` is INERT** — it cannot
fire under this code and has probably been inert since CRV skipping was added. With `14_43` out of
the default list, **the known-drift tolerance is now exercised by nothing**. A tolerance that cannot
fire is the mirror image of a check that cannot fail, and it was only caught by running the default
list instead of trusting the reasoning that chose it.

### The oracle's genus-0 branch was a one-bit check

`GuoYangQuotientOracle.m` compared genus-0 quotients by `HasRationalPoint` alone, so ANY two
POINTLESS conics MATCHED -- and 57 of its 170 comparisons (34%) take that branch. Now a real
`IsIsomorphic` on the conics, negative-controlled on `y^2 = -x^2-1` vs `y^2 = -x^2-3` (both
pointless, correctly distinguished).
⚠ **BE HONEST: it changed no verdict.** Every genus-0 quotient at all 20 oracle bases is a POINTED
conic, and pointed conics over Q are all isomorphic to P^1, so the old check was ACCIDENTALLY
equivalent; all 57 still match in the same 4.8 s. It is not equivalent in general -- 73 of the 281
genus-0 entries in `data/models/` ARE pointless (`6_5 6_7 6_83 82_1 93_1`) -- so this guards the
first such oracle base rather than discovering anything.

### ⚠ A LATENT SILENT-CORRUPTION RACE IN CONCURRENT RUNS (found, checked, NOT yet fixed)

`BorcherdsForms.m:180` writes every Normaliz solution to a DETERMINISTIC SHARED PATH
`polymake/polymake_solution_<M>_<n>_<m>`, and reads it back as `FileExists` -> `eval Read`.
`nmzsolve.py` writes that file NON-ATOMICALLY. So two concurrent Magma processes needing the same
UNCACHED triple can have one read a PARTIALLY WRITTEN file -- a valid-looking but truncated point
list, i.e. exactly the "a partially-cached base returns a wrong answer rather than an error" mode
`CLAUDE.md` flags as critical.

⚠ **This session ran up to five Magma processes at once, so it was exposed.** Checked rather than
assumed: NO `polymake_solution_*` was written during any of it and `polymake/nmzsolve.err` does not
exist, so every solve hit the committed cache and no race occurred. The results stand.
✅ **FIXED**: `nmzsolve.py` now writes via a pid-suffixed temp file and `os.replace` (atomic on
POSIX), at both write sites. Validated: byte-IDENTICAL output to the old writer on the same point
list, no temp file left behind, and exercised END TO END by a real `nmzsolve.py` invocation (not
just the cache-read path, which is all a passing test would have touched).
⚠ SHARED-PATH FILE -- **merge it down to the campaign branch.**

### ⚠⚠ AND A SECOND, WORSE ONE FOUND WHILE TESTING THAT: THE SOLUTION CACHE KEY IS INCOMPLETE

The cache key is `(M, n, m)` ONLY. It omits `k`, `sq_disc` and `cuspidal` -- **and all three change
the answer.** Measured at `(M,n,m) = (8,1,0)`, varying only the omitted parameters:

    k24=12 sq_disc=1 cuspidal=0  ->   4 points
    k24=12 sq_disc=0 cuspidal=0  ->   4 DIFFERENT points
    k24=24 sq_disc=1 cuspidal=0  ->  10 points
    k24=12 sq_disc=1 cuspidal=1  ->   0 points

So two call paths asking for the same `(M,n,m)` with different parameters means the second silently
gets the FIRST one's point set -- correct arithmetic about the wrong object, no error anywhere. And
the two call sites DO differ: `HolomorphicEtaQuotients` (`BorcherdsForms.m:194`, live, reached from
line 289) passes `sq_disc := true` pinned at `(M,0,0)`, while the Borcherds path (line 425) takes
the `sq_disc := false` default.

⚠ **LATENT, NOT MATERIALISED** -- measured, not hoped: of the 503 committed solution files **NONE is
a `*_0_0`**, so the `(M,0,0)` site has never cached anything and nothing can be mis-served today.
⚠ **DO NOT WIDEN THE KEY WITHOUT MIGRATING THE CACHE IN THE SAME COMMIT.** All 503 files are named
under the narrow key; widening makes them all invisible, and above the cached frontier a fresh solve
fails SILENTLY. That would turn a latent collision into a guaranteed silent regression everywhere.
Documented at the read site in `BorcherdsForms.m`.

⇒ **HOW IT WAS FOUND, because the method generalises:** regenerating a committed cache file to check
the atomic-write change gave a DIFFERENT point set. The tempting read was "my change broke it". The
actual cause was that the filename does not record the parameters, so I could not reconstruct the
original constraint system -- and that *is* the bug. A mismatch was evidence about the CACHE KEY,
not about the edit, and the writer had already been proven byte-identical independently.

### ⚠ THE BRANCH-DIVERGENCE INVARIANT IS RED, and not in the harmless direction

`CLAUDE.md`'s check -- `git diff origin/main origin/m0-theta-campaign --name-only --
':!vvdata/weyl-campaign/*'` -- "should print nothing but doc files". It prints **53**, including
`EquationsCovers.m`, `SchoferFormula.m`, `run_tests.m`, 15 model files and 30+ test files.

Direction checked, not assumed: **main-only 51 commits, campaign-only 180, and campaign is NOT an
ancestor of main.** So campaign carries real independent work (`rankcheck_gauge.py` on the
`rem:gauge` ambiguity, `cusp7.m`) AND is missing all 51 of main's recent commits -- which include
`EquationsByRebase`, the quotient oracle and the model fills.
✅ **RESOLVED the same day: `main` merged into `m0-theta-campaign`, no conflicts, both pushed.**
The invariant now prints **NOTHING AT ALL** (not even doc files), and main-only commits are **0** —
campaign contains everything on main. Sanity-checked by running from the campaign worktree itself,
which is the only thing that proves the point: `X0_38_1` passes in 8.5 s and the quotient oracle
makes its 170 comparisons there. Campaign keeps its own 182 commits of research work.
⚠ It will drift again the moment `main` moves. **Run the invariant, do not rely on discipline.**
⚠ NOT affected: `tools/regen-model.sh` runs campaign's `genmodels.m` but from the main checkout's
cwd, so `AttachSpec` loads MAIN's packages. Model regeneration is fine.

### A grep that read a fragment and generalised (again)

Building "which bases have an external oracle" by matching `<D, N,` tuples MISSED `93_1`, whose
Guo-Yang check is a bespoke `gy93_*` block at the END of `GuoYangEquations.m`. The count was rebuilt
searching for the model FILENAME and the `D_N` tag too. Caught only because the number contradicted
what the project already knew -- the same failure mode `CLAUDE.md` opens with.

## Handoff — 2026-09-09

    Guo-Yang published equations:   42 reproducible bases
    full curve stored:              38      <- 10_19 and 22_5 regenerated today
    remaining blockers:              4      95_1  119_1  159_1  69_1  -- still running on lovelace
    X0_*.m tests checking involutions:  34 of 34   (was 23 on 09-07)
    Guo-Yang quotient comparisons:     188 over 24 bases, 0 skipped, 0 mismatches

### The big change: the quotient oracle

⚠ **`GuoYangEquations.m` compares only the equations Guo-Yang PRINT** — usually the full curve
alone — so a base with fifteen cover keys got ONE external comparison. But they also print the
INVOLUTIONS, and every quotient follows from those:
`CurveQuotient(AutomorphismGroup(C,[w]))` is `X/W`. That turns one comparison per base into one per
cover key. `tests/GuoYangQuotientOracle.m` (generic, 20 bases) plus four hand-derived files for the
CRV bases now make **188 comparisons in a few seconds**.

**Three errors in Guo-Yang's tables are now determined**, each by evidence rather than preference:
* `93_1`: `-3t` is a typo for `-3s` (confirmed by the journal version).
* `14_5`: the table's `w_35` sign is wrong; **their own Example 36** has it right.
* `10_13`: the table **SWAPS `w_10` and `w_13`** — settled by Ogg's fixed-point rule, with the
  clincher internal to their paper: **their own CM table** puts disc `-52` at Hauptmodul `0` and
  `-40` at infinity, contradicting their involution table and agreeing with our pipeline.

### The pipeline change: `EquationsByRebase`

An EMPTY cover key is often a **Hauptmodul normalisation artefact, not an obstruction**. The
pipeline builds a genus-`g` curve only as a FIBRE PRODUCT, needing degree exactly `g+1` over a
shared base; which degree a quotient has depends on whether infinity is a branch point, which is
ours to choose. `t -> r + 1/u` at a RATIONAL ROOT fixes the degree profile. At `22_5` this
reproduces Guo-Yang's degree-12 polynomial VERBATIM.

Wired in as the last stage of `AllEquationsAboveCovers`; it is a **no-op unless some cover is
empty**, and only ADOPTS keys that were empty. The `ws` transport works because the rebase is
LINEAR on the weighted ambient: `psi = (r*x + z, y, x)`.

⚠ **Cost**: bases WITH empty covers get slower (`X0_10_11` 471 s, `X0_10_13` 872 s). Bases without
pay nothing.

### Traps that cost real time today, all now guarded

* **`run_tests.m` was SILENTLY TRUNCATING the suite.** It globbed every `tests/*.m`, including
  helpers; several end with `exit;`, which kills Magma. 72 of 79 files reported and no summary was
  printed — a truncated run looks like a clean one. Now mirrors the CI matrix. **Check the file
  count and that `Tests failed:` is present.**
* **A model entry may be `<genus, f, h>`, meaning `y^2 + h*y = f`.** Dropping `h` gives a DIFFERENT
  curve of the same genus; 9 entries across 7 files have one. This produced a false "defect" report
  against `models_87_1`, retracted.
* **An incomplete oracle cannot refute anything.** `GuoYangQuotients_10_19.m` was missing `w_10`
  and `w_95`, so "matches no quotient" really meant "matches none of the ones I computed" — that
  produced a wrong retraction of the rebase lever, since re-corrected.
* ⇒ **All three of the day's wrong verdicts were REFUTATIONS**, each correct about its arithmetic
  and wrong about its object. **A failing check needs its object verified as much as a passing one.**
* **A SKIP is a silent gap**: 13 oracle comparisons were being skipped behind a green summary line.
  All recovered; the cause was on the Guo-Yang side (`CurveQuotient` returns a plain `Crv`).
* **`tools/regen-model.sh`'s flag table had gone stale on both rows** — `CMNONCOPRIME` is a dead
  name (the code reads `CMCOPRIME`), and `Y2TWIST` would have produced models differing from the
  committed files and been read as drift. Table now empty; `51_1` and `22_5` verify IDENTICAL.

### Filed upstream

**[Magma-Maths/Magma#123](https://github.com/Magma-Maths/Magma/issues/123)**: `AutomorphismGroup`/
`CurveQuotient` fail for curves in weighted projective (toric) ambients — `IdentityMap` returns a
`TorMap`, not a `MapAutSch`. It blocks the oracle on exactly the `CRV` paired presentations, which
is why `10_19`, `22_5`, `10_13` and `26_3` each need a hand-derived oracle file.

### Later the same day — the fill, and an audit that found a systemic gap

* **ALL 18 remaining EMPTY cover keys filled**, across `6_29 6_31 6_37 10_11 10_13 10_23 14_5
  26_3`. **0 empty cover keys remain: 347 of 347 populated across 38 Guo-Yang bases.** Every filled
  key was checked against the quotient oracle **in a scratch directory BEFORE installing** — that
  ordering is what makes the data trustworthy, and it should not be inverted.
* ⚠ **`14_5` gained two cover keys that never existed in the file** (`[1,5,7,35]`, `[1,7,10,70]`).
  Its AL group has order 8, so there are 15 proper cover keys; the file had 13.
* ⚠ **Existing entries can come back RESCALED BY A SQUARE** (11 did at `10_19`). That is a
  re-presentation, not a regression. Verify entry-by-entry isomorphism; only a MISSING cover is a
  failure.
* **`10_13`'s labelling differs from Guo-Yang by a GROUP AUTOMORPHISM**, and only half is proven.
  Ours differs by `5 <-> 26` AND `10 <-> 13`, fixing `2, 65, 130`; the map is multiplicative so the
  swaps stand or fall together. `10 <-> 13` is PROVEN by fixed points with their own CM table as
  clincher. `5 <-> 26` CANNOT be: both quotients are genus 2 and Riemann-Hurwitz forces `r = 0`, so
  both involutions are FIXED-POINT FREE and Ogg's rule says nothing. We adopt our labelling for
  both; the second half is **inferred by consistency, not established**.
* ⚠ **A SECOND SILENT TRUNCATION, pre-existing**: `tests/test_weil_polynomial.m` ended with
  `quit;`, which kills Magma since `run_tests.m` evals every test in one process. It sorts
  second-to-last, so `trace_formula.m` never ran locally — which is probably why it was believed
  deliberately skipped. It is not slow: **2.4 s**. Fixed.
* ✅ **`tests/_offline/X0_87_1.m`'s long-standing failure DIAGNOSED AND FIXED**, and validated:
  **passes in 4081 s**. The cause was a **DROPPED h-TERM** — the model stores `[1,29]` as
  `<3, f, h>` with `h = x^3+x^2+1`, and the generator emitted only `f`, so the test compared a
  DIFFERENT curve of the SAME GENUS. It stays offline because it is slow, not broken.
* **`X0_206_1` went 1 -> 4 of 4 covers**, including its `h`-bearing `[1,103]`.

⚠⚠ **AND THE AUDIT THAT MATTERS MOST: the `X0_*` tests re-derive only 41% of the covers.**
**128 `cover_data` keys against 309 populated model keys.** Nine bases check 1 of 15
(`10_11 10_13 10_23 6_11 6_17 6_19 6_29 6_31 6_37`) and eleven check 1 of 4. The helper SILENTLY
SKIPS an absent key, so this is invisible unless counted.
⇒ The MODELS are well checked (~190 oracle comparisons over 25 bases against Guo-Yang). What is
thin is the **RE-DERIVATION** claim — that the pipeline reproduces them — which for most bases
rests on ONE cover. Closing it is mechanical but must handle `<genus, f, h>` entries and `CRV`
pairs, both of which have already caused defects, and it costs CI time.

### Still open

* **The 41% re-derivation gap above** is now the largest single opportunity: work the thin tests in
  order of missing covers, starting with the nine at 1-of-15.
* The four lovelace blockers are mid-FIRST-PHASE after ~3 days; weeks away, not days.
* `93_1` and `111_1` still have no `X0_*` re-derivation test (14-20 h per run).

## Older — Handoff 2026-09-07

**Supersedes** the 2026-07-17 handoff about producing cover models, archived as
`HANDOFF_2026-07-17.md`. That task is not dead, but it is gated on the blocker described below.

Everything here is committed and pushed. **`git pull` first — local `main` may be stale.**

**➡ For what to do next, see `PLAN.md`** — five tracks, a do-not list, and the recurring traps.
This file is the record of *what happened*; `PLAN.md` is the record of *what to do*. When the two
disagree about state, this file wins.

## Handoff — 2026-09-07 (the 09-06 section below is still accurate, just earlier)

    Guo-Yang published equations:  42 reproducible bases
    we now have a model for:       38      <- 111_1 recovered
    remaining blockers:             4      95_1  119_1  159_1  69_1   -- ALL RUNNING

* **`111_1` recovered** — 20.2 h, DEFAULT flags, another base the vx fix unblocked. Verified by
  **exact full-curve `IsIsomorphic`, true in 0.05 s**, which is cheap only because its `W={1}` is
  HYPERELLIPTIC. ⚠ Pinnable to one commit (`f87b0ae`, clean tree) — unlike `10_61`/`14_43`.
* **`93_1` and `26_3` upgraded to FULL-CURVE PROOFS** (were quotient-level). `IsIsomorphic` hangs on
  CRV pairs, so the isomorphism is CONSTRUCTED: Mobius map from the hyperelliptic quotient, both
  sides carried by constant squares, then `IsIsomorphism` certifies it. Hundredths of a second.
  ⇒ **The BASE chooses the `V_4`** (assaferan): `26_3` would not match until rebuilt with
  `base_label := 8103`, which is the `V_4` Guo-Yang use. When a CRV pair will not match, try
  another base before concluding anything about the curve.

### Three blind spots removed, each of which immediately exposed a real defect

* **`VerifyModelSet` skips every `CRV` entry** — so 21 paired presentations across 16 files had
  NEVER been checked. `tests/CRVStructure.m` found **5 storing their parent conic twice**
  (reducible schemes, not the genus-1 curves recorded). ROOT CAUSE: at `g = 1` the required degree
  `g+1 = 2` is also a conic's degree, so the conic could fill BOTH roles in the fibre product.
  Fixed; those covers now defer.
* **The `X0_*` helper silently skipped unmatched cover keys** — it could pass while verifying
  NOTHING. It now counts comparisons and errors on zero. That immediately turned CI red, correctly:
  **`X0_10_19` had been passing green in CI while making ZERO comparisons**, at 84 min a run.
* **CI never set `NORMALIZ_BIN`** — so polytope solves failed SILENTLY ("no solutions", not an
  error). Now installs `normaliz-bin` and exports the path. ⚠ Scope was MEASURED: every other
  `X0_*` job reported full coverage, so `10_19` was the only affected test.

### The coprime guard: no evidence it is needed

Full sweep of the 11 `N>1` `X0_*` tests (for `N=1` the filter is provably a no-op, which excludes
19 of 30 rigorously). **8 of 10 pass identically with `CMNONCOPRIME` on and off.** The 2 failures
(`10_13`, `6_17`) are both CRV tests whose PINNED COORDINATE MATRIX breaks under re-presentation —
not correctness. Removing that artifact is what `tests/_crviso.m` does.
⚠ A first sweep appeared to show `10_13` failing; that was MY OWN foreground timeout killing the
sweep's Magma, which then recorded a killed run as a failure. Retracted.

### Process notes

* **`nohup ... &` inside a background call reports completion for the WRAPPER**, not the job.
* **My own foreground timeout killed a background sweep** — a killed run and a failing run are
  indistinguishable in a one-line summary. Capture the error text before believing a regression.
* Bugs of mine caught only because a test could fail: an eager `AutomorphismGroup` (5x slowdown,
  870 s -> 71 min), an inverted conic scalar (`x/rg` for `rg*x`), hardcoded variable order, and
  image polynomials built in the wrong ring. Each was invisible in the first case tried.

## Handoff — 2026-09-07, later (test coverage; supersedes the earlier 09-07 block on these points)

**Re-derivation coverage went 4 → 11 Guo-Yang bases.** Passing an `X0_D_N.m` test IS reproduction,
the stronger claim than `GuoYangEquations.m`'s stored-model comparison. Now covered:
`51_1 55_1 57_1 14_5 14_3 26_3 21_2 15_2 22_3 22_5` in CI, `39_2` offline. 34 `X0_*` tests in CI.
* `X0_21_2` is the **first test that checks a CRV entry** — possible only because the helper now
  CONSTRUCTS those isomorphisms (`tests/_crviso.m`) instead of calling `IsIsomorphic`, which hangs.
* `MR_KNOWN_DRIFT` is down from 5 to **2**: only `14_43` (`INTSOL=1`) and `26_3` (deliberate
  `base_label := 8103`).

**⚠ `Y2TWIST` WAS THE WRONG SUSPECT, and I nearly flipped it on a confounded measurement.**
`PROVENANCE.md` had predicted for two days that making twist selection default was "the right
long-term fix" for `15_2`/`22_3`/`22_5`. Measuring `Y2TWIST=1` against the COMMITTED models showed
large gains — and those were the **coprime flip's**, from hours earlier the same day, because every
committed model predated it. The control is flag-on vs flag-off on the SAME code: run that way all
three bases are IDENTICAL and the deferral path logs zero messages. The selector never fires. The
flip was reverted; the mechanism is kept (it is sound: unique-or-defer) but not defaulted.
⇒ **Compare against a current baseline, never a committed artifact.**

**The real win was already sitting there.** Those three needed NO flag — their committed files were
simply STALE. Regenerated with the plain recipe: `22_5` 3 → 11 populated covers, `15_2` 12 → 15,
`22_3` 13 → 15, nothing lost, `GuoYangEquations` still passing.

**`X0_87_1` is the one known-broken test** and is under diagnosis. Established: the MODEL is fine
(`ModelRegen` reproduces it; `GuoYangEquations` matches its `W={1}`), and the test is well-formed
(expects exactly the model's 4 single-entry keys). So the failure is `assert is_isom`. Leading
hypothesis, which has bitten twice already: `ModelRegen` compares the AGGREGATED model while the
helper iterates EVERY BASE of every cover, so a second base with a different presentation fails only
the helper — fixed at `26_3` and `21_2` with a `base_label`.

## Handoff — 2026-09-06 (this session; supersedes the state notes below)

**Five models produced, and the Guo-Yang denominator was wrong.**

    Guo-Yang published equations:  42 reproducible bases (NOT 43 -- see below)
    we now have a model for:       37
    remaining blockers:             5   95_1  111_1  119_1  159_1  69_1   -- ALL RUNNING

* **`93_1`** — the vx fix unblocked it (default recipe, 14.1 h). It also **settles a typo in their
  table**: their `-3t` is `-3s`, determined by isomorphism from our own model against three refuted
  alternatives, and later **confirmed independently** by the journal version.
* **`26_3`** — recovered with `CMNONCOPRIME=1`, 189 s. Full `V_4` diagram matches; their conic
  `-8x^2-3` comes out coefficient for coefficient.
* **`15_4`** — a FOURTH provenance category: **literature-derived, not pipeline-produced**, and it
  never can be (see below). `a = -1` is confirmed by our point counts, `b = -1` by the full-curve
  trace-formula comparison.
* **`10_61`, `14_43`** — first two models out of the OBSTRUCTED class (41 h, 42 h). ⚠ Neither is a
  Guo-Yang base, so **no external oracle** — `ModelChecks` alone. Weaker evidence; quote it as such.

**⚠ 42, NOT 43.** `15_4` is outside the Guo-Yang method *by the authors' own statement* — their
published Remark 39 says the normalizer of the Eichler order strictly contains the Atkin-Lehner
group there, so the star quotient our pipeline forms is the wrong object. It is not a blocker; it
is out of scope.

**⚠⚠ WE HAD BEEN READING THE SUPERSEDED PAPER.** arXiv:1510.06193 has exactly ONE version (2015).
The paper of record is **Compositio Math. 153 (2017) 1-40**, substantially revised and NOT on
arXiv; our `ShimuraCurves-arxiv.tex` is the arXiv one. The journal fixes `93_1`'s equation and
`39_2`'s involutions, and adds Remarks 38 and 39. PDF is in the user's Dropbox. **Check the
journal, not just v1.** Tu (Pacific J. Math. 269 (2014) — also now in that Dropbox folder, free
from MSP, not on arXiv) confirms `15_4` and covers `26_3`, but supplies nothing for any other
non-squarefree base.

**Speedup shipped:** the q-expansion bootstrap, `qexp(t^j f) = qexp(t)^j qexp(f)`, both sides —
**up to 18.6x** on that step at `pole_order 800`. ⚠ NOT yet shown to help any base end to end.

**The coprime guard is under doubt.** Three bases (`39_2`, `14_3`, `26_3`) produce Guo-Yang-matching
models with it OFF, and `26_3` is the very base whose bad discriminants justify its existence. A
targeted sweep is running. ⚠ Only the 11 `N>1` `X0_*` tests can possibly show anything: for `N=1`,
`gcd(d,1)=1` makes the filter provably a no-op.

### Process lessons this session cost something to learn

* **Never `git pull` a clone that has jobs running from it.** I did, to lovelace, with eight jobs
  running from that directory. `AttachSpec` loads packages ON DEMAND, so a long run can compile
  source that changed under it. `10_61` and `14_43` cannot be pinned to a single commit because of
  it. Launch long runs from a COPIED tree.
* **My own foreground timeout killed a background sweep**, and the sweep recorded the killed run as
  a FAILURE. That produced a false "`X0_10_13` breaks under `CMNONCOPRIME=1`", since retracted. A
  killed run and a failing run are indistinguishable in a one-line summary.
* **`nohup ... &` inside a background call reports "completed" for the WRAPPER**, not the job. Check
  process state; do not trust the notification.
* **Corrections made:** `93_1` was first reported as "34 -> 35" (that is the CM-TABLE count, a
  different set); `EchelonForm` was called negligible when it has the STEEPEST growth (~`PO^4.4`);
  `15_4` was diagnosed as a squarefree-`N` code issue when the authors had stated the real reason
  in a version we had not read.

## ⇒ READ THIS FIRST — 2026-09-04, late

> ### Spend your effort on WHICH OBJECT the claim is about, not on whether the computation is right.

That is the single most useful thing this session produced, and it was learned the hard way. Nearly
every error made on 2026-09-04 had the same shape: **the arithmetic was correct and the object was
wrong.** A rank computed over MONOMIALS when the claim was about FORMS (twice, in two different
ways). A LaTeX parser emitting perfectly valid polynomials from silently truncated input. Three
`grep`s that read a fragment of a file and generalised from it. Every validation in place was of
the form "is this number computed correctly" — **none of them could catch "is this the right
number."**

Two habits did catch things, and are worth keeping:
* **Reproduce a KNOWN value before trusting a new one.** The rank result was only believable
  because the same script had to reproduce the paper's own rank-4 panel and a principal part known
  from `tests/M0Multiplier.m`. Both caught silent data corruption that had produced a
  plausible-looking right answer for the wrong reason.
* **Draft an edit instead of applying it.** The claim "the paper is wrong, fix `rem:gauge`" was
  retracted *while writing the diff*, because writing it forced a close enough read of `sec:exact`
  to notice it says FORMS where I had MONOMIALS. Applying directly would have degraded a correct
  argument in a paper heading for submission.

### State right now

**Three runs were left going on lovelace** (`~/shimura/models/*.genmodels.log`, `M0PROGRESS=1`).
**All three completed their a0 tables**, which is the result they were launched for:

    34_11   13 fallback points of 64      (~8.7 h elapsed)
    10_61   27 fallback points of 64      (~6.7 h)   <- previously DIED at gate 3
    14_43   22 fallback points of 64      (~6.7 h)   <- previously DIED at gate 3

**A completed a0 table means the two-point check never fired, i.e. GATE 3 IS CLEARED at `10_61`
and `14_43`** — the gate that killed both before. That retires the "924 gate-3 failures at 10_61,
zero near-misses" entry and, with it, the claim that `10_61` "is not runnable, it has a real
upstream defect". The fallback rates scale as expected: 1/64 at `15_2` (M=60), 8/64 at `58_5`
(M=580), 22-27/64 at M≈1200 — the threshold is M-scaled, so bigger bases trigger it more.
All three were still in FINAL ASSEMBLY (class-constancy, isotropic agreement, rational snap) when
this was written; **check those logs first — they may have finished, and `10_61`/`14_43` producing
models would change the sweep record's "two runnable candidates of 122" count.** lovelace is shared and busy again (~68 magma
processes, mostly other users) — check `uptime` before adding load. Both branches and lovelace's
clone are clean and in sync; housekeeping list is empty.

### What shipped

* **The per-coset `tau` fix** (`475e72b`) — the only thing that moved the mathematics. `34_11` went
  from failing to passing; validated three ways (`15_2` exact, `58_5` keeps its models, and both
  match the Prop 9.15 closed form 9/9 — `58_5`'s as a PRE-REGISTERED prediction). New
  `M0PROGRESS=1` diagnostic (WriteStderr, because buffered `printf` is lost when a run is killed).
* **`51_1` and `57_1`** — never blocked, just never run. Guo-Yang coverage 32 → **34 of 43**;
  81 model files, `ModelChecks` 8309 checks 0 failures.
* **34 Guo-Yang CM-value tables** as offline tests (`tests/_offline/`, 254 checks, 0 failures) and
  **`tests/GuoYangEquations.m`** comparing 7 committed models to the PUBLISHED equations — the
  first thing in CI that checks our output against the literature. Both negative-controlled.
* **Housekeeping**: branches 4 → 2, 10 archive tags, lovelace pulled up (was 52 behind), the
  t-shift fallback ported to `main` after 9 days missing, 109 stranded Normaliz solves harvested.

### ⚠ What was CORRECTED — do not re-litigate these

* **`rem:gauge` is CORRECT. Do not edit the paper.** It was claimed wrong twice, from measuring the
  wrong object. The `oo`-only model is valid for genuine Borcherds forms (39 of them) and FAILS on
  monomials (residual 2.08), so a rank over monomials says nothing about it.
* **`A_m` genuinely needs new mathematics.** `sec:determined` determines the CONSTANT TERM at
  isotropic cosets ("an indicator, not a phase"); the canonical representative is the SCALAR
  `-a_E`. All-`m` at a nonzero isotropic coset is absent. No shortcut by extraction.
* **`93_1`/`95_1`/`159_1` are the vx class**, not squarefree-`N` — all have `N=1`, which IS
  squarefree. Measured: Magma's own `GalFldFun.m:305 assert vx ge 0`. `genmodels.m`'s `vx_skip`
  list is INCOMPLETE.
* **`gtsweep`'s `FIRE` lever does nothing** (measured on all three bases it claimed to fix).
* **The `cusp7` "scoped implementation task" was falsified** — a monomial-pool coverage gap, not a
  dump bug. MAIN LINE's "provably resolvable" premise is dead with it.
* **Gate 4's "GENUINE 43% violation" was numerical**, not mathematics.

### Where to pick up

`PLAN.md` "Picking this up cold" is current. In short: wait on the three runs; then the 9 remaining
Guo-Yang blockers are correctly classified for the first time (1 structural, 3 vx, 1 nonintegral,
1 non-rational, 2 odd-`D` basis ceiling, 1 open anomaly — `26_3`'s exact `z -> z/(z-1)`
involution). `22_5` and `14_3` need full-curve models GENERATED, not transcribed.

## Update — 2026-09-05: Guo-Yang coverage re-measured, and two code changes

Six commits, `64d9316`..`36ac71e`, all pushed. **Two of them change code and only ONE of the two
is fully validated** — read the status column before building on either.

### ✅ TWO GUO-YANG EQUATIONS RECOVERED — 34 of 43 (`39_2` and `14_3`)

Both were blocked by the **coprime-to-level CM filter**, not by mathematics, and both were filed
under diagnoses that had gone STALE rather than been wrong. `CMNONCOPRIME=1` (env-gated, OFF by
default) unblocks them; the published equations are what make the results believable.

    39_2   filed NONINTEGRAL.   filter on: 3 CM points vs demand 19 -> "not enough points".
           filter off: 24 points, 15 keys 0 empty. W={1} genus-7 hyperelliptic,
           IsIsomorphic to Guo-Yang in 0.06 s. Pinned in tests/GuoYangEquations.m (9 bases).
    14_3   covers under-determined by default, W={1} EMPTY (6 keys, 3 populated).
           filter off: 16 keys, 0 empty. W={1} genus-3 CRV pair, IsIsomorphic in 6817 s.
           Pinned in tests/_offline/GuoYangCurve_14_3.m -- OFFLINE because ~2 h would wreck
           GuoYangEquations.m's ~97 s.

`ModelChecks` passes both independently (82 files, 8767 checks, 0 failures) via trace-formula
point counts rather than the path that produced them.

⚠ **IS EVERYTHING REPRODUCIBLE FROM COMMITTED CODE? NO — and that is partly deliberate.** Of the
34: **24** are verified BY re-derivation (the `X0_D_N.m` tests run the pipeline, so passing IS
reproduction); **10** are stored-model comparisons that never run it. Of those 10, four do NOT
regenerate by default, for two different reasons that should not be conflated:
* `22_3`, `15_2` — ACCIDENTAL: the y2-scale guard (`1768517`) postdates the files. `Y2TWIST=1`
  restores them. This is drift, and it went unseen because every test reads the stored model.
* `39_2`, `14_3` — DELIBERATE: they exist only under `CMNONCOPRIME=1`, which stays off because it
  has no theoretical guarantee. Their justification is the published equation, not regeneration.
  Enabling the flag to make them "reproducible" would trade a documented gap for an undocumented
  risk on every base.
⇒ The target is NOT "everything regenerates by default", it is "every non-reproducing file has a
recorded reason and an independent validation" — which is what `MR_KNOWN_DRIFT` and the model
headers now encode. (`14_5`, `55_1`, `21_2`, `87_1` were still unchecked when this was written.)

### ✅ X_0^39(2) RECOVERED — the first of the two

**The one coverage result of the session.** `39_2` was filed as the NONINTEGRAL malformed-form
base; that was wrong. It is starved by the **coprime-to-level CM filter**: with the filter on it
sees 3 CM points against demand 19 and dies with "Could not find enough points"; with
`CMNONCOPRIME=1` it sees 24 and builds cleanly (15 keys, 0 empty). Its `W={1}` genus-7 curve is
**`IsIsomorphic` to Guo-Yang's published equation** (0.06 s), now pinned in
`tests/GuoYangEquations.m` (9 bases); `ModelChecks` passes it independently (82 files, 8573 checks,
0 failures, via trace-formula point counts rather than the path that produced it).

⚠ **The flag is NOT safe by default and is not enabled.** The `p | gcd(d,N)` local factor has no
live implementation (`kappaminuszero` is dead code), and at `26_3` two non-coprime discriminants
give provably wrong values. What makes `models_39_2.m` trustworthy is the INDEPENDENT ORACLE, not
the flag — **any further base produced this way must clear the same bar before being committed.**
⚠ That file does not regenerate by default; recorded in its header and in `ModelRegen`'s
`MR_KNOWN_DRIFT`, with the reason distinguished from the three y2-guard entries.

`26_3` is the same story but NOT yet closed: a model is produced, its `[1,78]` cover is verified
against Guo-Yang, but the genus-5 `W={1}` comparison had not returned when this was written.

### The Guo-Yang picture, measured rather than inherited

    43   published equations                    (see the counting trap below)
    34   we reproduce today, with a test        (24 pipeline + 9 GuoYangEquations + 14_3 offline)
     9   the gap: 8 with no model, + 22_5

`PLAN.md`'s COVERAGE section had `24` tested and a `10`-base transcription gap; both were stale.
**TIER 1' is finished as a transcription task** — `57_1` was the last transcribable base
(`64d9316`, a paired presentation like `21_2`; 8 bases, 112 s, still dominated by `21_2`'s ~100 s).
`14_3` and `22_5` are NOT transcribable: we do not possess the object to compare, so they are
model-GENERATION items. `10_19` looks like a gap in the stored models but is not — its
`X0_10_19.m` re-derives the curve via `AllEquationsAboveCovers` instead of reading a model file.

⚠ **Counting trap: the obvious grep for the 43 bases returns 41.** Two rows write the label
without braces round `D` (`$X^6_0(17)$`, `$X^6_0(29)$`), so a pattern anchored on `X^{D}_0(N)`
drops exactly those two and yields a plausible 41. Cross-check on the equation cell instead:
`multirow{1}{*}{\text}` occurs 43 times. Also, `6_17`/`6_29` appear ONLY in CM-value captions
elsewhere — having a `tests/X0_6_17.m` does not imply a published equation — and `15_1` has a test
but is not a GY equation base at all.

### Two code changes

| commit | change | status |
|---|---|---|
| `d9b52d0` | `BorcherdsForms`: shift the oo-side basis by its own valuation, not the 0-side `n0` | fix landed, **NOT yet shown to unblock any base** |
| `36ac71e` | `EquationsCovers`: `Y2TWIST=1` prototype, decide the unpinned twist instead of dropping the cover | works, and **yields 0 new GY equations** |

**The vx fix (`d9b52d0`).** The "vx class" crash is Magma's own `assert vx ge 0`
(`GalFldFun.m:305`) reached from `AbsEltseq` on a deep Laurent pole (`93_1`: `q^-60`). Cause is a
wrong-object normalisation at `BorcherdsForms.m:771`: `ech_fs_oo` holds the **oo**-expansions of the
**ZERO**-side etas `ech_etas_0`, but the shift applied was `n0`, which comes back from
`WeaklyHolomorphicBasis(... : Zero, n0 := n0)` and bounds the 0-side, not the pole at oo. Shift by
`max(n0, -min valuation)` instead, and carry the same `n_oo` into BOTH places that must agree:
`coeffs_to_divisor_matrix(-n_oo, ...)` (the shift DEFINES the column↔exponent mapping) and
`min_m := Minimum(min_m, -(n_oo + k - 1))` (`relevant_ds` must stay a superset of
`relevant_ds_0_oo`).
⚠ **Both of those were learned by running it, not by reading it.** Missing the `min_m` one made
`95_1` clear the assert and then die at `:891` with `column index not in [1..37]`; and `n_oo` was
unassigned for even `D` (the block computing it is odd-`D` only), caught by the regression.
**Safety property:** where every oo-pole already fits within `n0` — every base that currently works
— the maximum IS `n0` and the change is a literal no-op.
✅ **SWEEP DONE, AND THE FIX IS EXONERATED — do not re-open this.** The 8-base sweep came back
**6 IDENTICAL / 2 DIFFERS** (`22_3`, `15_2`). ⚠ The `DIFFERS` does NOT falsify the no-op claim, and
the instruction that stood here ("if any base says DIFFERS the commit needs revisiting") would lead
you to exactly the wrong conclusion. The discriminating test is **fresh-vs-fresh**: regenerate at
`d9b52d0~1` and compare to regenerating at HEAD. Both came back `IDENTICAL`, so the fix changed
nothing; the two bases differ because their COMMITTED models are stale (see below). No-op verified
on 7 of 7 testable bases, including odd `D` (`51_1`).

**The `Y2TWIST` prototype (`36ac71e`).** `find_y2_scales` cannot always pin the y2-scale from
sparse CM data, so `EquationsOfCovers` force-defers the cover (issue #36, `1768517`) and back-fill
usually cannot recover it. But the twist is decidable by machinery INDEPENDENT of the
Borcherds/Schofer path that produced the equation — the Eichler-Selberg point count that
`ModelVerification.m` runs as check [4]. Env-gated, off by default, and it accepts only when
exactly one squarefree twist survives at 3+ good primes, so it never trades a deferral for a guess.
Ground truth: at `22_5` it recovers `W={1,2,5,10}` at `d=1` and reproduces the committed
polynomial coefficient-for-coefficient. **So `models_22_5.m` is regenerable from current code with
the twist VERIFIED rather than trusted.** Scope, honestly: 1 cover of 4 withheld, 0 new equations,
`W={1}` still empty; the other three failed their SOLVES, so they are under-determined like
`14_3` — a different problem that this does not touch.

### THREE COMMITTED MODELS DO NOT REGENERATE — and now there is a test for it

Found while validating the vx fix, not looked for. `models_22_5.m`, `models_22_3.m` and
`models_15_2.m` do not reproduce from current code, all for one reason: the unpinned-y2-scale
guard (`1768517`, 2026-08-24 19:20) POSTDATES all three, so regeneration withholds covers they
contain. ⚠ **They are NOT wrong** — all three pass `ModelChecks` and their Guo-Yang comparisons.
They are *unreproducible*, which is a different failure and one nothing in the suite could see:
`ModelChecks` and `GuoYangEquations` read STORED models and never run the pipeline, and the
`X0_D_N.m` tests run it for only ~25 bases, each needing hand-written cover/AL data.

**`tests/_offline/ModelRegen.m` (`afa0412`) closes that gap** — auto-discovering, no per-base
authoring, works for all 81 models; regenerates and checks each committed cover is still produced
and still ISOMORPHIC. The three are listed in `MR_KNOWN_DRIFT`, reported rather than asserted away.
Two traps it cost: matching must be a **multiset** match (the first draft passed `22_3` clean while
it had LOST a cover, because two committed entries matched the same survivor), and selection must
be an **env var** (`MODELREGEN_BASES`) because `run_tests.m` `eval`s test files and a `name:=value`
argument is invisible there — it silently runs the default list instead.

**What `Y2TWIST=1` restores** (measured against the committed files): `22_5` FULLY (3/3,
coefficient-for-coefficient), `15_2` FULLY (12/12 keys, one cover differing in presentation but
`IsIsomorphic`), `22_3` 13/14. ⇒ The residual gap is **GENUS 0 by construction** —
`select_y2_twist` skips `X`g lt 1` because `HyperellipticCurve` needs degree >= 3, and both `22_3`
losses are conics. **Extending twist selection to conics is the next concrete step**, with `22_3`
as its regression target. This revises the scope note above: `Y2TWIST` is a REPRODUCIBILITY fix for
the model corpus, not the single-cover curiosity the commit message describes. Still 0 new
Guo-Yang equations.

### `26_3`: the Mobius anomaly is an `s` <-> `s~` SWAP

At discs `-267` and `-708` Guo-Yang's `s` sits in **our `s~` row**; the other 12 of 14 are correct,
and `s + s~ = 1` holds at every disc. The exact `z -> z/(z-1)` is just how an `s -> 1-s` swap looks
after the checker's cross-ratio normalisation — the involution was the shadow, not the cause.
⚠ NOT a CM-point selection ambiguity (the old framing): both values are the same point, and each
disc appears exactly once in our table. **Root cause: `s + s~ = 1`, the relation used to pin the
pair, is SYMMETRIC under exchanging them**, so it cannot resolve the ordering; the signs are forced,
only the labelling is free. Deliberately NOT fixed — what pins the ordering at the other 12 discs is
unidentified, and a tie-break without that invariant is a guess. Memory: `26-3-hauptmodul-swap`.

### Corrections made this session — do not re-derive these

* **`22_5` is not "unreproducible".** I claimed a fresh run drops `[1,2,5,10]` because the target
  cover sets differ and `08ce5fa` came from a lost path. **False** — `{1,2,5,10}` (label 7584,
  g=1) is in `GetHyperellipticCandidates()` and `Xstar`CoveredBy` today. The cover is withheld
  **on purpose** by the y2 guard, which postdates the model file by 12 hours.
* **And the committed entry is CORRECT**, measured: `VerifyModelSet` passes it 24/24 and
  discriminates the twist — six twists `d = -1,2,-2,5,-5,11` all fail with 3-5 failures. So that
  guard is CONSERVATIVE, which is what motivated `Y2TWIST`. "Do not overwrite `models_22_5.m`"
  still stands; the reason changed.
* **The vx crash is at `BorcherdsForms.m:771`, not `ShimuraQuotients.m:1116`.** My first candidate
  was the unnormalised `denom` in `IsHyperelliptic`; the traceback never reaches it. `:615` is
  exonerated by its own `assert minval eq -Minimum(...)`.
* **`111_1`/`119_1` were attempted 2026-09-03, AFTER the 66x speedup (`04f1d7b`, 08-29).** So
  "re-run them now that the basis step is faster" is NOT free progress.
* **`cmsupply`'s `CMVERD OK` does not apply to the full curve.** It iterates `Xstar`CoveredBy`
  (`ShimuraQuotients.m:1526`), the immediate covers, whose genera top out at 1 (`14_3`) and 2
  (`22_5`) — while GY's published curves there are genus 3 and 5. `OK margin 0` means "adequate for
  the easy targets, zero slack", nothing about `W={1}`.

### Runs in flight when this was written

**On lovelace, and note there are now TWO checkouts there.** The long runs use
`~/shimura/ShimuraCurveALQuotients` (at `05471c8`, behind `main`); the new work uses a SEPARATE
clone `~/shimura/vxfix` pinned to `36ac71e`, deliberately, because `AttachSpec` loads packages
lazily and pulling under a 17-hour run could swap code mid-flight. **Do not `git pull` the first
one while those jobs are alive.**

* `34_11` — ~17 h, inside `AllEquationsAboveCovers` (4 ambiguous-sign points, 16 combinations),
  RSS plateaued ~22 GB. Past `M0MultiplierExact` and `ValuesAtCMPoints` entirely. This is the run
  that would give a SECOND base ever to produce models.
* `10_61`, `14_43` — ~15 h, still in the absolute-values phase. No gate failures anywhere.
* `93_1`, `95_1`, `159_1` — the vx bases, in `~/shimura/vxfix`, output to `~/shimura/vxout`.
  ⚠ Ran locally first; that was a mistake — the Mac has 48 GB and this class peaked at 40.6 GB on
  `119_1`. Use lovelace for these.
* the 8-base regeneration sweep, `~/shimura/regen/`.

⚠ **Clearing the vx assert is necessary, not sufficient** — these bases may still die downstream,
and `95_1`/`159_1` sharing `93_1`'s cause is inherited from memory, not measured.

## Update — 2026-09-04: `tier1-models` is RETIRED; `main` is the only code branch

`main` was fast-forwarded to `tier1-models` (`475e72b`) — a clean FF, `main` was a strict
ancestor 40 commits behind with nothing of its own. **`main` is no longer "code only": it now
carries `paper/` too**, so every description of the split below is historical. The open SHIP
question "merge `paper/`, or write down that the split is deliberate" is answered by merging.

What that merge carried, beyond the paper: the per-coset `tau` fix in `M0MultiplierExact`
(`34_11` goes from failing to passing, and its multipliers match the Prop 9.15 closed form 9/9 —
see `PLAN.md` REPAIR), `tests/KudlaYangLocal.m`, and 34 Guo-Yang CM-value tables as offline tests
under `tests/_offline/` (verified after the merge to emit zero CI targets).

**Then `tier1-models` was retired**, closing that question: deleted local and remote (it was a
strict ancestor of `main`, so `git branch -d` accepted it — nothing lost), and
`worktrees/mainport` removed as redundant. The layout is now just:

    .                    main  (this checkout)
    worktrees/campaign   m0-theta-campaign

⚠ **Everything below describing a `main` / `tier1-models` split is HISTORICAL.** `main` is not
"code only" any more; commit to it directly and do not recreate `tier1-models`.

Also confirmed 2026-09-04 while checking: the branch housekeeping this file and `PLAN.md` list as
open is already DONE — only `main`, `m0-theta-campaign` and `whbasis-speedup` exist, and the nine
retired branches (`non_optimal`, `odd_DN`, `pointlessconics` and the six `SCRATCH:` ones) are all
preserved as `archive/<name>` tags **on `origin`**, so `fix-15-2-find-signs`'s 3 never-pushed
commits are safe.

## Update — 2026-09-03 (evening)

Not a full rewrite; the record below (08-30) still stands for the model-pipeline arc. This adds
what happened on the `A_m` theorem (`MAIN LINE`) and two pieces of housekeeping.

**Worktree layout changed.** `worktrees/campaign` and `worktrees/mainport`, nested inside this
checkout — not sibling directories (`-campaign`, `-mainport`) as everywhere below still calls
them. See `CLAUDE.md` for the convention on new worktrees. Committed `827266b`.

**New CI test, `tests/KudlaYangLocal.m` (`a2e888c`).** Extends the existing Prop 5.4 check
(`mu = 0`) to Prop 5.5 (nonzero isotropic, `N`-only-supported coset): the level-prime local
Whittaker *value* there is the constant polynomial `1` for every `m` tested, every base tested
(1440 checks, including `N = 2`). This is a load-bearing negative result, not a tidy-up: it rules
out "insert KY Prop 5.4/5.5 into a level-`N` analogue of Theorem 8.1" as the source of `A_m`.

**Three derivation routes for `A_m` closed, all with reasons now on record (see `PLAN.md`, MAIN
LINE, for the full writeup):**
1. KY Prop 5.4/5.5 insertion — refuted above.
2. Schwagenscheidt's oldform relation (`eq:oldform`, `sec:ident`) — **the paper itself already
   tried this and documents why it fails** at weight 3/2 (needs analytic continuation, hits a
   non-holomorphic term; numerically witnessed as `-5` where the identity forces `0`). Missed on
   first read this session; re-read `sec:ident` before re-attempting anything like it.
3. The `s`-law / genus-theta closed form (`sec:slaw`) — this is the derivation *behind*
   `prop:closedcoef` (`-a_E(m)`), i.e. the already-refuted scalar route under a different name.

**The `rem:gauge` ambiguity is real, and provably resolvable — but the resolution needs data that
doesn't exist yet.** `-a_E(m)` and "Table A" disagree at 6 indices on `X_0^{15}(2)`
(`m=1,2,3,10,15,30`) while both fit the 9-form panel (rank 4 of 6, 2-dim kernel; `gauge152.py`).
The 158-monomial data already computed in `cusp7_15_2.out` (campaign branch) gives that same
6-index matrix **full rank 6** — this specific ambiguity is not one of `sec:exact`'s "50
unremovable" directions. But a naive per-monomial solve is invalid: `c_{eta*}(0)` sums over every
cusp class, and individual eta-monomials — unlike genuine Borcherds principal parts, protected by
[GY, Lemma 24] / `prop:nohalf` — can carry a nonzero constant term at *intermediate* cusps that an
`(A_m, B_j)`-only model never sees. Restricting to combinations with zero constant term at the
intermediate classes `cusp7_15_2.out` actually recorded (`g=3,4,5,12,15,20`, from its `PP` dump)
still gives rank 6, but the numeric solve against real `c_{eta*}(0)` remains inconsistent — because
that dump is missing 4 more intermediate classes (`g=2,6,10,30`), never captured by `cusp7.m`'s
first-encountered-per-class logic. **Next step: re-run a `cusp7.m`-style pass that guarantees every
intermediate class gets dumped, then redo the reachable-subspace solve.**

**Worktree/branch housekeeping, surveyed but not executed** (see `PLAN.md`, HOUSEKEEPING, for the
full list): the three stale remote-only branches (`non_optimal`, `odd_DN`, `pointlessconics`) are
still there, plus six more local branches missed by the 09-02 cleanup because they're genuinely
unmerged (mostly `SCRATCH:`-prefixed, predate 09-02) — `fix-15-2-find-signs` notably has 3 commits
that never even reached its own remote. One untracked stray file in `worktrees/campaign`
(`polymake/nmzsolve.err`, harmless).

    main               1c53865   speedup + zero-skip + hoist + IntegralSolution + Targets
                                 + slash-constant tolerance + pointless-conics guard
                                 + models_58_5.m + 3 new tests + cache (394 files)
    tier1-models       (merged)  carries the paper work; `origin/main` merged IN on 2026-09-02,
                                 so it now has the full code side too. The two had diverged
                                 27/27 on a clean split: paper/ on this side, all code on main.
    m0-theta-campaign  9059bb0   research branch: triage results, probes, predictors
    odd-d-zeroskip     b7067c3   MERGED to main; branch kept as the CI-green record
    odd-d-invariant-hoist 4c29d1e MERGED to main (afb80b2); correctness/clarity, ~2%
    intsol-optin       969fa85   MERGED to main; CI green
    fix-pointless-conics-empty 6071772  MERGED to main (133de9c); CI green
    whbasis-speedup    624b68e   MERGED (cherry-picked as 04f1d7b); branch kept likewise

### Where the model pipeline actually stands (2026-09-02)

**Four verified models exist** — `data/models/models_58_5.m`, `ModelChecks` 48/0. That is the
*only* base that has produced models. Everything else attempted since failed:

    34_11   gates 1-3 OK, then class-constancy (GENUINE, 43% of scale -- see below)
    74_5    gates 1-3 OK, capped at 5 h in the CM-value stage, no verdict
    74_3    no longer crashes (guard merged), but yields 0 keys
    10_61   slash-constant check, at the 1e-15 calibration
    14_43   slash-constant check, at the 1e-15 calibration

**⚠ Gate 3 is NOT fully solved.** The absolute→relative fix unblocked `58_5` decisively, but
`10_61` and `14_43` still fail at `1e-15` — the siblings' own calibration. **Do not make a third
tolerance change**: either their constants genuinely disagree (real mathematics) or something
else is wrong. A `GATE3B=1` measurement on `10_61` was running on lovelace when this was written;
its output lands at `~/shimura/models/10_61.gate3b.log` there.

**Gate 4 is a genuine violation** and is characterised: `vvdata/weyl-campaign/gate4/` on the
campaign branch. Deviation tracks the CLASS, pointing at the cusp-class partition.

**Both cheap-predictor routes are closed by measurement** for the 81 unclassified bases: route A
still costs a full `WeaklyHolomorphicBasis` (20-24 min and rising), route B's `k = 3/2` phase is
known wrong. `deficit.m` has been REPAIRED and validated (`38_5` -> 1 at every pole order).

Worktrees: `-campaign`, `-mainport` (main), `-spanprobe` (**THROWAWAY**, carries the live
instrumentation — the template for the next profiling pass). `-whspeed`, `-oddd` and
`-diagnostic` were removed on 08-30; see "Housekeeping" below for what was rescued first.

---

## What landed on 2026-08-30

**1. A CI failure, fixed.** `95bd502` ("PROTOTYPE (do not merge): integrality as an acceptance
criterion") was an ancestor of the campaign tip. It replaces the divisor solve's solution and
rejects triples, so `fs[-1]` becomes a *different* form and the reference comparisons in
`SchoferIsometry.m` (Guo–Yang Table 45) and `VectorValuedForm.m` (15_2 multiplier) fail.
Reverted in `badfe5d`; preserved as `vvdata/weyl-campaign/intsol-acceptance-criterion.patch`.
The prototype's *finding* stands (33_2 does go integral under it) — what was reverted is
shipping it as unconditional pipeline behaviour.

**2. "Failed to find all Borcherds forms" is the genuine BORCHERDS OBSTRUCTION.** Not a bug, not
a too-small space. At `38_5` the rank deficit is exactly 1 and **invariant** under deepening the
pole order (bump 0→8: rows 164→172, cols 36→38, rank 35→37); the annihilator φ is stable under
enlargement (hence a fixed modular object) and φ(target) = −22 ≠ 0. Bases that **succeed** have
deficit 0 at every key and find all forms on the first triple (`34_3`, `38_7`) — so the models
in CI never meet this because they have no obstruction space at all. **Not** a level threshold:
`38_7` is larger than `38_5` in every dimension and still surjective.
⇒ **Neither more divisor triples nor deeper poles can ever help an obstructed base.**

Untested escape hatch: φ(target) is **even** and gcd(φ) = 1, and a double cover depends on its
branch divisor only mod 2, so an even correction kills the pairing without changing the cover.
Unresolved: whether that introduces an unramified quadratic twist, and the exact-divisor
`assert` would need relaxing.

**3. A 66× speedup, merged to `main`.** Two changes in `WeaklyHolomorphicBasis`: select a
spanning row subset mod p before the echelon, and skip zero terms in the basis reconstruction
(the latter was the real bulk). At `38_5`: echelon 755 s → 1.5 s, basis 810 s → 12.3 s,
**end-to-end `ppint.m` 856 s → 30 s with the verdict unchanged.**

**4. Wave-4 triage: 11 of 18 previously-unmeasured bases recovered.**

    9 form-failure  10_43 10_47 142_3 22_23 46_11 74_7 82_5 86_5 94_5
    1 NONINTEGRAL   14_37        1 assertion  115_2
    7 still TIMEOUT 65_2 6_73 77_2 85_2 91_2 119_2 146_3

Backlog tally (wave 3 → now): form-failure 19 → **28**, NONINTEGRAL 20 → 21,
TIMEOUT 18 → **7**, assertion 4 → 5, INTEGRAL 4, CM-starved 1.
Results in `vvdata/weyl-campaign/triage-wave4/` on the campaign branch.

**5. The solution cache is unified on `main`** (`f90c441`): 43 new Normaliz solves from wave 4
extend the frontier from **M ≤ 1212 to M = 2260**, plus the two M = 1236 files that had been
committed to `m0-theta-campaign` only. Cache on main: 333 → 376 files. It *was* split across
branches — with the Normaliz backend an uncached level silently costs a full solve rather than
erroring, so a split cache produces confusing "why is this slow" sessions.

---

## 6. The odd-D branch, PROFILED — and one line was 93% of it

`BorcherdsForms.m:817` on `main` was the **last** `T[i][j]*` recombination in the file lacking
the `| T[i][j] ne 0` guard (others: ~437, ~480, ~618). It sits in the **odd-D-only** 0-side
block — exactly why the 66× speedup helped even D and left odd D untouched.

    stage             65_2    77_2    85_2
    oo basis           9.3     8.4     4.0
    0-cusp basis       3.6     2.4     4.0    <- REFUTED as the cause, <=3%
    CM points          0.10    0.11    0.09   <- REFUTED, negligible
    everything after  1787    1789    1792

`T` there is a pure SELECTION matrix (measured `nnz(T) = Nrows(T)`, one 1 per row), so the
unguarded sum did 79–426× more `EtaQuot` arithmetic than needed. **Fixed in `b7067c3`, now
MERGED to main as `619051a`, measured ~140×**: `etarecomb` 1.515 s → 0.0106 s per call;
`zside` 1539 s / 3 passes → 45.0 s / 4 passes. Tests green (incl. `15_2`, odd D).

**It does not make `65_2` / `85_2` complete.** The dominant cost is now
`basis_of_weakly_holomorphic_forms(... : Zero)` — real work, steep in pole order: 2.25 s @130,
26.45 @325, 72.85 @455, **556.12 @845**; `65_2`'s last m implies pole order 8450. Constant
factor removed, ceiling unmoved.

**~~The larger win, NOT yet done~~ — DONE, and REFUTED as a lever** (`odd-d-invariant-hoist`,
`4c29d1e`). The invariance is real: the 0-side block was recomputed 336× per `m_idx` pass at
`65_2` (= 8·7·6 triples, one key each before the break) and 210× at `85_2`, and it now runs
once. **But what repeats is cheap.** Instrumenting the *pre-hoist* code directly at `65_2`:

    pass 1   336 executions    1.75 s   (mean 0.0052)
    pass 2   336 executions    5.97 s   (mean 0.0178)
    pass 3   336 executions   10.02 s   (mean 0.0298)

≈ 17.7 s over three passes against an 1800 s cap, versus ≈ 0.9 s hoisted — **about 2%.**

**Why the "336× lever" claim was wrong, and the lesson.** That judgment was formed when the
block cost 1539 s / 3 passes — but *that* cost was the unguarded recombination line, and the
zero-skip removed it (140×), collapsing the block to 1.75 s/pass. The redundancy framing
outlived the fix that made it irrelevant, because nobody re-measured the block after changing
it. **A multiplier (336×) is only a lever when multiplied by something expensive; re-measure
the multiplicand after any fix that touches it.**

**The `T`-shadowing worry was vacuous.** The ∞-side `T` was never *read* — on odd D the 0-side
kernel overwrote it immediately, on even D nothing below touches `T` at all. It was dead code,
now deleted; the 0-side matrix is renamed `T_ker0` so the question cannot recur. One real trap
found while moving it: the hoisted lines must stay together, since `ech_etas_0` is sliced out of
`ech_etas_all_0` and then *replaced in place* by its own recombination — hoisting the
recombination without the slice recombines an already-recombined list on the second key, a
silent wrong answer rather than a crash.

⇒ **The odd-D constant factors are now exhausted. Everything left is the
`basis_of_weakly_holomorphic_forms(... : Zero)` ceiling above.** Do not spend more time here.

**Reclassify**: `133_2` is **not** a TIMEOUT — it fails an assertion at 168 s.

---

## Housekeeping — 2026-08-30 (later)

`odd-d-zeroskip` went CI-green (1h37m, 0 failures) and is **merged to main as `619051a`**.
Verified before the merge: all four `T[i][j]` recombination sites (437, 480, 618, 826) now
carry the zero-skip guard, and the merge brought one commit touching one file — the
`tier1-models` merge trap did not apply, because this branch was cut from `main`.

Three worktrees were retired. Everything single-copy in them was rescued to
`vvdata/weyl-campaign/` on the campaign branch first (`4b752d8`, `9059bb0`):

    bfprof.m  dsmall.m  diag_15_2.m         drivers that existed nowhere else
    bfprof-instrumentation.patch            BFPROF/BFINV timers  -- see note-probes.md
    valuesatcmpoints-characterization.patch the non-rationality probe        "
    MISSING_TARGETS.txt                     351 bases -- see note-missing-targets.md
    note-probes.md  note-missing-targets.md the caveats, which matter more than the code

**Two caveats worth carrying forward** (both in `note-probes.md`): the bfprof patch does **not**
apply to current main — its hunks mix the timers with a superseded inline prototype of the
speedup, so lift the timers by hand; and the characterization probe **cannot have run as
written**, since it patches the two-argument `ValuesAtCMPoints` at `SchoferFormula.m:1498`,
which has no `Xstar` in scope while the added lines reference `` Xstar`N ``. Any conclusion
attributed to that probe is unevidenced.

**`-whspeed`'s 61 uncommitted files were discarded, having been shown worthless**: 43 were
polymake solutions byte-identical to main's, and the 17 `data/curves_after_*.dat` were an
*incomplete* pipeline re-run — same 18379 records, but **strictly fewer** `IsHyp`/`IsSubhyp`
determinations at every stage (−190 UpdateByGenus, −340 UpdateCurves1, −95 UpdateCurves6).
`git diff HEAD origin/main -- data/` was empty, so main was never affected.
*Method note:* `grep -v TestInWhichProved` does **not** strip the attribution — the string sits
on a continuation line, which makes that diff look like ~58k lines of content change when it is
almost entirely attribution. Parse by splitting on `*])` and keying on `CurveID`.

Also harvested: two M = 532 Normaliz solutions left behind in `-spanprobe` (`9cc771e`). The
other two at that level were already tracked, so the cache was partial exactly there — the
silent-full-resolve mode `f90c441` set out to close. Cache on main: **378 files**.

---

## NEXT — in this order

**~~1. The invariant hoist~~ — DONE and merged-pending on `odd-d-invariant-hoist` (`4c29d1e`).
Worth ~2%, not the lever. See section 6. Nothing further to do on odd-D constant factors.**

**~~1. Decide the even-correction escape hatch~~ — MEASURED, and BLOCKED on open theory.**
Full account and tooling: `vvdata/weyl-campaign/even-correction/` on the campaign branch
(`e8d68f5`). Three results:

* **The precondition holds 28/28.** `φ(target)` is EVEN with `gcd(φ) = 1` at *every* obstructed
  base. Nothing is out of reach on parity grounds. The probe aborts at the first failing key
  (the deficit is invariant across triples), turning each base from a 900–1700 s exhaustive
  failure into one key — `38_5` reproduced exactly in **29.8 s instead of 860 s**.
* **CORRECTION to "the deficit is exactly 1"**, which was measured at `38_5` alone: `166_3`,
  `22_19` and `74_7` have a **2-dimensional** obstruction space and need a simultaneous
  2-condition solve, not a single shift. Both values are even in all three.
* **The correction is constructible but unusable.** A positive control at `34_3` (whose baseline
  reproduces the committed model exactly) builds forms with `div_f` exactly `ram + <disc,2>` at
  every key, then dies in `ValuesAtCMPoints`. Diagnostic: baseline **0** non-rational cells,
  perturbed **17**.

**Why it is blocked.** The mechanism is the `KNOWN DEFECT` at `SchoferFormula.m:589` — `Kappa0`
returns a zero log-`N` coefficient at firing discriminants where it should return `A_m`. **The
preprint does NOT supply `A_m`**: `prop:closedcoef` gives the *scalar* `a_E(m)` for all `m`, but
it reproduces only **1 of 13** measured `A_m`, and structurally so — `a_E` carries the embedding
support rule and vanishes exactly where `A_m` is nonzero (`15_2` m=2; `21_2` m=2,6,18). Per
[[b-eisenstein-coefficients-solved]] the relation is `A_r = -b^{η*}_0(r)/4` with `b` the
**vector-valued** coefficient at a **nonzero isotropic coset** (support `N | r`) — a different
object, and one that no product of local densities reproduces under any convention.

⇒ **The next theorem is: general `m` at a nonzero isotropic coset.** The preprint has `m = 0`
there (`prop:kappa0`) and all-`m` for the scalar (`prop:closedcoef`); `A_m` needs the
intersection. Until that exists this hatch cannot be finished, so **do not re-attempt it as an
implementation task** — and note the same defect is what the `coprime_to_level` filter
(`ShimuraQuotients.m:1420`, self-described as "a blunt instrument") already works around.

**1. The NONINTEGRAL class is mapped, and blocked at a THIRD gate.** Full account and tooling:
`vvdata/weyl-campaign/intsol/` on the campaign branch (`fe39373`, `c27707c`, `144b20f`); code on
branch `intsol-optin` (`6a6267c` the opt-in parameter, `4cdf1fb` the `Targets` threading), CI
clean. What was established:

* `IntegralSolution := false` makes the reverted August `intsol` finding usable — the prototype
  only ever failed because it shipped *unconditionally* and changed `fs[-1]` on working bases.
  It rescues **7 of 18** measured bases, so **≈39% of the class was a choice artifact** (which
  point `Solution` returned from `sol + Kernel`), not a divisor defect.
* **But none of the seven yields a model.** Six die of CM starvation, and **`cmsupply.m`
  predicted every one** at `ppint` cost — a 7/7 validation. **Run `cmsupply` FIRST on this
  class.**
* Genus-capping via `Targets` **clears** the CM gate (`58_5`, `74_5` ran 18 and 37 min instead of
  dying at it), but then `34_11`, `58_5` and `74_5` all fail the **`M0MultiplierExact`
  slash-constant two-point check** — three bases by two independent routes.

    integrality  →  CM supply  →  M0MultiplierExact slash-constant check

* **GATE 3 IS FIXED, and it was a miscalibrated tolerance** (`79d4e89` on main). The
  slash-constant check compared two evaluations with an **absolute** `1e-30` while the other four
  guards in `M0MultiplierExact` are relative — the one site the merged
  `m0exact-relative-tolerance` work never reached. Since `absdiff = reldiff * |k|` and `|k|`
  spans ten orders, it failed on LARGE constants whose agreement was unchanged.

⇒ **FOUR VERIFIED MODELS EXIST**: `data/models/models_58_5.m` (`3a50fa6`), the first output of
this entire line of work. `ModelChecks`: **48 checks, 0 failures**. `X_0(58,5)*` needed all three
fixes together — `IntegralSolution`, the `g ≤ 2` genus cap, and the tolerance — and fails without
any one. It is a **partial set by construction** (4 of 7 covers; the header says so).

**⇒ THE NEXT GATE IS 4, NOT 3.** `34_11` clears gates 1–3 and then fails
**class-constancy** (`dev 0.0186, scale 0.0430` — a 43% deviation). That is **NOT** another
tolerance: the check's own comment records roundoff at `1e-22` for the deepest known base and
states that *a genuine violation is O(scale)*. Loosening it would manufacture a multiplier wrong
by 43%. Treat as a real defect in the m=0 assembly at that base.

    integrality → CM supply → slash-constant (FIXED) → class-constancy (open, real)

**THE LESSON, and I got this wrong twice.** *Achievable precision is base-dependent*: at the same
`Prec := 80` the two evaluation points agree to **33 digits at `58_5` but only 18 at `34_11`**
(longer eta products, more accumulated rounding). I first made the guard relative but kept
`1e-30`, calibrated on `58_5` alone — that passed `58_5` and still blocked `34_11` and `74_5`,
which were otherwise ready. The siblings' `1e-15` is the right calibration and still catches what
the guard is for (a wrong constant differs at O(1), not in the 19th digit).
**Do not re-tighten a tolerance on the evidence of one base.**

Also worth carrying: losing 60 of 80 digits at `34_11` is real precision attrition — the first
place to look if a model from these bases ever appears suspect.

Two caveats from the earlier work still stand: `IntegralSolution` is **not monotone** (`69_2`
gets four orders of magnitude worse) so it must stay per-base opt-in; and `74_3`'s failure is
**unlocalised** because the driver truncated the error and had no verbosity — a bug in the
`Targets` threading is not excluded there.

**~~2. Re-run wave 4b (the 122 never-started bases)~~ — DONE.** Predictor sweep on **lovelace**
(256 cores, idle; `galois`/`verne` irrelevant, `legendre` busy and has no Normaliz, `lava` needs a
jump host). Full account: `vvdata/weyl-campaign/sweep122/` on the campaign branch (`a5d018f`).

     81 CAPPED-1h   21 OBSTRUCTED   13 VX-ASSERT   3 ASSERT
      2 NONINTEGRAL  1 CM-STARVED    1 INTEGRAL     = 122

* **TWO runnable candidates of 122**: `10_61` (INTEGRAL, CM OK margin 0) and `14_43`
  (NONINTEGRAL — the *fixable* gate — CM OK margin 0). Both at margin 0, the position `34_11`
  was in when it cleared gates 1–3.
* **⚠ THE OBSTRUCTED CLASS IS 49, NOT 28** — 21 more bases fail "Failed to find all Borcherds
  forms", which neither more triples nor deeper poles can help. **This raises the priority of the
  theory item (5) substantially: it is now worth 49 bases.**
* **`ppint`/`cmsupply` are NOT cheap predictors on large bases.** 81 of 122 gave no verdict in a
  full hour — `ppint` must build Borcherds forms before it can speak. My earlier advice "run
  `cmsupply` first, it is `ppint`-cost" holds only for bases earlier waves already reached. A
  600 s cap was strictly worse than useless: **the cap was bounding the measurement itself.**
  The remaining 81 need a genuinely cheaper predictor, not more wall-clock.

**Remote-run notes** (lovelace): Magma 2.29-9, Normaliz 3.10.2 at `/usr/bin/normaliz` — verified
to produce lattice points identical to local 3.11.1. Clone at `~/shimura/ShimuraCurveALQuotients`.
`pkill -x magma` matches NOTHING there (the binary is `magma.exe` behind a wrapper) — and
verifying a kill with the same pattern used to kill reports false success, which caused a
double-launch here. **Verify with a different pattern.**
**Pre-solve the cache first**: of the 351 bases in `MISSING_TARGETS.txt`, 328 sit inside the
committed M ≤ 2260 frontier and 23 do not — but those 23 share only **11 distinct M** (the cache
key), so it is ~22–33 solves, a bounded batch to run *ahead* of the wave rather than a silent
per-base tax inside it. That cohort is also the high-genus tail (g up to 17, CM demand
`max(2g+5)` = 39), so run `cmsupply.m` over it first — see `note-missing-targets.md`.

**3. Re-run the 7 remaining TIMEOUT bases.** Low value, and expect it to *confirm* rather than
clear them: both odd-D constant-factor fixes are in and the ceiling is untouched.
`basis_of_weakly_holomorphic_forms(... : Zero)` is steep in pole order (556 s @ 845; `65_2`'s
last m implies 8450) and `77_2` is structurally out of reach.

**4. Route B's k = 3/2 phase.**

**5. The theory item, if the paper is the priority:** state and prove the general-`m` analogue of
`prop:kappa0` — the vector-valued weight-3/2 Eisenstein coefficient at a nonzero isotropic coset.
It is the one object standing between the obstructed class (**49 bases** — see item 2) and a model, and the exact
values are already known at `15_2`, `6_5`, `10_3`, `21_2` as a regression set. **No longer low
priority: at 49 bases this is the largest single blocker in the backlog.**

---

## TRAPS — recorded so they are not repeated

* **Wave 4b's numbers are confounded; discard them.** The cap was halved (2400 → 1200 s) *and*
  concurrency raised to 8 streams (load 30 on 14 cores), so its near-total timeout rate measures
  the scheduling as much as the bases.
* **REFUTED: "the odd-D eta-quotient explosion is the blocker."** Proposed from a parity pattern
  without checking the mechanism. Measured: `nsol` ≈ 12k on odd D exactly as on even (`133_2`
  11964, `65_2` 14346, `38_5` 12784) and `t_ip` is *instant*. No explosion.
* **`deficit.m` is EVEN-D ONLY** (guard now in the file). For odd D it omits the 0-side block
  joined into `coeffs_trunc`, so deficits are **overestimates**. Tell-tale: the value *drifts*
  with pole order instead of staying invariant (`65_2`: 5, 6, 6, 9) where a real one is constant
  (`38_5`: 1 everywhere). The three validated results are all even D and stand.
* **MERGE TRAP.** `whbasis-speedup` branched off `tier1-models`, so merging it into `main` would
  have brought **26 commits** including the entire unmerged paper rewrite
  (`level-prime-kappa.tex` +1244, the PDF, `gtsweep.m`). Cherry-pick instead, and always run
  `git log origin/main..origin/<branch>` before merging anything off `tier1-models`.
* `git -C <repo> worktree add <relative-path>` resolves the path against the **repo**, not your
  cwd — it will silently create a worktree *inside* the repo. Use absolute paths.
* `magma | head` / `| tail` can hang; redirect to a file instead.
* **The `PROBESPAN` printf in `-spanprobe` ran two `Rank()` calls per key.** Any timing taken in
  that worktree before 2026-08-30 is inflated. Now gated behind `PROBE_SPAN=1` (default off).
* **A multiplier is only a lever when the multiplicand is expensive.** The "336x redundancy"
  claim was formed when the 0-side block cost 1539 s / 3 passes, survived the zero-skip that
  collapsed it to 1.75 s/pass, and was still being quoted as "the biggest remaining lever" a
  session later. Re-measure the multiplicand after any fix that touches it.
* **Measure the thing you changed, not the whole pipeline.** The 817 fix was first tested
  end-to-end with a 2400 s cap: both bases timed out before *and* after, so the test could not
  have detected the 140× win it actually produced. For partial speedups use the `BFPROF`
  per-stage timers in `-spanprobe` against the recorded baselines.
* `ppint.m`'s first `printf` fires only *after* `BorcherdsForms` returns, so an empty log tells
  you nothing about where a run is — instrument if you need progress.

## Tools (campaign branch, `vvdata/weyl-campaign/`)

    spanprobe.m  deficit.m  matrank.m  dsize.m    route A (measured deficit)
    weildim.m  weildim2.m  dsmall.m              route B (Weil rep) + its #disc_grp table
    bfprof.m                                     per-stage odd-D profiler — USE THIS to
                                                 measure the invariant hoist
    diag_15_2.m                                  non-rationality characterization driver
    span-obstruction-probe.patch                 instrumentation — THROWAWAY WORKTREE ONLY
    bfprof-instrumentation.patch                 BFPROF/BFINV timers — DOES NOT APPLY to main
    valuesatcmpoints-characterization.patch      probe — CANNOT HAVE RUN as written
    note-probes.md                               the caveats on both patches. Read first.
    MISSING_TARGETS.txt  note-missing-targets.md 351-base target list + cache-frontier analysis
    retriage.sh  wave4_*.txt                     triage driver + stream lists
    ppint.m  cmsupply.m  genmodels.m  backlog.m  earlier triage tooling

`matrank.m` records a refuted shortcut: `coeffs_to_divisor_matrix` has **full column rank**, so
the deficiency is pure Borcherds duality, not a property of the divisor matrix.
`weildim2.m`'s O(d) trace formulas are cross-checked against explicit Weil matrices at `6_1`
(all ten traces agree), but its k = 3/2 phase is wrong by very nearly **−d/6** — **do not tune
that constant to fit**; get the half-integral convention right, then check it against the
measured deficits (`38_5` → 1, `38_7` → 0, `34_3` → 0).
