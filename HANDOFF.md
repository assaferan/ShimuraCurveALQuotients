# Handoff — 2026-08-30

**Supersedes** the 2026-07-17 handoff about producing cover models, archived as
`HANDOFF_2026-07-17.md`. That task is not dead, but it is gated on the blocker described below.

Everything here is committed and pushed. **`git pull` first — local `main` may be stale.**

**➡ For what to do next, see `PLAN.md`** — five tracks, a do-not list, and the recurring traps.
This file is the record of *what happened*; `PLAN.md` is the record of *what to do*. When the two
disagree about state, this file wins.

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
