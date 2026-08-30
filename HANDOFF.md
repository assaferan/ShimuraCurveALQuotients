# Handoff — 2026-08-30

**Supersedes** the 2026-07-17 handoff about producing cover models, archived as
`HANDOFF_2026-07-17.md`. That task is not dead, but it is gated on the blocker described below.

Everything here is committed and pushed. **`git pull` first — local `main` may be stale.**

    main               f90c441   66x WeaklyHolomorphicBasis speedup (04f1d7b) + unified solution cache
    tier1-models       a924e1a   unchanged; carries the unmerged paper work
    m0-theta-campaign  eb700b1   research branch: triage results, probes, predictors
    whbasis-speedup    624b68e   the speedup branch (CI green: 57 jobs, 0 failures)

Worktrees: `-whspeed` (speedup), `-spanprobe` (**THROWAWAY**, carries all instrumentation —
this is the template for the next profiling pass), `-campaign`, `-diagnostic`.

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
unguarded sum did 79–426× more `EtaQuot` arithmetic than needed. **Fixed on branch
`odd-d-zeroskip` (`b7067c3`), measured ~140×**: `etarecomb` 1.515 s → 0.0106 s per call;
`zside` 1539 s / 3 passes → 45.0 s / 4 passes. Tests green (incl. `15_2`, odd D).

**It does not make `65_2` / `85_2` complete.** The dominant cost is now
`basis_of_weakly_holomorphic_forms(... : Zero)` — real work, steep in pole order: 2.25 s @130,
26.45 @325, 72.85 @455, **556.12 @845**; `65_2`'s last m implies pole order 8450. Constant
factor removed, ceiling unmoved.

**The larger win, NOT yet done**: the entire 0-side block is invariant across triples *and*
keys — checksummed, one distinct signature across all 336 triples of `65_2` — yet recomputed
per key per triple. That is 336× redundancy at `65_2`, 210× at `85_2`. Before hoisting, resolve
that the ∞-side `T` (~line 844) is shadowed by the 0-side kernel matrix (~line 885); nobody has
verified the shadowing is harmless.

**Reclassify**: `133_2` is **not** a TIMEOUT — it fails an assertion at 168 s.

---

## NEXT — in this order

**1. Merge `odd-d-zeroskip` once CI is green.** Then re-run the 7 remaining TIMEOUT bases.

**2. The invariant hoist (336×)** — the biggest remaining lever on odd D. Resolve the `T`
shadowing first.

**3. Decide the even-correction escape hatch.** Worth it now: the obstructed class is 28.
The twist question and the exact-divisor assert are the two unknowns.

**4. Re-run wave 4b (the 122 never-started bases)** at ≤ 4 streams and the original 2400 s cap.

**5. Route B's k = 3/2 phase.** Low priority; route A measures directly.

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
* **Measure the thing you changed, not the whole pipeline.** The 817 fix was first tested
  end-to-end with a 2400 s cap: both bases timed out before *and* after, so the test could not
  have detected the 140× win it actually produced. For partial speedups use the `BFPROF`
  per-stage timers in `-spanprobe` against the recorded baselines.
* `ppint.m`'s first `printf` fires only *after* `BorcherdsForms` returns, so an empty log tells
  you nothing about where a run is — instrument if you need progress.

## Tools (campaign branch, `vvdata/weyl-campaign/`)

    spanprobe.m  deficit.m  matrank.m  dsize.m    route A (measured deficit)
    weildim.m  weildim2.m                        route B (Weil representation)
    span-obstruction-probe.patch                 instrumentation — THROWAWAY WORKTREE ONLY
    retriage.sh  wave4_*.txt                     triage driver + stream lists
    ppint.m  cmsupply.m  genmodels.m  backlog.m  earlier triage tooling

`matrank.m` records a refuted shortcut: `coeffs_to_divisor_matrix` has **full column rank**, so
the deficiency is pure Borcherds duality, not a property of the divisor matrix.
`weildim2.m`'s O(d) trace formulas are cross-checked against explicit Weil matrices at `6_1`
(all ten traces agree), but its k = 3/2 phase is wrong by very nearly **−d/6** — **do not tune
that constant to fit**; get the half-integral convention right, then check it against the
measured deficits (`38_5` → 1, `38_7` → 0, `34_3` → 0).
