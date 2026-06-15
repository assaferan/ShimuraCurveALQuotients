# Target-restricted construction, CM-point pre-check, and per-group timeout for the ratpts sweeps

Date: 2026-06-14
Tables in scope: Table 1 (genus 1), Table 6 (genus 0), Table 10 (genus-2 biellipticity).
Table 7: out of scope (no driver, dropped for now).

## Problem

The three `ratpts_table*.m` drivers each call `EquationsOfCovers(Xstar, curves)` once per
`(D,N)` group. That call is the expensive step (`WeaklyHolomorphicBasis` + polymake lattice
enumeration + Borcherds forms + Schöfer values at CM points). Two pain points:

1. **Expensive failures.** Some groups run for minutes and then fail:
   - `SchoferFormula.m:805` — `error "Could not find enough points, sorry!"`
     (e.g. `D=26,N=5` after 876 s; `D=34,N=5` after 37 s).
   - `BorcherdsForms.m:641` — `require #pts ge 3 : "Could not find enough rational CM points!"`
     (e.g. `D=58,N=5` after 945 s).
   Both come down to: the number of CM points available is less than the number the run demands.

2. **Non-termination / OOM wedges the whole sweep.** A hung or OOM-killed polymake call
   (`Killed: 9`) is a C-level call Magma cannot catch internally, so it takes down the entire
   `magma` process and aborts the remaining groups in that run.

## Findings from the lattice spike (2026-06-14, `spike_lattice.m`)

Empirically established by reading each curve's `Covers`/`CoveredBy` (no Borcherds work):

- **X\* is genus 0** in every case checked. Its `CoveredBy` (the immediate, index-2 covers)
  is a flat "star": e.g. `D=34,N=5` → `1×g0, 3×g1, 3×g2`, all directly over X\*.
- **Every table target is an immediate cover of X\***. A target's `Covers` (curves below it) is
  **exactly `[X*]`** — nothing between X\* and the target. Higher-genus curves sit *above* the
  target (`CoveredBy`), not below. Verified for `6,29` (g2 target), `10,7` (g0 target),
  `34,5` (two g2 targets), `26,5` (g2 target).
- Therefore there is **no genus-0→1→2 prerequisite chain for the targets**, and `EquationsAboveP1s`
  (which builds a curve above an already-built cover) does **not** apply: a genus-2 target is a
  double cover `y² = f(t)` of the P¹ hauptmodul base X\*, and `f` (degree 5–6) must be solved
  directly from `2g+5 = 9` CM-point values. No climb shortcut exists.
- **`MaxNum` is set globally:** `EquationsOfCovers(Xstar, curves)` (EquationsCovers.m:178) uses
  `num_vals = Maximum([2*g+5 : g in (genera of ALL of Xstar`CoveredBy)])`. So a high-genus
  *sibling* we don't care about inflates the demand for *every* target in the group.
- **The CM-point count is cheap-ish:** `RationalandQuadraticCMPoints(Xstar : bd:=8)` for
  `D=34,N=5` ran in **27 s** and returned `#rat=3, #quad=1` (total 4) — vs the demand of 9.
  (The slow failures 26,5/58,5 burned 876 s/945 s before failing; the guard skips that.)
- **No `timeout`/`gtimeout` binary** is installed on this machine.

### Consequences for the design

- Restricting `num_vals`/the solve set to the **target** covers helps exactly when a target is
  **lower genus than its siblings**: Table 6 (genus-0 targets) and Table 1 (genus-1 targets)
  with genus-2 siblings currently fail at `num_vals=9` even though the target needs only 5 or 7.
  Dropping the siblings rescues those.
- For a Table 10 target that **is** the max-genus (genus-2) sibling, restriction gives a
  **speedup** (don't solve the other siblings) but **no rescue** — if the CM points aren't there
  (34,5: 4 < 9) the target is genuinely unsolvable by this method.
- The genuinely-unsolvable cases are handled by the **cheap pre-check**, not by reduction.
- The climb / prerequisite-closure idea is **dropped** (the spike shows it does not apply).

## Goals

1. **Construct only the target covers** (restrict `num_vals` + solve set), reducing the
   CM-point demand to what the targets need (rescues Table 6 / Table 1; speeds up Table 10).
2. **Predict-and-skip** the cases that still cannot be met, cheaply, before paying the minutes —
   scoped to the target covers.
3. **Per-group OS timeout** so a single hang/OOM cannot abort a sweep.
4. **Run the tractable groups across Tables 1, 6, 10** under the above, launched as background
   sweeps, recording success/failure of the new path distinctly.

Non-goals: `EquationsAboveP1s` / lattice-climb (does not apply — targets are immediate covers);
any refactor of `EquationsOfCovers` beyond the target restriction below; Table 7.

## Design

### Section 1 — Target-restricted construction (rescue for Tables 6/1, speedup for Table 10)

Add an optional `Targets` parameter:

```
intrinsic EquationsOfCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot]
                            : Prec := 100, Targets := {}) -> SeqEnum, Assoc, SeqEnum
```

`Targets` is a set of `W` subgroups (each a set of AL involutions, as produced by
`AllALsFromGens(gens, D*N)`) identifying the immediate covers the caller wants.

When `Targets` is non-empty (else behaviour is exactly as today — all of `Xstar`CoveredBy`):

- **Restrict `num_vals` to the targets** (EquationsCovers.m:178):
  ```
  target_keys := [i : i in Xstar`CoveredBy | curves[i]`W in Targets];
  genus_list  := [curves[i]`g : i in target_keys];
  num_vals    := Maximum([2*g+5 : g in genus_list]);
  ```
  This is the only place `MaxNum` for `AbsoluteValuesAtCMPoints` is set; lowering it is what
  reduces the CM-point demand.
- **Restrict the solve set to the targets.** The per-cover solve in
  `EquationsOfCovers(schofer_table, all_cm_pts)` (EquationsCovers.m:118) iterates
  `K_idxs := [i : i->k in keys_fs | k gt 0]` (cmtables.m:22) over *all* covers. A non-target
  sibling left underdetermined by the smaller CM-point set would throw `require #B eq 1` /
  `require #coeffs eq 1`. So restrict the solve to target covers — e.g. filter `K_idxs` (or the
  loop) to indices whose `keys_fs[i]` curve has `W in Targets`. The returned `keys`/`crv_list`
  then contain exactly the targets.
- `BorcherdsForms` and the hauptmoduls are unchanged and still correct; only the CM-point demand
  and the solve set shrink. (Cutting sibling forms out of `BorcherdsForms` is a possible later
  speedup, not required.)

Each driver passes `Targets := { AllALsFromGens(gens, D*N) : gens in gensets }` for the row.

**Correctness:** the equation produced for each target is identical to the full run's equation
for that same cover (same Schöfer constraints; we just stop over-collecting CM points for, and
stop solving, the siblings). Verified by re-running a settled case with `Targets` and confirming
the identical model/verdict (tests below).

### Section 2 — Cheap predict-and-skip guard (scoped to targets)

A helper intrinsic:

```
intrinsic EnoughCMPointsForTargets(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], Targets::SetEnum)
          -> BoolElt, RngIntElt, RngIntElt
{ Returns (available ge required), required, available. }
```

- `required := Maximum([2*g+5 : g in (target genera)])` — same quantity Section 1 uses.
- `available := #rat + #quad` from
  `RationalandQuadraticCMPoints(Xstar : coprime_to_level := true, bd := 8)`.
  bd 8 matches the highest bd the pipeline reaches (SchoferFormula.m:795), so the count is
  faithful and conservative.

Drivers call this **after** the existing `#div(M) < DIV_CUTOFF` guard and **before**
`EquationsOfCovers`. If it returns false: log `SKIP-insufficient-CM(need=k,have=m)` for each `W`
in the row and continue — no `WeaklyHolomorphicBasis` / polymake / Borcherds work.

Measured cost: ~27 s for `D=34,N=5`. Net win is large for the slow failures (26,5: 876 s;
58,5: 945 s) and break-even for already-fast failures (34,5: 37 s); always run once per group.

### Section 3 — Per-group timeout via a portable runner

`run_table.sh <table-number> <timeout-seconds> [lo] [hi]`:

- Iterates group indices (`1..N` for the chosen table, or `lo..hi`).
- For each `i`, runs the group under a **portable timeout** (no coreutils `timeout` on this box):
  ```
  perl -e 'alarm shift; exec @ARGV' <T> magma idx:=i ratpts_table<table>.m
  ```
  (`alarm` delivers SIGALRM after `<T>` s, killing the `magma` child; the wrapper exits non-zero.)
  stdout is tee'd to a per-run log; the driver appends its own verdict row to the results file.
- **Resume:** before launching group `i`, check whether that group's `(D,N,W)` rows are already
  present in the results file; if so, skip (consistent with the existing "do not re-run a settled
  `(D,N,W)`" rule).
- On timeout (wrapper non-zero / SIGALRM), append a `TIMEOUT(>Ts)` verdict row for the group's
  `W`s so it is recorded and not retried.
- A hang/OOM kills only that one `magma` process; the loop proceeds to the next group.

Add `idx:=` single-group support to `ratpts_table6.m` (it currently selects via a `CANDIDATES`
list, whereas `ratpts_table1.m`/`ratpts_table10.m` already accept `idx:=`). Give `ratpts_table6.m`
an ordered `TABLE6` list indexable by `idx`, so the runner is uniform across all three.

### Section 4 — Extend coverage + launch

- With the timeout net + reduced demand + cheap guard, the sweeps attempt every tractable group
  (`#div(M) < DIV_CUTOFF`, minus insufficient-CM skips) across Tables 1, 6, 10.
- Re-attempt the current Table 10 failures first: `34,5` and `26,5` should now log
  `SKIP-insufficient-CM` quickly; `58,5` likewise via the guard or its own CM requirement.
- Re-run Table 6 / Table 1 rows whose genus-0/1 targets were previously inflated by genus-2
  siblings; these are the cases the reduction can newly rescue.
- Launch the three sweeps as background runs; report verdicts as they land.

### Logging / recording (all three tables)

Each driver's results file gains rows that distinguish the new path's outcomes:

- `BIELLIPTIC` / `NOT-bielliptic` / model-with-point / no-point — as today (success).
- `SKIP-insufficient-CM(need=k,have=m)` — Section 2 cheap skip.
- `TIMEOUT(>Ts)` — Section 3 OS timeout.
- `FAILED:<reason>` — any residual exception (as today).
- Successful runs that went through the **reduced** path are tagged (e.g. a
  `via=targets(num_vals=k)` note in the row or run log), so we can measure what the reduction
  bought versus the old global `num_vals`.

## Components and interfaces

| Unit | Location | Responsibility | Depends on |
|------|----------|----------------|------------|
| `EquationsOfCovers(... : Targets)` | EquationsCovers.m:173 | Restrict `num_vals` to target genera | `BorcherdsForms`, `AbsoluteValuesAtCMPoints`, `EquationsOfCovers(schofer_table,...)` |
| target-key restriction | `EquationsOfCovers(schofer_table,...)` EquationsCovers.m:118 (+ `K_idxs`, cmtables.m:22) | Solve/emit equations only for target covers | `SchoferTable` |
| `EnoughCMPointsForTargets` | ShimuraQuotients.m (near `RationalandQuadraticCMPoints`, line 1403) | Cheap available-vs-required CM-point predicate | `RationalandQuadraticCMPoints` |
| `idx:=` support + ordered `TABLE6` | ratpts_table6.m | Uniform single-group interface | — |
| `Targets`/guard wiring + verdict rows | ratpts_table1.m, ratpts_table6.m, ratpts_table10.m | Pass row targets; call guard; log distinct verdicts | Sections 1–2 |
| `run_table.sh` | repo root | Portable per-group timeout, resume, record TIMEOUT | the drivers |

## Testing / verification

0. **Genus-0 reduction rescues (Table 6):** find a Table 6 row whose genus-0 target shares its
   group with a genus-2 sibling; confirm the old run demands `num_vals=9` (and may fail) while
   the `Targets`-restricted run demands 5 and produces the conic/model. Confirm the model is
   identical to any previously-settled value for that target.
1. **Reduction is behaviour-preserving (Table 10 anchor):** re-run `D=6,N=29` (settled
   `BIELLIPTIC`) with `Targets := {AllALsFromGens({3,29},174)}`; confirm identical model and
   `BIELLIPTIC` verdict, and that logged `num_vals` ≤ the old global value.
2. **Guard is faithful & cheap:** on `D=34,N=5`, confirm `EnoughCMPointsForTargets` returns
   `false` with `required=9, available=4` in ~tens of seconds, and the driver logs
   `SKIP-insufficient-CM(need=9,have=4)` instead of running to the `error`.
3. **Timeout:** force a tiny `T` on a known-long group and confirm the runner records
   `TIMEOUT(>Ts)`, the `magma` process is killed, and the next group still runs.
4. **Resume:** re-invoke `run_table.sh` and confirm already-logged groups are skipped.

## Risks / open items

- The solve-side restriction (filtering `K_idxs`/the solve loop in
  `EquationsOfCovers(schofer_table,...)`) must keep each target's own solve well-determined;
  verify against tests (0)/(1) before trusting any new verdicts.
- The CM-count guard's value is large for slow failures, break-even for fast ones; it is always
  net-positive but note per-case timing in the log.
- `run_table.sh` resume logic must parse the results file's `(D,N,{gens})` rows correctly to
  avoid both re-running settled groups and skipping un-run ones.
- The `perl alarm` wrapper kills the process group's `magma` child on SIGALRM; confirm no
  orphaned polymake subprocess survives (kill the whole process group if needed).
