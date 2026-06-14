# Target-restricted requirement, CM-point pre-check, and per-group timeout for the ratpts sweeps

Date: 2026-06-14
Tables in scope: Table 1 (genus 1), Table 6 (genus 0), Table 10 (genus-2 biellipticity).
Table 7: out of scope (no driver, dropped for now).

## Problem

The three `ratpts_table*.m` drivers each call `EquationsOfCovers(Xstar, curves)` once per
`(D,N)` group. That call is the expensive step (`WeaklyHolomorphicBasis` + polymake lattice
enumeration + Borcherds forms + Schöfer values at CM points). Two pain points:

1. **Expensive failures.** Some groups run ~15 min and then fail:
   - `SchoferFormula.m:805` — `error "Could not find enough points, sorry!"`
     (e.g. `D=26,N=5` after 876 s; `D=34,N=5` after 37 s).
   - `BorcherdsForms.m:641` — `require #pts ge 3 : "Could not find enough rational CM points!"`
     (e.g. `D=58,N=5` after 945 s).
   Both ultimately come down to the number of CM points available vs. needed.

2. **Non-termination / OOM wedges the whole sweep.** A hung or OOM-killed polymake call
   (`Killed: 9`) is a C-level call Magma cannot catch internally, so it takes down the entire
   `magma` process and aborts the remaining groups in that run.

### Why the requirement is inflated

`EquationsOfCovers(Xstar, curves)` (EquationsCovers.m:173) solves **all** immediate (index-2)
covers of `X*` at once. It sets

```
genus_list := [curves[i]`g : i in Xstar`CoveredBy];   // ALL immediate covers
num_vals   := Maximum([2*g+5 : g in genus_list]);      // CM points demanded
```

and passes `MaxNum := num_vals` to `AbsoluteValuesAtCMPoints`. So a high-genus sibling cover
we do **not** care about inflates `num_vals`, which is what trips the "not enough points"
failures. Each table row only cares about a specific set of target subgroups `W` (the `gens`
sets in `TABLE1` / `CANDIDATES` / `TABLE10`).

## Goals

1. **Reduce the CM-point requirement to the targets we actually want** (rescue + speedup).
2. **Predict-and-skip** the residual cases that still cannot be met, cheaply, before paying
   the ~15 min — scoped to the target covers, not all covers.
3. **Per-group OS timeout** so a single hang/OOM cannot abort a sweep.
4. **Run the tractable groups across Tables 1, 6, 10** under the above, launched as background
   sweeps, recording success/failure of the new path distinctly.

Non-goals: touching `EquationsAboveP1s` / the lattice-climb path; refactoring
`EquationsOfCovers` end-to-end to build only target covers (the requirement/solve restriction
below is sufficient); Table 7.

## Design

### Section 1 — Target-restricted requirement (the rescue + speedup)

Add an optional `Targets` parameter to the high-level intrinsic:

```
intrinsic EquationsOfCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot]
                            : Prec := 100, Targets := {}) -> SeqEnum, Assoc, SeqEnum
```

`Targets` is a set of `W` subgroups (each a set of AL involutions, as produced by
`AllALsFromGens(gens, D*N)`) identifying the covers the caller wants.

When `Targets` is non-empty:

- Restrict `genus_list` to the target covers only:
  ```
  target_keys := [i : i in Xstar`CoveredBy | curves[i]`W in Targets];
  genus_list  := [curves[i]`g : i in target_keys];
  num_vals    := Maximum([2*g+5 : g in genus_list]);
  ```
  (When `Targets` is empty, behaviour is unchanged: all of `Xstar`CoveredBy`.)
- The downstream per-cover solve (`EquationsOfCovers(schofer_table, all_cm_pts)`,
  EquationsCovers.m:118) must only emit/require equations for the target covers, so that a
  non-target cover left underdetermined by the smaller CM-point set cannot throw the
  `require #B eq 1` / `require #coeffs eq 1` failures. Plumb the target key set through (e.g.
  restrict `k_idxs` in the `SchoferTable`, or filter the solve loop to target keys).
- `BorcherdsForms` and the hauptmodul computation are unchanged and still correct; only the
  CM-point demand (`MaxNum`) and the solve set shrink. (Optionally, passing fewer covers into
  `BorcherdsForms` could also cut form-completion work, but that is a follow-up optimisation,
  not required for correctness.)

Each driver passes `Targets := { AllALsFromGens(gens, D*N) : gens in gensets }` for the row.

**Correctness note:** restricting `MaxNum` only lowers the number of CM points gathered; the
equations produced for the target covers are identical to what the full run would produce for
those same covers (same Schöfer constraints, just not over-collecting points for siblings).
This is verified by re-running an already-settled case (e.g. `D=6,N=29`, expected
`BIELLIPTIC`) with `Targets` set and confirming the same model/verdict.

### Section 2 — Cheap predict-and-skip guard (scoped to targets)

A helper intrinsic:

```
intrinsic EnoughCMPointsForTargets(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], Targets::SetEnum)
          -> BoolElt, RngIntElt, RngIntElt
{ Returns (available >= required), required, available. }
```

- `required := Maximum([2*g+5 : g in (target genera)])` — same formula Section 1 uses.
- `available := #rat + #quad` from `RationalandQuadraticCMPoints(Xstar : coprime_to_level := true, bd := <bd the pipeline reaches>)`.
  Match the bd the pipeline actually uses (it escalates to 8 in `AbsoluteValuesAtCMPoints`,
  SchoferFormula.m:795), so the guard is faithful and conservative.

Drivers call this **after** the existing `#div(M) < DIV_CUTOFF` guard and **before**
`EquationsOfCovers`. If it returns false: log
`SKIP-insufficient-CM(need=k,have=m)` for each `W` in the row and continue — no
`WeaklyHolomorphicBasis` / polymake / Borcherds work performed.

**Caveat to verify before committing the guard:** that this CM-count call is genuinely cheap
relative to a full run. It is arithmetic/class-number based, so expected to be seconds, but it
will be timed on `D=34,N=5` (a known fast-ish failure) first. If it turns out to be a large
fraction of the full run, keep it (it is still run once and saves the rest) but note the cost
in the results log.

### Section 3 — Per-group timeout via a shell runner

`run_table.sh <table-number> <timeout-seconds> [lo] [hi]`:

- Iterates group indices (`1..N` for the chosen table, or `lo..hi`).
- For each `i`, runs `timeout <T> magma idx:=i ratpts_table<table>.m`, tee-ing stdout to a
  per-run log and letting the driver append its own verdict row to the table's results file.
- **Resume:** before launching group `i`, check whether that group's `(D,N,W)` rows are
  already present in the results file; if fully present, skip (do not re-run settled groups —
  consistent with the existing "do not re-run a settled (D,N,W)" rule).
- On `timeout` exit code 124, append a `TIMEOUT(>Ts)` verdict row for the group's `W`s so the
  failure is recorded and not retried.
- A hang/OOM kills only that one `magma` process; the loop proceeds to the next group.

Add `idx:=` single-group support to `ratpts_table6.m` (it currently selects via a `CANDIDATES`
list / env, whereas `ratpts_table1.m` and `ratpts_table10.m` already accept `idx:=`). Give all
three drivers the same `idx:=i` interface so the runner is uniform. `ratpts_table6.m` will gain
an ordered `TABLE6` list (mirroring the others) indexable by `idx`.

### Section 4 — Extend coverage + launch

- With the timeout net and reduced requirement, the sweeps attempt every tractable group
  (`#div(M) < DIV_CUTOFF`, minus insufficient-CM skips) across Tables 1, 6, 10.
- The 5 current Table 10 failures (`26,5`; `34,5`; `58,5`) are re-attempted first under
  `Targets`-restricted `MaxNum`; record whether the reduction rescues them.
- Launch the three sweeps as background runs and report verdicts as they land.

### Logging / recording (applies to all three tables)

Each driver's results file gains rows that distinguish the new path's outcomes:

- `BIELLIPTIC` / `NOT-bielliptic` / model-with-point / no-point — as today (success).
- `SKIP-insufficient-CM(need=k,have=m)` — Section 2 cheap skip.
- `TIMEOUT(>Ts)` — Section 3 OS timeout.
- `FAILED:<reason>` — any residual exception (as today).
- Successful runs that went through the **reduced** `MaxNum` path are tagged (e.g. a
  `via=targets(num_vals=k)` note in the row or the run log) so we can measure what the
  reduction bought versus the old global `num_vals`.

## Components and interfaces

| Unit | Location | Responsibility | Depends on |
|------|----------|----------------|------------|
| `EquationsOfCovers(... : Targets)` | EquationsCovers.m:173 | Restrict requirement + solve set to targets | `BorcherdsForms`, `AbsoluteValuesAtCMPoints`, `EquationsOfCovers(schofer_table,...)` |
| target-key restriction in `EquationsOfCovers(schofer_table,...)` | EquationsCovers.m:118 | Only emit/require target-cover equations | `SchoferTable` |
| `EnoughCMPointsForTargets` | ShimuraQuotients.m (near `RationalandQuadraticCMPoints`, line 1403) | Cheap available-vs-required CM-point predicate | `RationalandQuadraticCMPoints` |
| `idx:=` support + `TABLE6` | ratpts_table6.m | Uniform single-group interface | — |
| `Targets`/guard wiring | ratpts_table1.m, ratpts_table6.m, ratpts_table10.m | Pass row targets; call guard; log distinct verdicts | Sections 1–2 |
| `run_table.sh` | repo root | Per-group `timeout`, resume, record TIMEOUT | the drivers |

## Testing / verification

1. **Reduction is behaviour-preserving:** re-run `D=6,N=29` (settled `BIELLIPTIC`) with
   `Targets` set; confirm identical model and verdict, and that the logged `num_vals` is ≤ the
   old global value.
2. **Guard is faithful & cheap:** on `D=34,N=5`, confirm `EnoughCMPointsForTargets` returns
   false quickly (seconds), matching the eventual real failure, and that the driver now logs
   `SKIP-insufficient-CM` instead of running ~37 s+ to the `error`.
3. **Rescue check:** for each of the 5 current failures, record whether target-restricted
   `MaxNum` now succeeds; if still insufficient, the guard should pre-empt it.
4. **Timeout:** force a tiny `T` on a known-long group and confirm the runner records
   `TIMEOUT(>Ts)`, the `magma` process is killed, and the next group still runs.
5. **Resume:** re-invoke `run_table.sh` and confirm already-logged groups are skipped.

## Risks / open items

- Plumbing the target-key restriction through `SchoferTable` (Section 1, solve side) is the
  most delicate change; verify against test (1) before trusting any new verdicts.
- If the CM-count call is not cheap, Section 2's value drops to "skip slightly earlier"; still
  net-positive but note it.
- `run_table.sh` resume logic must parse the results file's `(D,N,{gens})` rows correctly to
  avoid both re-running settled groups and skipping un-run ones.
