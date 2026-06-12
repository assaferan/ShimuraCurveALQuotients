# Handoff: Table 10 geometric-bielliptic check

**Date:** 2026-06-12. Pausing because this box (`lovelace`) is heavily shared and
slow — multiple other users' Magma jobs are running. Resume on the new machine.

## Goal
For each curve `X = X_0^D(N)/W` in **Table 10** (`Table10_OanaFreddy.txt`, the 469
genus-2 curves that are NOT bielliptic via an Atkin–Lehner involution, authors
unsure if geometrically bielliptic), decide **bielliptic vs not**.

## Method (confirmed)
1. Build the genus-2 model `C` exactly as in `ratpts_table6.m` / `ratpts_table1.m`:
   - `Xstar` = the star curve for `(D,N)`;
   - `crv_list, ws, keys := EquationsOfCovers(Xstar, curves);`
   - target group `W := AllALsFromGens(gens, D*N)` where `gens` = the AL subscripts
     in the table's `W` column (e.g. `<w3,w29>` → `gens = {3,29}`, `D*N = 174`);
   - find `k` with `curves[k]`W eq W`, then `C := crv_list[Index(keys,k)]`.
2. Run **CHIMP** `HeuristicDecompositionFactors(C)` on `Jac(C)`:
   - factor dims `[1,1]` ⇒ Jacobian splits geometrically into two elliptic
     curves ⇒ **X is (geometrically) BIELLIPTIC**;
   - dims `[2]` ⇒ Jacobian stays 2-dim/simple ⇒ **NOT bielliptic**.
   - CHIMP attaches via `AttachSpec("/home/sachihashimoto/CHIMP/CHIMP.spec");`
     intrinsic verified present at
     `CHIMP/endomorphisms/endomorphisms/magma/heuristic/Buttons.m:328`.

## What I built
- **`ratpts_table10.m`** — the driver. Groups Table 10 by `(D,N)` (228 groups,
  429 squarefree-N rows) so `EquationsOfCovers` runs **once per (D,N)**, then
  tests every `W` in that group. Sorted by `D*N` ascending. Run modes:
  - `magma idx:=15 ratpts_table10.m`   → only group #15 (the D=6,N=29 anchor)
  - `magma lo:=1 hi:=20 ratpts_table10.m` → groups #1..#20
  - `magma maxdn:=400 ratpts_table10.m`   → all groups with `D*N <= 400`
  - `magma ratpts_table10.m`              → prints the indexed table, computes nothing
  **Always wrap in OS `timeout`** on a shared box.
- **`ratpts_table10_results.txt`** — verdict log; driver appends one tab-separated
  line per `(D,N,W)`. Header carries the running notes. Don't re-run settled rows.

## Resource ordering (project heuristics)
- N must be squarefree (40 rows excluded — non-squarefree N out of scope here).
- Larger `N` harder; larger `D*N` harder.
- Monotone rays: fixing `D` and increasing `N`, or fixing `N` and increasing `D` —
  if the first `(D,N)` is intractable, the next almost surely is; stop along that
  ray and record the reason instead of re-running.
- `N=1`/`N=2` groups tend to be sparse-CM / huge-LP (seen in Table 6) — deprioritized
  implicitly since they only appear at the small-`D*N` top when `D` is large.

## Suggested first session on the new machine
1. Smoke-test CHIMP timing on a quiet box: `magma idx:=15 ratpts_table10.m` under
   `timeout 1800`. This both (a) validates the whole pipeline end-to-end and (b)
   should print **BIELLIPTIC** for D=6,N=29 (user-anchored expectation).
2. If fast, sweep `maxdn:=400` then widen. Watch wall time per group; the cost is
   dominated by `EquationsOfCovers` (Borcherds/CM-point step, "slow") and the CHIMP
   endomorphism computation.
3. Record every settled row OR failure reason in `ratpts_table10_results.txt`.

## Open items / unverified
- **CHIMP wall time not yet measured.** First smoke test (`y^2=x^6+1`) did not finish
  within 300s here, almost certainly due to machine load — re-measure on the new box
  before trusting any per-curve budget.
- Driver was syntax-checked in table-only mode (see `ratpts_table10_results.txt`
  notes / your run log) but no `(D,N,W)` verdict has been computed yet.
- Possible intrinsic-name collisions between `ShimuraQuotients.spec` and
  `CHIMP.spec` when both attached — watch for Magma "redefining intrinsic" warnings;
  harmless unless a signature actually clashes.
