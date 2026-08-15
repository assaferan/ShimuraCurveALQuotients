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
2. Take **CHIMP** `HeuristicDecomposition(C)[4][1]` — the GEOMETRIC decomposition
   descriptor over Qbar: a list of `[dim, exp]` pairs. Flatten `exp` into a dimension
   multiset `dims`:
   - `dims = [1,1]` (incl. `[[1,1],[1,1]]` and `[[1,2]]`=E²) ⇒ Jac splits geometrically
     into two elliptic curves ⇒ **X is (geometrically) BIELLIPTIC**;
   - `dims = [2]` (`[[2,1]]`) ⇒ Jacobian stays 2-dim/simple ⇒ **NOT bielliptic**.
   - **Do NOT use the library `HeuristicDecompositionFactors`**: it (a) already extracts
     `tup[1]` internally (so the old plan's `[tup[1]:tup in facts]` double-indexed and
     crashed), and (b) decomposes over the BASE field Q (`[3]`), missing splits that only
     appear over an extension. Use entry `[4]` (over the closure) directly.
   - CHIMP attaches via `AttachSpec("/Users/sachihashimoto/Repos/CHIMP/CHIMP.spec");`
     (new Mac path); `HeuristicDecomposition` at
     `CHIMP/endomorphisms/endomorphisms/magma/heuristic/Buttons.m:309`.

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

## FEASIBILITY WALL — `#div(M)` (ported from Table 6 handoff, 2026-06-12)
The `EquationsOfCovers` → polymake step is the **same** machinery as Table 6, so its
tractability law applies verbatim. Cost is driven by `#Divisors(M)`, where
`M = 4·D0 = 4·(D·N)/2^v₂(D)` — **NOT** by the pole-order bound `n` (BorcherdsForms.m's
`LP_SIZE_CUTOFF`, default 10000, only bounds `n`):
- `#div ≤ 12` (M = 4·p·q): completes reliably.
- `#div = 16–20` (N=2 / M = 8·p·q): only at small `n`; large-`n` cases hit
  `LP_SIZE_CUTOFF` and skip. (E.g. D=51,N=2 is the Table-6 `n≈78M` intractable case.)
- `#div ≥ 24` (≥3 odd prime factors of M): **OOM-kills polymake (`Killed:9`), which
  aborts the whole magma process** — not catchable by the driver's try/catch.

`ratpts_table10.m` now has a `DIV_CUTOFF := 24` guard (mirrors `ratpts_table6.m`) that
skips `#div(M) ≥ 24` groups BEFORE calling EquationsOfCovers → clean skip, not a crash.
**Impact on Table 10's 228 squarefree-N groups: only 58 are `#div ≤ 12` (reliably
tractable, 69 W-tests); 10 are `#div 16–20`; 160 are `#div ≥ 24` (hard wall on this
hardware).** A realistic full sweep is the ~58 tractable groups, not 228. A
bigger-memory box could push the `#div ≥ 24` wall; nothing else here will.

## Suggested first session on the new machine
1. Smoke-test CHIMP timing on a quiet box: `magma idx:=15 ratpts_table10.m` under
   `timeout 1800`. This both (a) validates the whole pipeline end-to-end and (b)
   should print **BIELLIPTIC** for D=6,N=29 (user-anchored expectation).
2. If fast, sweep `maxdn:=400` then widen. Watch wall time per group; the cost is
   dominated by `EquationsOfCovers` (Borcherds/CM-point step, "slow") and the CHIMP
   endomorphism computation.
3. Record every settled row OR failure reason in `ratpts_table10_results.txt`.

## Status on new Mac (2026-06-12) — pipeline VALIDATED
- Migrated & working. Three blockers fixed: CHIMP path → `/Users/.../Repos/CHIMP`;
  restored `data/curves_after_UpdateCurves8.dat` from git `03b6593^` (untracked — it was
  deliberately deleted in `03b6593`); verdict logic rewritten (see Method §2).
- **Timing measured: ~4.5 min/group** (EquationsOfCovers ~180s; CHIMP endo step <10s).
- **Anchor D=6,N=29 → BIELLIPTIC**, logged. Spec attach: no intrinsic clashes observed.
- Next: sweep the **58 `#div ≤ 12` tractable groups** (the DIV_CUTOFF guard makes a
  `maxdn`/full sweep skip the 160 OOM-doomed groups cleanly). Batch #1–#14 in progress.
