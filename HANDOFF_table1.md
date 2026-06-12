# Handoff: Table 1 genus-1 rational-point search

## Goal
For each genus-1 curve `X = X_0^D(N)/W` in **Table1_OanaFreddy.txt** (25 curves the
authors are unsure have a rational point), compute the cover model `C` and test
`IsEllipticCurve(C)`. If `C` is elliptic it **has** a rational point (resolves
`X(Q) != empty`); the hyperelliptic involution = image of Fricke `w_{DN}`.
We will NOT finish all 25 — pick off the feasible ones and record why the rest fail.

## What exists
- **[ratpts_table1.m](ratpts_table1.m)** — the genus-1 driver. Built as the analog of
  the already-working genus-0 **[ratpts_table6.m](ratpts_table6.m)**.
  - `EquationsOfCovers(Xstar, curves)` once per `(D,N)` star curve → pull model `C`
    for the target `W = AllALsFromGens(gens, D*N)` → assert `Genus(C) eq 1` →
    `IsEllipticCurve(C)`. Records min model + conductor on success, or "no rational
    point found" / error reason on failure.
  - `TABLE1` is the full 25 entries as `<D, N, gens, note>`, **ordered by feasibility**.
  - Run one: `magma idx:=N ratpts_table1.m`  (N = index into TABLE1, 1-based).
    Run all: `magma ratpts_table1.m`. Output appends to **ratpts_table1_output.txt**.
- **[ratpts_table1_output.txt](ratpts_table1_output.txt)** — running log of results/aborts.

## Method constraints / resource principles (from the user)
- **N must be squarefree** (all Table 1 N already are).
- Bigger `N` is harder; bigger `D*N` is harder.
- Monotonic: fix `D`, increase `N` — if the first fails, later ones fail too.
  Fix `N`, increase `D` — same. Use this to stop early within a `D`- or `N`-family.
- Tiny `N=1,2` are NOT easy despite small `D*N`: they hit the **sparse-CM / huge-LP**
  wall (Table 6: `D=51,N=2` LP ~78M points = intractable; `D=34,N=3` = not enough
  CM points). These are deprioritized to the END of TABLE1.
- **Cost**: `EquationsOfCovers` is the expensive step. The genus-0 analog took ~940s
  for `D*N=70` on a quiet machine; genus-1 needs more CM-point values (`2g+5=7`) so
  expect more. Run **serially, one at a time** — do not parallelize.

## CRITICAL: machine
The previous run was on a **shared, overloaded server** (load avg ~18; a 128-way
`FilterByWeilPolynomial` pipeline + many other users' Magma jobs were running).
Our single job was CPU-starved. **Move to a quieter/faster machine before resuming.**
Only ever kill our OWN processes — many other users share this box.

## Status of runs
| D | N | W (gens) | D*N | result |
|---|---|---|-----|--------|
| 6 | 35 | {10,42} | 210 | **ABORTED** after ~29 min @ 99.9% CPU, still inside `EquationsOfCovers`. NOT a math/resource verdict — abandoned only due to server load. **Retry first; do not mark intractable.** |

Everything else: **not yet attempted.**

## Recommended order to resume (priority batch, N>=5, moderate D)
1. `6,35  {10,42}`   (resume here — was in progress)
2. `34,5  {10,34}`   DN=170
3. `10,21 {5,21}`    DN=210
4. `15,14 {2,5,7}`   DN=210
5. `34,7  {2,17}`    DN=238
6. `10,33 {2,55}`    DN=330
7. `74,5  {10,74}`   DN=370
8. `10,39 {3,5,13}`  DN=390
9. `10,51 {2,15,51}` DN=510
10. `6,85 {2,15,51}` DN=510
11. `21,26 {2,21,39}` DN=546
12. `10,61 {10,122}`  DN=610
13. `6,115 {5,6,46}`  DN=690
14. `14,57 {2,21,57}` DN=798
15. `10,93 {3,10,62}` DN=930
16. `6,161 {2,21,69}` DN=966

Deprioritized (likely intractable; N=1/2 or huge D*N): `210,1 {7,15}`,
`330,1 {2,33}`, `330,1 {3,10}`, `462,1 {11,14}`, `798,1 {2,3,19}`,
`1230,1 {3,10,82}`, `1722,1 {6,14,41}`, `119,2 {7,17}`, `210,19 {6,7,10,19}`.
(These match TABLE1 indices 17–25 in ratpts_table1.m.)

## Logging rule (from the user)
Once a `(D,N,W)` is **successfully run** OR you know **why it fails** (too slow /
code errors / not enough CM points), write a one-line note (the MODEL on success,
or the reason on failure) into ratpts_table1_output.txt and the TABLE1 entry's
`note`, and **never re-run it**. An *aborted-for-load* run (like 6,35) is NOT a
verdict — it stays eligible for retry.
