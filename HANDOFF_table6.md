# Handoff: Table 6 rational-points search

## Goal
We are improving on Table 6 of OanaFreddy.pdf: the 48 genus-0 curves
X = X_0^D(N)/W (W ≤ W_0(D,N) nontrivial) for which the authors couldn't decide
whether X(Q) = ∅. We find explicit genus-0 models (conics) using the repo's
machinery and test for rational points. Bigger D·N is harder; large N is also hard.

- Full table: `Table6_OanaFreddy.txt`
- Driver script: `ratpts_table6.m`  (run with `magma ratpts_table6.m`)

## How the driver works
For each `<D, N, [generator-sets]>` candidate it:
1. finds the star curve for (D,N),
2. computes immediate-cover equations via `EquationsOfCovers` (slow step;
   internally calls polymake to enumerate lattice points of an LP),
3. for each target W (given by AL subscripts) builds the genus-0 model and calls
   `HasRationalPoint(Conic(C))`.

Only **squarefree N** works (method requirement). The 7 non-squarefree-N table
rows (N = 4, 8, 9, 25, 49) are out of scope.

## Results so far

| D | N | D·N | Status |
|---|---|-----|--------|
| 10 | 7 | 70 | ✅ SOLVED — all 3 groups have rational points (models below) |
| 34 | 3 | 102 | ❌ FAILED — "not enough CM points" (hard wall) |
| 26 | 5 | 130 | ❌ FAILED — "Could not find enough points" (hard wall; user confirmed skip) |
| 51 | 2 | 102 | ⛔ INTRACTABLE — polymake LP size n≈78M, far above cutoff |
| 21 | 10 | 210 | ⏸️ INCOMPLETE — was running >17 min on a slow server, killed mid-computation |

### D=10, N=7 models (all genus 0 ≅ conic ≅ P¹, all HAVE rational points)
- W=<{2,5}>:   `27/16*x^2 + 47/64*x*z + y^2 + 5/64*z^2 = 0`,  pt (-6/31 : 7/248 : 1)
- W=<{5,7}>:   `27*x^2 - 22*x*z + y^2 - 5*z^2 = 0`,            pt (-5/27 : 0 : 1)
- W=<{10,14}>: `27/64*x^2 + 5/64*x*z + y^2 = 0`,               pt (-5/27 : 0 : 1)

## Resource guards in place
- `LP_SIZE_CUTOFF := 10000` is set at the top of `ratpts_table6.m` and read by
  `get_integer_prog_solutions` in `BorcherdsForms.m`. LP instances with n above
  this skip polymake and return [] (instead of hanging). 24*n is the bounding
  polytope coefficient; max LP solved historically is n=499.
- Two distinct failure modes to expect:
  - "not enough CM points" → fundamental for that curve, don't retry.
  - LP too large → bails fast now thanks to the cutoff.

## Next steps
The candidate list in `ratpts_table6.m` (the `CANDIDATES` default) already
contains every squarefree-N row, ordered roughly by D·N. Resume from **D=21, N=10**
onward (D·N = 210, 330, 462, …). Just run `magma ratpts_table6.m` on a faster box.

Remaining untried squarefree cases include: (21,10), (15,14), (14,15), (10,21),
(6,35), then the D·N=330 tier (33,10),(22,15),(15,22),(10,33),(6,55), etc., up to
(390,7) which is almost certainly too large (LP cutoff will bail).

## Gotchas
- **Kill jobs with** `pkill -9 -f "magma.exe ratpts_table6"` — target `magma.exe`,
  NOT the bash wrapper. Killing the wrapper (`pkill -f "magma ratpts_table6"`)
  orphans the magma.exe child, which keeps running and writing to the same output
  file. We accidentally spawned 4 concurrent orphans this way and garbled the
  output; always verify with `ps aux | grep magma.exe` that exactly one (or zero)
  remains.
- Output goes to `ratpts_table6_output.txt` (overwritten each run).
- Each case took ~15 min on the (slow) `lovelace` server; faster elsewhere.
