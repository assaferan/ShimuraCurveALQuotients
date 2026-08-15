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
| 21 | 10 | 210 | ⏸️ INCOMPLETE — was running >17 min on a slow server, killed mid-computation. Tried again but it seems like the polymake keep getting killed because of too much memory? |

### D=10, N=7 models (all genus 0 ≅ conic ≅ P¹, all HAVE rational points)
- W=<{2,5}>:   `27/16*x^2 + 47/64*x*z + y^2 + 5/64*z^2 = 0`,  pt (-6/31 : 7/248 : 1)
- W=<{5,7}>:   `27*x^2 - 22*x*z + y^2 - 5*z^2 = 0`,            pt (-5/27 : 0 : 1)
- W=<{10,14}>: `27/64*x^2 + 5/64*x*z + y^2 = 0`,               pt (-5/27 : 0 : 1)

## Resource guards in place
- `LP_SIZE_CUTOFF := 10000` is set at the top of `ratpts_table6.m` and read by
  `get_integer_prog_solutions` in `BorcherdsForms.m`. LP instances with n above
  this skip polymake and return [] (instead of hanging). 24*n is the bounding
  polytope coefficient; max LP solved historically is n=499.
- `DIV_CUTOFF := 24` (added top of the `CANDIDATES` loop). Skips any star curve
  whose polymake level **M = 4·D0 = 4·(D·N)/2^v₂(D)** has #Divisors(M) ≥ 24,
  BEFORE calling EquationsOfCovers. See the #div law below.
- **Three** distinct failure modes to expect:
  - "not enough CM points" → fundamental for that curve, don't retry.
  - huge n (sparse CM, esp. N=2) → LP too large, bails fast via `LP_SIZE_CUTOFF`.
  - huge #Divisors(M) → polymake OOMs (`Killed: 9`); now a clean skip via `DIV_CUTOFF`.

## The #div(M) tractability law (learned 2026-06-12)
The polymake step enumerates lattice points of a polytope of dimension
#Divisors(M), M = 4·D0. Cost is driven by **#div(M), not by n** (the OOM-killed
M=420 case had n=145 ≪ 10000, so `LP_SIZE_CUTOFF` never fired). From completed
vs `Killed:9` polymake artifacts:
- #div ≤ 12 (M = 4·p·q): completes reliably — proven to M=1212, n=499.
- #div = 16–20 (N=2 / M = 8·p·q): only at small n (M=408 done at n=79, killed at n=143).
- #div ≥ 24 (3 odd primes in M): OOMs even at its minimum forced n. **Hard wall.**

## Next steps — essentially exhausted
Computing #div(M) over the remaining `CANDIDATES`: **all 25 are ≥ #div 16**, and
**21 of them are ≥ #div 24** (OOM wall). The only non-OOM rows are the four N=2
cases (55,2),(87,2),(95,2),(111,2), all #div 16 — but those are the deprioritized
sparse-CM family, and their cousin (51,2) already died on huge-n. So there are
**no clearly tractable Table 6 levels left** on this hardware. `DIV_CUTOFF` just
makes a `magma ratpts_table6.m` sweep skip the doomed rows cleanly instead of
dying. A bigger-memory box could push the #div-24 wall; nothing else will here.

## Gotchas
- **Kill jobs with** `pkill -9 -f "magma.exe ratpts_table6"` — target `magma.exe`,
  NOT the bash wrapper. Killing the wrapper (`pkill -f "magma ratpts_table6"`)
  orphans the magma.exe child, which keeps running and writing to the same output
  file. We accidentally spawned 4 concurrent orphans this way and garbled the
  output; always verify with `ps aux | grep magma.exe` that exactly one (or zero)
  remains.
- Output goes to `ratpts_table6_output.txt` (overwritten each run).
- Each case took ~15 min on the (slow) `lovelace` server; faster elsewhere.
