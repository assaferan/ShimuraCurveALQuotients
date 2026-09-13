# The obstructed class re-run against post-vx-fix code (2026-09-10)

**WHY.** Every "Failed to find all Borcherds forms" verdict on record was taken 2026-09-01/02
(`sweep122/`, the triage waves, the span probes). `BorcherdsForms.m` had six commits after that,
including **`d9b52d0` (09-05) -- the vx fix**, a CORRECTNESS fix to the very stage that raises the
error. So the 49-base figure justifying `A_m`'s priority in `PLAN.md` rested on pre-fix verdicts.

**METHOD.** `spanprobe.m` at `PROBE_BUMP=0` (inert; reproduces production), run from the campaign
worktree AFTER `main` was merged down, so it uses current shared-path code. Each verdict is scored
three ways -- SUCCESS / obstructed(unchanged) / OTHER -- so a base failing for a DIFFERENT reason
cannot be silently counted as "still obstructed".

**RESULT: 49 re-run, 49 still obstructed, 0 flipped, 0 OTHER.** Runtimes 18 s to 1349 s.

**THE 49 WAS RECOVERED, NOT ASSUMED.** `sweep122/SUMMARY.txt` names only the 21 NEWLY obstructed and
refers to "the known 28" without listing them. Harvesting every obstructed verdict across
`sweep122/SUMMARY.txt`, `triage-wave*/verdicts*.txt` and all `*.log` unions to EXACTLY 49 distinct
bases -- an independent confirmation of the "28 + 21" figure. `bases49.txt` is that union.

**TWO METHOD NOTES WORTH KEEPING.**
* `38_5` reproduces its RECORDED GEOMETRY, not merely its failure: `pole_order=190 pool=164` here
  against `poleord 190 rows 164` in `spanprobe-38_5-2026-08-29/`. That is what makes it the same
  computation. It also ran in **18 s against 901 s recorded** (~50x; the q-expansion bootstrap).
* ⚠ **LEVEL DOES NOT PREDICT COST**: 18 s to 1349 s, uncorrelated with `M = 2DN`. `106_5` was
  chosen first by smallest `M` and is ~5x slower than `38_5`, which additionally had a known value
  to reproduce. Rank candidates by the sweep's own log mtimes, not by `M`.
* The rank detail (`PROBESPAN`, deficit exactly 1) needs `span-obstruction-probe.patch`, which
  **no longer applies** -- `BorcherdsForms.m` moved under it. Port it if a verdict ever flips.

---

## ⚠⚠ THE 49 IS A LOWER BOUND, NOT THE SIZE OF THE CLASS (2026-09-13)

A routine backlog sweep launched `142_1` and `166_1` — neither in `bases49.txt` — and BOTH came
back `Failed to find all Borcherds forms`, i.e. obstructed.

    142_1  (D = 2*71, N = 1)   NEW, obstructed
    166_1  (D = 2*83, N = 1)   NEW, obstructed

⇒ **The obstructed class is larger than 49.** The figure is a property of WHICH BASES HAVE BEEN
SWEPT, not of the obstruction. `bases49.txt` remains a correct union of every obstructed verdict on
record — it is the interpretation as "the obstructed class" that is wrong.

⚠ **AND IT BREAKS THE `N` PATTERN.** All 49 have `N >= 3` (odd prime; no `N = 1`, no `N = 2`).
These two have **`N = 1`**. So "obstruction needs nontrivial level" was a SAMPLE ARTIFACT — the
same confound already flagged for `D`: within `sweep122` every base is either (even `D`, odd prime
`N`) or (odd `D`, `N = 2`), so `D`-parity and `N` were never separable there. Now they are, and the
`N` half is refuted.

⇒ The extent of the obstruction in `(D, N)` is genuinely unknown. Do not quote 49 as the size of
the class; quote it as "49 known".

### ⚠ AND THE (odd `D`, odd `N`) CELL IS UNREACHABLE BY CONSTRUCTION, not merely untested

For odd `D` we have `v_2(D) = 0`, so `M = 4DN/2^{v_2(D)} = 4DN = 2^2 · DN`. An odd `D` is a product
of an even number (>= 2) of odd primes, and an odd `N > 1` coprime to it adds at least one more, so
`DN` carries >= 3 distinct odd primes and `#div(DN) >= 8`. Hence

    #div(M) = 3 · #div(DN) >= 24     ALWAYS, for every (odd D, odd N>1)

All **74** such targets therefore sit behind the `#div >= 24` wall, which a direct Normaliz probe
confirms is real (`vvdata/weyl-campaign/normaliz-wall-probe.md`). ⇒ That is WHY `sweep122` contains
no such base — not an oversight in the sweep design. The `D`-parity confound in the obstructed
class **cannot be separated without first beating the wall**, so it is not a cheap experiment and
should not be listed as one.
