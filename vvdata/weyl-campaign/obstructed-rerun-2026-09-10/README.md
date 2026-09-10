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
