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

    142_1  (D = 2*71,  N = 1)    158_1  (D = 2*79,  N = 1)
    166_1  (D = 2*83,  N = 1)    214_1  (D = 2*107, N = 1)
    6_97   (D = 6,     N = 97)

**FIVE new, so 54 known obstructed** (2026-09-13/14), all found incidentally by a routine backlog
sweep rather than by looking for them.

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


## ⚠⚠ AND A PATTERN I FLOATED FROM THESE FIVE IS REFUTED — by data already in `data/models/`

Four of the five are `D = 2p`, `N = 1` with `p >= 71`, which looked like a size effect. **It is not.**
Restricting to `D = 2p`, `N = 1`:

    HAVE MODELS   p = 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 73 89 97 101 103   (23)
    OBSTRUCTED    p = 71 79 83 107                                                          (4)

`71` is obstructed but `73` builds; `79` and `83` are obstructed but `89 97 101 103` all build;
then `107` is obstructed. **Interleaved, not a threshold.** Several of the builders are even
EXTERNALLY corroborated — `14_1`, `34_1`, `46_1` against Gonzalez-Rotger, `82_1` against Guo-Yang —
so `D = 2p, N = 1` is one of the BEST-reproduced shapes in the repo, not a suspect one.

This is the same conclusion the record already held for a fixed `D` (`38_7` is larger than `38_5`
and fully surjective); it now holds across `D` too. ⇒ **The obstruction is not predictable from the
shape of `(D, N)`.** Every `D = 2p` family contains both: `6_5 ... 6_83` build while `6_97`/`6_103`
are obstructed; `10_3 ... 10_61` build while `10_43`/`10_47`/`10_59` are obstructed.

⇒ That is a direct argument for the **DEFICIT PREDICTOR** (`vvdata/weyl-campaign/deficit.m`): the
deficit is `Ncols(mat) - Rank(ech_basis * mat)`, independent of any divisor choice, so obstruction
is computable from the Weil representation WITHOUT the CM points or the 96-triple search. Five new
obstructed bases found by accident, with no predicate in sight, is evidence that guessing from
`(D, N)` will not work.

⚠ `deficit.m` is **EVEN `D` ONLY** (for odd `D` the 0-side block is missing and the deficit drifts
with pole order instead of staying invariant). All five new bases have even `D`, so they are a
VALIDATION SET the script was never calibrated on — its recorded ground truth is only
`38_5 -> 1`, `38_7 -> 0`, `34_3 -> 0`.
