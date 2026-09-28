# Handoff — 2026-09-28 — open the PR for the pipeline data rerun

**For an agent with GitHub access (`gh` authenticated).** lovelace's `gh` is not logged in, so the
branch was pushed from there but the PR could not be opened.

## What exists

- Branch **`data/pipeline-rerun-2026-09-26`** on `origin`, head **`98c33b8d`**, two commits on top
  of `origin/main` (`8af0c552` at the time):
  1. `d7dc64da` — every `curves_after_<Stage>.dat` from the full `run_pipeline.sh` rerun on lovelace
     (2026-09-26 to 09-28), exactly as produced (from `data/par-2026-09-26c` there).
  2. `98c33b8d` — the 4 HHProposition1 labels replaced by SpecialFiberIsomorphism (PR #48 made
     HHProposition1 check-only; same verdicts, same source curves).
- The commit messages carry the full provenance and checks; the PR body below summarises them.

## Task

1. Find the open issue that tracks the missing pipeline data, e.g.
   `gh issue list --repo assaferan/ShimuraCurveALQuotients --state open --search "data"` (also try
   "missing", "curves_after"). Call its number `N`. **If you cannot identify it with confidence, ask
   the user — do not guess.**
2. `git fetch origin` and check `git log --oneline origin/main..origin/data/pipeline-rerun-2026-09-26`
   shows exactly `98c33b8d` and `d7dc64da`. If `main` has moved on, do **not** rebase or force-push;
   just open the PR.
3. Save the body below (between the two `----8<----` lines) as `pr_body.md`, with `#<issue>`
   replaced by `#N`, then:

       gh pr create --repo assaferan/ShimuraCurveALQuotients \
         --base main --head data/pipeline-rerun-2026-09-26 \
         --title "Data: full pipeline rerun 2026-09-26/28, with the HHProposition1 relabel" \
         --body-file pr_body.md

4. Do **not** merge it, and do not edit, rebase or push the branch. Report the PR URL.

## PR body

----8<----
## Summary

This PR commits the output of a full `run_pipeline.sh` rerun on lovelace (2026-09-26 to 09-28), the pipeline data that was missing. It closes #<issue>.

Two commits:

1. **`d7dc64da`: the files exactly as produced.** Every `curves_after_<Stage>.dat` from the run: 22 replace the old files, 9 are new stage files (the automorphism-group and twisted filters and their updates), and 4 are byte-identical to the old ones. `curves_after_D1Oracle.dat` is not a pipeline output and is unchanged.
2. **`98c33b8d`: HHProposition1 relabel.** Since #48 HHProposition1 is check-only. Its 4 star verdicts are now made by SpecialFiberIsomorphismStar: **same verdict, same source curve**, only the label changes (X_0^*(194), X_0^*(546), X_0^*(205,3), X_0^*(1995,2)). I checked this by rerunning the star stages up to SpecialFiberIsomorphismStar with `main`'s code, and by loading every file in Magma to confirm that only those 4 labels differ.

## Result (`curves_after_UpdateCurves8.dat`, 18379 curves)

| | this run | before |
|---|---|---|
| non-subhyperelliptic | 11884 | 11702 |
| subhyperelliptic | 5847 | 5847 |
| undecided | **648** | 830 |

- 0 contradictions with the old file and 0 old verdicts lost.
- 182 curves newly decided: twisted trace 68, twisted Weil about 38, automorphism group 8, special fiber the rest.
- 0 failed chunks in any parallel stage.

## Provenance

The code was `integration` at `7b9b69e8`: #48 through `34a3df9e`/`47b63bb4`, plus the merged fix/gc-v3-certificates, fix/weil-poly-at-2 and fix/nonal-involution-guards. That is **not** exactly `main`. Since then `main` has gained #48's follow-ups (HHProposition1 check-only, which the second commit covers, and `1c193b09`) and `bd7f08fc` / `ddb01bc9` / `e363db4b` in the NonAL and GeneralizedComplicated filters. **Whether those change any output has not been checked yet.**

## Timing

About 36 h wall time from the relaunch, on 128 workers:
- FilterByWeilPolynomial 16.0 h, almost all of it on one (21,20) g=5 curve
- FilterByWeilPolynomialStar 8.7 h
- FilterByTwistedTrace 4.6 h
- FilterByTwistedTraceStar 3.8 h

#50 (skip star curves in the non-star stages) removes the re-tests of star curves inside those stages.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
----8<----

## Still open (not part of the PR task)

- Check `main`'s NonAL / GeneralizedComplicated changes (`bd7f08fc`, `ddb01bc9`, `e363db4b`) against
  this data: rerun those stages with `main`'s code on the run's inputs and compare. Not done yet.
- `docs/RUNNING_PIPELINE.md:35` example still points at `data/par`, which holds 12 committed
  pre-PR snapshots; `run_pipeline.sh` skips existing outputs, so it silently reproduces old results.
