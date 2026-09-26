# Running the classification pipeline

This page covers a full rerun of the hyperellipticity pipeline in parallel, on a Linux server.
`run_pipeline.sh` runs the stages of `FILTER_STAGES` in `workingcode.m`, plus a Weil-polynomial
pass on the star curves (`FilterByWeilPolynomialStar`), and splits each `FilterBy*` stage across
many Magma processes. The sequential alternative is
`GetHyperellipticCandidates(:recompute_data)`; it takes about four days on one core.

## Requirements

* **Magma**, on the `PATH` as `magma` (the wrapper that sets `MAGMAPASSFILE`).
* **GNU parallel**. `run_parallel_filter.sh` uses it to run the chunks.
* **GNU `timeout`** (coreutils). The scripts do not set timeouts themselves. Use `timeout` for the
  optional per-job cap below, and to bound any hand-run stage. macOS has no `timeout`, which is
  one reason to run this on Linux.
* **`NORMALIZ_BIN`**. No filter stage runs a polytope solve, but the repo rule is to set it
  whenever the package is loaded, because a missing value fails silently:

      export NORMALIZ_BIN=/path/to/normaliz-3.11.1/normaliz

* **`CLASS_GROUPS_FAST_DIR`**. This is the class-number cache that the trace-formula stages use.
  It can grow to hundreds of GB. Point it at a large local disk. The default in
  `run_pipeline.sh` is `/var/tmp/class-groups-fast`.
* **Memory.** Most stages are light. `FilterByTwistedTrace` and `FilterByTwistedWeilPolynomial`
  hold a modular symbols space of level D·N and its Hecke operators, which reaches several GB
  per worker at the largest levels (D·N in the thousands); see "Cost hot spots".

## The full parallel rerun

From the repo root, on a clean checkout of the branch you want to measure, **into an empty data
directory**:

    export NORMALIZ_BIN=/path/to/normaliz
    export CLASS_GROUPS_FAST_DIR=/big/disk/class-groups-fast
    nohup ./run_pipeline.sh 64 data/par 1024 < /dev/null > pipeline.log 2>&1 &

The arguments are `[num_workers] [data_dir] [num_chunks]`:

* `num_workers` (default 128) is the number of Magma processes that run at once. Set it to the
  number of cores you may use. If memory is tight, lower it for the two twisted stages.
* `data_dir` (default `data/par`) receives every intermediate `curves_after_<stage>.dat`. The
  committed `data/` directory is never written. **For a rerun, start from an empty directory:
  wipe an existing `data/par` (or pass a fresh name such as `data/par-2026-10`).** A stage is
  skipped whenever its output file exists, so a stale later stage left over from an older run
  would be reused as it is, and silently drop everything the new stages decide before it.
* `num_chunks` (default 1024) is the number of pieces each parallel stage is split into. Many more
  chunks than workers lets the cost-aware ordering isolate the few heavy curves.

`MAGMA_CMD` overrides the Magma command (default `magma -b`). To put a per-job cap on the parallel
stages:

    MAGMA_CMD="timeout 72h magma -b" nohup ./run_pipeline.sh 64 data/par 1024 < /dev/null > pipeline.log 2>&1 &

A chunk that times out fails its stage (`parallel --halt now,fail=1`), and no merged file is
written. That makes it visible; it never produces silently incomplete output.

**Restarting.** Every stage is skipped if its output file already exists. To resume after a crash
of *this* run, rerun the same command. To redo a stage, delete its `curves_after_<stage>.dat` and
every later one; never mix files from different runs or code versions in one directory.
Per-chunk logs are in `data/par/parallel_chunks_<stage>/results/`, and the job log is
`joblog.txt` there.

**One stage by hand:**

    ./run_parallel_filter.sh FilterByTwistedTrace 64 data/par 1024 < /dev/null
    magma -b stage:=UpdateCurvesAfterTwistedTrace input_dat:=data/par/curves_after_FilterByTwistedTrace.dat \
          output_dat:=data/par/curves_after_UpdateCurvesAfterTwistedTrace.dat run_sequential_stage.m < /dev/null

**A stage that stops on purpose.** The twisted stages raise an error, and so stop their stage,
if the W-fixed D-new modular symbols of a curve do not have dimension 2g (`BADDIM` in the
message). That means the genus or W of the curve is inconsistent with the space, and no verdict
of the filter could be trusted; it was never met in the sweeps. Investigate the named curve
rather than skipping it.

**When it finishes,** the final file is `data/par/curves_after_UpdateCurves8.dat`. Compare it with
the committed data before promoting it, for example by counting verdicts and listing every curve
whose verdict changed. Then copy `data/par/curves_after_*.dat` into `data/`. Tests and
`GetHyperellipticCandidates()` read the files in `data/`.

## Stage order

The star-curve block in `run_pipeline.sh` reads:

    FindPairs, UpdateGenera, UpdateByGenusStar
    FilterByTraceStar                         parallel
    FilterByTwistedTraceStar                  parallel, by level  NEW
    HHProposition1                            (+ VerifyHHTable2, VerifyHHProposition1)
    SpecialFiberIsomorphismStar
    FilterByWeilPolynomialStar                parallel
    FilterByTwistedWeilPolynomialStar         parallel, by level  NEW
    FilterStarCurvesByFpAutomorphisms         parallel
    FilterByNonALInvolutionsStar              parallel
    GetQuotientsAndGenera + UpdateByGenus     (+ VerifyFHTheorem3)

`FILTER_STAGES` in `workingcode.m` has no star Weil-polynomial stage, so there the star twisted
Weil stage comes directly after the star twisted trace. On a star curve W is the full AL group,
so the only twists are V2 (8 | N), V3 (9 || N, and then 9 is in W) and V2 V3. A star level with
neither has nothing to test and its modular symbols are never built. As with
`FilterByNonALInvolutionsStar`, the star determinations are carried onto the full-W entries by
`GetQuotientsAndGenera`, and from there to the covers by the closures. The group filter does not
run on the star curves.

The star twisted stages decide some D = 1 curves of [HH] Table 2 (for example X_0^*(396), by V3
at q = 5). So `VerifyHHTable2` and `VerifyHHProposition1` are now run on HH's own input, the
`FilterByTraceStar` snapshot with HH Proposition 1 applied, and no longer on the pipeline state.
The check is the same, exact reproduction of the table.

The all-quotients block (after `GetQuotientsAndGenera`) now reads:

    UpdateCurves1
    FilterBySpecialFiber                      parallel
    FilterByALFixedPointsOnQuotient           parallel
    UpdateCurves2, Genus3CoversGenus2, UpdateCurves3
    FilterByDegeneracyMorphism                parallel
    UpdateCurves4
    FilterByComplicatedALFixedPointsOnQuotient    parallel
    FilterByGeneralizedComplicatedFixedPoints     parallel
    UpdateCurves5                             (+ VerifyFHTable3)
    FilterByAutomorphismGroup                 parallel            NEW
    UpdateCurvesAfterAutomorphismGroup                            NEW
    FilterByTrace                             parallel
    UpdateCurves6
    FilterByTwistedTrace                      parallel, by level  NEW
    UpdateCurvesAfterTwistedTrace                                 NEW
    FilterByWeilPolynomial                    parallel
    UpdateCurves7
    FilterByTwistedWeilPolynomial             parallel, by level  NEW
    UpdateCurvesAfterTwistedWeilPolynomial                        NEW
    FilterByNonALInvolutions                  parallel
    UpdateCurves8

The new closure stages are called `UpdateCurvesAfter<Stage>` rather than renumbering the existing
`UpdateCurves6..8`. Stage names are also file names, and tests and `GetHyperellipticCandidates()`
read `curves_after_UpdateCurves<N>.dat`, so every existing name keeps its meaning.

**Splitting by level.** Most parallel stages deal curves to chunks one at a time, in descending order of
`CurveCostProxy`. The twisted stages (all four) compute the modular symbols of level D·N once and use them
for every curve at that level. `parallel_filter_worker.m` therefore deals whole **levels** to the
chunks for them, again heaviest first. The merge (`parallel_merge.m`) is unchanged: every chunk's
curves are tagged with their original index.

## Cost hot spots

* **`FilterByWeilPolynomial`**: about 5 h on a single curve at the top end, dominated by
  class-number lookups at depth 4·Qmax·p^g. The heavy curves are dispatched first. The makespan
  of this stage is roughly the slowest single curve.
* **`FilterByTwistedTrace`**: modular symbols of level D·N, plus T_p for every good p < 4g².
  Small levels take seconds (level 1530 takes about 30 s), but the largest levels (D·N from about
  2000 up to 15330) take **hours each**, and the sweep's top levels ran for more than 4 h.
  Expect the stage's makespan to be set by the largest level, not by the total. There is **no
  level cap**: every level is run. The star version is lighter, since only star levels with V2
  or V3 (8 | N or 9 || N) do any work.
* **`FilterByTwistedWeilPolynomial`**: the same modular symbols, but only T_p at the table primes
  (g = 3: p ≤ 23; g = 4: p ≤ 5; g = 5, 6: p = 2), and only for g ≤ 6. It is cheaper than the
  twisted trace, but still hours at the largest levels. No level cap either.
* **`FilterByAutomorphismGroup`**: cheap coset arithmetic, except when a subgroup would violate the
  lemma and contains a non-AL central involution. Its quotient genus then comes from
  `TraceDNewQuotient`, which takes about 30 s per involution at level ~1000 and grows with the
  level.
* **`FilterByNonALInvolutions`** and **`FilterByGeneralizedComplicatedFixedPoints`** skip levels
  above 3000 (`NonALModSymMaxLevel`, `GeneralizedComplicatedMaxLevel`).

## After the rerun

`analysis_stages.m`, `make_latex_tables.m` and `reconstruct_attribution.m` (and so
`make_pipeline_summary.sh`) list the new stages in pipeline order. A stage whose file does not
exist yet is skipped with a note, so on data from before the new stages they reproduce the old
outputs exactly. After the rerun, regenerate them and update the paper's counts;
`make_latex_tables.m` cites `sec:autgroups`, a label the paper has to define.
