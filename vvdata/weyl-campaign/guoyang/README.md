# Guo-Yang CM-value tables: extraction and verification record

Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), has 42 appendix
tables titled "CM-values of X_0^D(N)" -- Schofer's formula applied to many bases, tabulating
the primary hauptmodule (and, for level > 1, one or two further hauptmoduls of partial
Atkin-Lehner quotients) at CM points. Before this pass, exactly ONE of the 42
(`X_0^15(2)`, "Table 45") had been turned into a test (`tests/SchoferIsometry.m`), by hand,
with a manually-derived normalisation constant. This directory holds the machinery that
turned the rest into tests.

## Files here

    ShimuraCurves-arxiv.tex   the arXiv v1 TeX source (`curl -sL arxiv.org/e-print/1510.06193`,
                               despite the .tar.gz-looking URL it's a single gzipped .tex file)
    parse_tables.py           extracts all 42 tables into tables.json
    tables.json               structured data: per table, D, N, column headers/relations, and
                               per row the discriminant + each column's value (kind: rat/oo/other)
    gen_drivers.py            emits throwaway Magma verification drivers (used to MEASURE which
                               bases pass before committing anything as a permanent test)
    gen_final_tests.py        emits the permanent tests/_offline/GuoYang_D_N.m files for the
                               bases that passed

## Method

Only the PRIMARY hauptmodule column (full star quotient X_0(D,N)/W_{D,N}) is covered so far --
`ncols` in tables.json shows some tables have 1-2 more columns for partial AL-quotients ("other
functions"), left for a follow-up (would need selecting the matching intermediate ShimuraQuot
object from `curves` by its W-subgroup, and most of those secondary columns are radical-valued,
not plain rationals, which the current parser only records as raw LaTeX strings).

The check is a cross-ratio (Mobius-invariant) comparison, not a hand-derived normalisation
constant: pick 3 discriminants the table lists with distinct values, use OUR computed values
there as one Mobius frame (z0,z1,z2 -> 0,inf,1) and Guo-Yang's published values there as
another, then compare cross-ratios at every other listed discriminant. This needs no assumption
that our hauptmodule shares Guo-Yang's zero/pole choice -- see `tests/_offline/GuoYangCheck.m`
on tier1-models for the implementation (generalises `tests/ExternalCMValues.m`'s single-value
version of the same idea).

Run on `lava` (32 cores, jump through `lovelace`) at 16-way parallelism after `lovelace` itself
turned out to be at load 320/256 from other users -- NOT idle, contrary to the "idle" note this
memory carried from 2026-09-02. `lava` needed a fresh clone (`git clone --branch main --depth 1
https://github.com/assaferan/ShimuraCurveALQuotients.git`); no prior checkout existed there.

## Result, 2026-09-03

    all 42 bases attempted

**34 PASS clean** (254 total individual CM-value checks, 0 failures):
`10_11 10_13 10_23 134_1 14_3 14_5 146_1 15_2 194_1 206_1 21_2 22_3 22_5 26_1 35_1 38_1 39_1
39_2 51_1 55_1 57_1 58_1 6_11 6_17 6_19 6_29 6_31 6_37 62_1 74_1 82_1 86_1 87_1 94_1`
-- now permanent tests, `tests/_offline/GuoYang_D_N.m` on tier1-models. `15_2` reproduces 18 of
27 rows exactly (9 rows dropped: discs where `ValuesAtCMPoints`'s `Keep` didn't admit them --
not investigated, they may not satisfy this order's local splitting condition), superseding the
hand-normalised subset in `SchoferIsometry.m` (keep that test too -- it also checks isometry
invariance and the cusp-support invariant, which this one doesn't).

**26_3: 2 of 11 FAILED, and it is NOT noise.** Discs -267 and -708 both fail with `got` and
`want` related by EXACTLY the Mobius involution `z -> z/(z-1)` (25/8 <-> 25/17, 49/11 <-> 49/38,
both verified to the last digit). That map fixes 0 and 2, swaps 1 and infinity -- the signature
of a two-way CM-point ambiguity (our pipeline and Guo-Yang consistently picking different
members of a conjugate pair at exactly these two discriminants), not an arithmetic error. Worth
a real look -- NOT added as a test, and NOT yet understood. Driver output preserved:
`gy_verify_26_3.log` alongside this file on lava (`~/gy_verify/logs/verify_26_3.log`) --
pull it before it's cleaned up.

**10_19: EXCLUDED, arXiv v1's own appendix table looks stale relative to the paper's worked
Example 37.** `tests/ExternalCMValues.m` already has a PASSING check on this exact base using
reference discs -8 (0), -40 (infinity), -3 (1) from Guo-Yang's Section 4.2 prose, correctly
getting s(-760) = 32/5. But the arXiv v1 appendix table parsed here has NO -40 row at all --
instead it lists disc -52 TWICE with two different values (infinity and -4), which this
extraction correctly refuses to use (ambiguous). Re-running the cross-ratio check with the
table's own reference triple (-3, -8, -67 after the -52 exclusion) fails uniformly at every
other row by an exact constant factor (5/23, checked to the last digit across all 6 remaining
rows) -- the signature of a reference-point/table mismatch, not a per-discriminant bug. Given
the codebase's own independent check on the SAME base already passes with data this table
doesn't even contain, the conclusion is that this specific appendix table, in the v1 preprint,
predates whatever correction produced the published Example 37 (arXiv has only one version,
v1 -- no way to fetch a newer one from arXiv itself). Not added as a test. If revisited: track
down the actual published (Compositio/Proc. LMS, doi:10.1112/S0010437X16007739) PDF's appendix
table for X_0^10(19) rather than trusting this extraction.

**4 bases hit known pipeline boundaries, not table-check failures:**
- `69_1`: `Runtime error in 'RationalNumber': s does not represent a rational number!` -- the
  non-rational-value class (memory: [[embedding-selection-root-cause]], 15_2's own issue).
- `93_1`, `95_1`, `159_1`: `Assertion failed` -- matches the documented squarefree-N method
  boundary (memory: [[assert-failed-is-squarefree-n]]), not investigated further here.

**`111_1` and `119_1`: KILLED at ~17.5 h CPU each, stuck building the basis.** Neither ever
returned from `ValuesAtCMPoints` (the driver's first `printf` comes after it, and both logs stayed
empty). Diagnosis, from `/proc` on lava rather than from the logs:

* **not** a polytope solve -- no `normaliz` child process on either, and the `M = 444` (111_1) and
  `M = 476` (119_1) solutions are already in the committed cache;
* `119_1` reached **VmPeak 40.6 GB** (8.2 GB resident at kill, 572 s system time -- heavy
  allocation churn), far past the "Magma dies ~11 GB on pools >~2000 vectors" threshold;
  `111_1` peaked at 3.4 GB;
* both are **odd D at level 1** (`111 = 3*37`, `119 = 7*17`) -- exactly the class whose documented
  remaining ceiling is `basis_of_weakly_holomorphic_forms(... : Zero)`, steep in pole order
  (556 s at pole order 845; these are far higher).

⚠ Localised to "inside `BorcherdsForms`, not the polytope stage" -- NOT proven to be the
`Zero`-side basis specifically. A verbose re-run (`SetVerbose` + `WriteStderr`) would pin it.
Note lava runs `earlyoom`, so a 40 GB excursion is also a candidate for being killed from outside.

Neither is in the sweep122 or wave-4 triage records (only in `MISSING_TARGETS.txt`), so these two
were genuinely uncharacterized before this. Of the **10** level-1 bases with no `X0_D_1.m`
equations test (`51 55 57 69 87 93 95 111 119 159`), this batch splits **4 pass / 4 error /
2 basis-ceiling** -- the trouble concentrates exactly in the previously-untested ones.

**Timing, for planning future batches:** wildly variable, 26s (`38_1`) to ~7000s (`87_1`); the
heaviest 6 bases (`6_29`, `6_31`, `6_37`, `10_19`, `10_23`, `87_1`) account for the bulk of the
wall-clock, and `111_1`/`119_1` never finished at all. 16-way parallel on 32 cores finished the
other 34 well within the time the stragglers took.

## Not done here (follow-ups)

* The secondary "other functions" hauptmodule columns (partial AL-quotients) -- data already
  extracted into `tables.json` (raw LaTeX strings, many radical-valued), needs a curve-selection
  step (matching each column's W-subgroup to the right `ShimuraQuot` in `curves`) and, for
  radical entries, a symbolic (not just rational) comparison.
* `26_3`'s involution anomaly -- real, unexplained.
* `10_19` -- needs the actual published table, not the arXiv v1 draft.
* The 4 pipeline-boundary bases (`69_1`, `93_1`, `95_1`, `159_1`) -- already have named root
  causes elsewhere in memory; not re-diagnosed here.
* `111_1`, `119_1` -- killed at the basis stage; a verbose re-run would confirm which basis call
  they sit in, and whether the odd-D `Zero`-side ceiling is the whole story.
