# Reproducing the classification: the proofs of the last 34 curves as pipeline stages

Status: design, not yet implemented (2026-11-05). Branch `classification-pipeline` = `main` plus every open
PR except the three based on the campaign branch (#132, #136, #137) and the draft tooling PR #47.

## The goal

The pipeline decides, for every Atkin–Lehner quotient of a Shimura curve X_0(D,N), whether it is
hyperelliptic. Its last committed run (PR #117) leaves 61 quotients undecided. Three further filters bring
these to 34, and each of the 34 has since been settled by a proof written for that curve: 33 are
hyperelliptic and one, X_0(141,2)/{1,6,94,141}, is not. Those proofs live in per-curve scripts, offline
tests and notes. The goal is a pipeline run that ends with **no undecided quotient**, where every decision
names the argument behind it, and where every argument is code that a reader can rerun.

## Where things stand

**The stages.** `workingcode.m:20-95` (`FILTER_STAGES`) and `run_pipeline.sh:14-55` list the same 42
steps (`tests/PipelineStages.m` checks they agree). After the dual graph at non-squarefree levels
(`FilterByDualGraph`, `DualGraph.m:679`), the last deciding stages are the twisted trace and Weil
polynomial filters, Noether's count on stored models (`FilterByCanonicalQuadrics`,
`CanonicalQuadrics.m:129`) and the non-Atkin–Lehner involutions. A decision is recorded on the curve
record (`alquots.m:3`): `IsSubhyp`/`IsHyp` are set, and `TestInWhichProved` names the stage; a curve is
undecided while `IsSubhyp` is unassigned.

**From 61 to 34 is not in the data.** `data/stage_counts.txt:40` still reads 61 undecided. The cut to 34
(22 by the metric dual graph of #116, 3 by special fibres, 2 by models at 6_35) is recorded only in
`HANDOFF_2026-10-19.md:320-328`, and the 34 are listed only in `paper/hyperelliptic-candidates.tex`
(table `tab:candidates`). The committed data also predates Criterion H (`FilterByComponentInvolution`).

**A rerun would already credit the stored models.** `FilterByCanonicalQuadrics` marks a quotient
hyperelliptic as soon as a stored model y² + hy = f of the right genus agrees with the trace formula at
three primes (`CanonicalQuadrics.m:100-127`). That includes the models of #128 that its own note calls
*candidate* (6_385, 714_5) or found in a *window* (15_38). Those curves are proved hyperelliptic by other
means (#151–#154), so the verdicts are right, but the label would name the wrong argument.

**The 34th curve has no code.** That X_0(141,2)/{1,6,94,141} is not hyperelliptic is asserted in two
sentences (`paper/heegner-height-models.tex:61-63`, `paper/hyperelliptic-candidates.tex:309-310`): the
enumeration of #128 finds one cover consistent with the trace formula, and it is a plane quartic. No test,
script or provenance entry reproduces this.

## The 34 curves by the shape of their proof

A label `D_N` names the curve X_0(D,N)/W with the W listed. "Fast" means the test runs in the normal suite
(minutes); "offline" means a test in `tests/_offline/` that needs a precomputed basis.

**A. Genus three, W the full group: enough degree-two CM divisors that are pairwise torsion-related.**
If there are more of them than the rational torsion of the Jacobian has elements, two are linearly
equivalent and give a g¹₂ (Theorem 1.2 of `paper/torsion-counting-hyperelliptic.tex`).

| curve | D·N | relation certified by | torsion bound | PR | test |
|---|---|---|---|---|---|
| 46_31 | 1426 | six banked Borcherds forms | 2 | #133 | fast |
| 6_379 | 2274 | dual Borcherds certificate at level 4548 | 1 | #146 | fast + offline (~1 h) |
| 21_58 | 1218 | dual test, 2 ∥ N variant, level 2436 | 2 | #147 | fast + offline (~40 min) |
| 77_10 | 770 | dual test, 2 ∥ N, searched | 2 | #148 | fast + offline (~25 min) |
| 39_34 | 1326 | dual test, 2 ∥ N, both dictionaries | 4 | #149 | fast + offline (~35 min) |
| 4130_1 | 4130 | dual test at level 8260 | 1 | #157 | fast + offline (basis 61 h) |
| 6_649 | 3894 | dual test at level 7788 | 6, lowered to 2 by Φ₃ | #156, #158 | fast + offline (basis 50 h) |

**B. Genus four over a genus-two star: the branch pair is a fibre of the star's g¹₂.** The quotient is a
double cover of the star branched at two points, and is hyperelliptic exactly when those two points form a
fibre (`genus4-over-genus2` normal form, #140/#142).

| curve | D·N | certified by | PR | test |
|---|---|---|---|---|
| 2002_1 | 2002 | dual certificate at level 4004 | #139, #140 | fast + offline (~1 h) |
| 2030_1 | 2030 | three banked Borcherds forms, level 4060 | #142 | fast |
| 2418_1 | 2418 | three banked Borcherds forms, level 4836 | #142 | fast |
| 6_403 | 2418 | three banked Borcherds forms, level 4836 | #142 | fast |
| 6_341 | 2046 | Borcherds forms, branch pair up to 2-torsion, then the star identified | #145 | fast |
| 3390_1 | 3390 | dual test at level 6780 | #155 | fast + offline (basis 43 h) |

**C. Double cover of an elliptic star: the branch pairing is torsion, then zero.** The dual test shows the
pairing is torsion; the reduction of the CM points to the Čerednik–Drinfeld fibre (Molina) and the
component groups Φ_p show it is zero; the Prym fixes the twist (`paper/component-groups.tex`, #158).

| curve | D·N | torsion step | zero step | PR | test |
|---|---|---|---|---|---|
| 15_38 | 570 | dual test, level 1140 | Φ₃ = Z/8, Φ₅ = Z/4 | #151 | fast + offline (~10 min) |
| 6_385 | 2310 | dual test, level 4620 | Φ₂ = Z/20 | #152 | fast + offline (basis 19 h) |
| 714_5 | 3570 | dual test, level 7140 | Φ₃ × Φ₇ × Φ₁₇ | #153 | fast + offline (basis 99 h) |
| 6_323 | 1938 | Gross–Zagier (Yuan–Zhang–Zhang) | reduction at 2 | #154 | fast |

**D. The sixteen models from Heegner heights (#128, note #141).** Each is a double cover of an elliptic
curve E branched at CM points whose heights the Gross–Zagier formula of Cai–Shu–Tian gives; the
branch divisors are enumerated, every double cover is built exactly, and the trace formula selects one.
The model is hyperelliptic, so the curve is.

1326_1, 210_13, 39_11, 55_7, 6_301, 6_455, 6_247, 6_253, 6_185, 6_187, 22_51, 141_2 ({1,3,94,282}),
21_34, 57_10, 15_46, 690_7. (#128 also models 15_38, 6_385, 714_5, 6_323, with caveats; group C proves
those four independently.) Tests: `tests/HeegnerHeights.m` reruns the route end to end only for
1326_1 and 210_13; `tests/ModelChecks.m` compares every stored model with the trace formula.

**E. The one that is not hyperelliptic.** X_0(141,2)/{1,6,94,141}: the same enumeration as D, read in the
other direction. Here completeness of the enumeration *is* the proof, so every input that could lose a
candidate (see below) must be justified. Nothing runs it today.

## The building blocks

Most of what the stage needs exists; the main gap is that the dual test is a script copied into each
proof rather than an intrinsic.

| block | today | to do |
|---|---|---|
| CM dictionary (Z(n) as a sum over d = −4n/r², weights 1/e) | inline `HasCM`/`Stab` in 9 tools and 11 tests | one intrinsic, with the 2 ∥ N plus-space and odd-D variants |
| dual test: a CM divisor certified torsion against W | `tools/proof-*/withbasis.m`, `stepH-*.m` | one intrinsic: input the curve, a divisor, a basis; output certified/refused, dim W, the Hecke check if determined |
| weight-3/2 bases | `tools/proof-*/basisS.m`, files of 70–330 MB | a cache keyed by (M, B) with checksums, stored outside git |
| torsion bound | `JacobianTorsionBound` (`JacobianTorsion.m:61`) | none |
| component groups, Molina | `ComponentGroupData`, `CMSpecialisation`, `ClassInComponentGroup`, `TorsionOrderBoundOfCMDivisor` (`ComponentGroups.m`) | none (check the squarefree-N statement in its header against its `require`) |
| Borcherds certificates | `WriteBorcherdsForms`, `BorcherdsFormsFromFile` (`BorcherdsForms.m:1981-1995`), files in `data/borcherds/`, `data/dual/` | none |
| exact covers, heights | `DoubleCoversWithBranchDivisor`, `MatchCoverToCounts` (`ExactDoubleCover.m`), `ModelByHeegnerHeights` (`HeegnerHeights.m:217`) | record the inputs of each run (curves of the isogeny class tried, windows) |
| Noether's count | `IsHyperellipticFunctionField` (`HyperellipticDetection.m:14`) | none |

Measured basis costs on lovelace (single thread, 100–125 GB): level 4620 19 h, 6780 43 h, 7788 50 h,
8260 61 h, 7140 99 h, roughly dim^4 in the dimension of S_{3/2}; precision is cheap and must reach every
target index (B ≥ |d|/4).

## The proposed stage

**Verify, do not search.** Each of the 34 gets a *certificate file*, `data/certificates/<label>.m`, read
the way banked Borcherds forms are read: it names the shape (A–E above) and holds what the proof needs —
the CM divisors or fibres, the banked forms or principal parts, the levels and precisions of the dual
tests, the primes for the component groups, and for D/E the run's inputs. A new stage,
`FilterByHyperellipticCertificates`, placed after `FilterByCanonicalQuadrics`, loads the certificate of each
undecided curve, rechecks it, and records the verdict with a label naming the shape and the certificate
(for example "Certificate A: 4 divisors against torsion 2, data/certificates/46_31.m"). The searches that
*found* the divisors (the Nelson screen, the scans in `tools/proof-77-10/`, the Borcherds searches) stay as
tools that write certificate files; the pipeline never depends on them.

**The bases.** A dual test needs the weight-3/2 basis of its level. The stage reads it from a cache
directory (`DUALW_CACHE`, default `~/dualw`), refuses a file whose checksum differs from the one in the
certificate, and leaves the curve undecided, with a message naming the missing basis, when it is absent.
Recomputing is opt-in. The bases are hosted outside git (Zenodo or a release asset), with their checksums
in the repository.

**Cost.** With the bases present, every check is minutes; the slowest fast tests today are 6_649 (~10 min)
and the component groups (~6 min). Without them, five levels cost about 270 CPU-hours.

**Registering a stage** touches: `FILTER_STAGES` (`workingcode.m`), `run_pipeline.sh` (a `run_seq` line and
an `UpdateCurvesAfter…` step), `run_sequential_stage.m`, `tests/PipelineStages.m`, `analysis_stages.m`,
`reconstruct_attribution.m`, `make_latex_tables.m`, `make_latex_figures.m` (a token group) and
`docs/RUNNING_PIPELINE.md`; data files are written with `SetColumns(0)`.

## Order of work

1. **Land 61 → 34 in the data.** Rerun the pipeline on this branch (Criterion H, the metric dual graph, the
   special fibres, the 6_35 models), commit the stage counts, and add `data/open_34.m` with a test that the
   undecided set equals it.
2. **The dual test as an intrinsic.** One CM dictionary and one dual test, with the 2 ∥ N and odd-D
   variants; the basis cache with checksums. Acceptance: it reproduces, from the saved bases, every
   certified relation and every refused control of #146–#157.
3. **The certificate format and the stage, one shape at a time:** A (46_31 first: no basis needed), then
   B, then C, then D. Each shape is accepted when the stage decides its curves from the certificates alone
   and the per-curve tests still pass.
4. **#128 made exact** (decisions 1 and 3): exact L-ratios and recognition, no windowed place, the oldform
   multiplier; first for the 34th, then for the sixteen curves of group D.
5. **Labels of stored models.** Decide how `FilterByCanonicalQuadrics` treats a stored model: credit it only
   when its provenance says the model is proved, or move the crediting of the #128 models into certificate D
   so that the label names the argument.
6. **The full run** on lovelace from this branch, ending with no undecided quotient; tables and figures
   regenerated.

## Decisions (Eran, 2026-11-05)

1. **The 34th is proved by making #128's enumeration exact for that curve.** L-values enter as exact
   ratios from modular symbols instead of zero-tests at 1e-8, heights and coordinates are recognised
   exactly instead of by `BestApproximation`, and the certificate shows that no branch place needs the
   window |n| ≤ 12. Then "one cover consistent with the trace formula, a plane quartic" is a proof, rerun
   by a test.
2. **6_323: the Gross–Zagier step is checked exactly.** The certificate computes the root numbers of
   E^(−8) and E^(969) and the rational number L(E^(969), 1)/Ω by modular symbols, and the stage checks it
   is 0; with the signs this gives order of vanishing ≥ 3. The theorem of Yuan–Zhang–Zhang is cited as an
   input, as Molina's theorems are for the component groups.
3. **The sixteen curves of group D: #128 made exact, the same way as for the 34th.** Point counts make the
   Jacobians isogenous, not the curves isomorphic, so the positive verdicts also rest on the enumeration
   missing nothing. Each certificate records its run's inputs (the curves of the isogeny class tried, the
   Mordell–Weil proofs, Manin constant 1) and the oldform multiplier is derived, or shown unneeded, at every
   place where it is used (15_46 at d = −15, p = 2 first).
4. **The bases go to Zenodo**, with checksums in the repository and a script that downloads them into the
   cache; the deposit is created from Eran's or Sachi's account.
5. **The list of 34 is a data file**, `data/open_34.m`, with a test that the curves left undecided before the
   new stage are exactly these; certificates are keyed by the same labels.

## Branch

`classification-pipeline` merges the open PRs as they stand on 2026-11-05. Conflicts were resolved by keeping
both sides in the provenance log, the paper README, `.gitignore` and the spec file; in `SchoferFormula.m`
the two fixes to the field of definition (no Artin map when the field is K; fall back to the slow routine
when the Artin map refuses the twist ideal) are combined, as #66 already does. When a PR here changes, merge
its branch again rather than editing it on this branch.

Checked on the merged branch (2026-11-05): 13 of 14 fast tests pass (CM fields of definition, the divisor
point at 21_2, dual graph, Kodaira, component groups, plane quartics, 15_38, 46_31, 2030_1, 6_341, 4130_1).
`tests/HyperellipticCRV.m` (#108) fails on a curve that #112 adds: X_0(38,7)/{1,7}, genus 7, stored as a
pair over a conic that has a rational point, so it has a model y^2 = f(x) and should be stored as one, as
#108 does at 10_53. Neither PR is wrong alone; #112's entry needs converting.
