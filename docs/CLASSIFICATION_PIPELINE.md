# Reproducing the classification: the proofs of the last 34 curves as pipeline stages

Status (2026-11-09): implemented. The certificate stage decides all 34 and the five further quotients; the
pipeline ends with no quotient undecided. Branch `classification-pipeline`
= `main` plus every open PR except the three based on the campaign branch (#132, #136, #137) and the draft
tooling PR #47.

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
polynomial filters, Noether's count on proved stored models (`FilterByCanonicalQuadrics`,
`CanonicalQuadrics.m:129`), the non-Atkin–Lehner involutions and, last, the certificate stage
(`FilterByHyperellipticCertificates`, `HyperellipticCertificates.m`). A decision is recorded on the curve
record (`alquots.m:3`): `IsSubhyp`/`IsHyp` are set, and `TestInWhichProved` names the stage; a curve is
undecided while `IsSubhyp` is unassigned.

**The data (rerun of 2026-11-06 on this branch).** Before the certificate stage 39 quotients are undecided
(`data/stage_counts.txt`): the 34 of `paper/hyperelliptic-candidates.tex`, listed in `data/open_34.m`, and
five more that the second list of that file names with what each is missing (`tests/OpenThirtyFour.m`
checks both lists against the data). Of the 27 that the earlier notes settle beyond the 61→34 cut, the
stages now decide 22; the five are:

* the two quotients of X_0(6,35): their stored models are not hyperelliptic (Noether's count), but the
  models rest on the m = 0 term of Schofer's formula at a composite level. Both are proved not hyperelliptic
  without them, by shape V (2026-11-09, `paper/classification-certificates.tex`, Section 7): each is a
  Klein-four cover of the genus-zero star with three elliptic intermediate quotients, and over
  E = X/{1,3,5,7,15,21,35,105} its branch places Z_-420, Z_-84 have different fields, so they are fibres over
  the star and a function of degree two would give the star a genus-zero double cover inside the curve.
  X_0(6,35)/{1,3,5,15} has a second proof, shape T (the dual graph at 3);
* X_0(14,51), X_0(38,15), X_0(158,3) with the groups in the file: at p = 3 the special fibre is the genus-3
  curve Z = X_0(D, N/3)/W′, which the dual graph at another prime shows is not hyperelliptic over Q. That
  does not rule X out: X has good reduction at 3 with X mod 3 ≅ Z mod 3, and Z mod 3 is hyperelliptic
  (2026-11-08, `tools/special-fibre-3/zmod3.m`: for 14_17 and 38_5 the certified Heegner heights leave one
  cover of the star; it reduces mod 3 to a smooth hyperelliptic curve with the trace formula's L-polynomial
  of Z at 3). All three are in fact **hyperelliptic**, by the certificate stage: 14_51 and 38_15 by shape C
  (dual tests at 1428 and 1140), 158_3 by shape O (`paper/classification-certificates.tex`, Section 6).

**Status (2026-11-09): all 34 are certified, and the five further quotients.** The certificate stage
decides every one of the 34 (33 hyperelliptic, the 34th not), the three special-fibre quotients above
(hyperelliptic) and the two quotients of X_0(6,35) (not), so the pipeline ends with no quotient undecided
(`data/stage_counts.txt`: 12,493 not hyperelliptic, 5,886 hyperelliptic, 0 undecided;
`tests/OpenThirtyFour.m`). The 39 curves by the shape of the certificate the stage uses: A 7, B 5, B2 1
(6_341), C 15 (three of group C, ten of group D, 14_51, 38_15), C4 1 (6_323), D 6 (the rest of group D), N 1
(the 34th), O 1 (158_3), V 2 (the 6_35 pair); 6_35/{1,3,5,15} also has a shape-T certificate, kept as a second
proof. See `paper/classification-certificates.tex`.

**What is left.** No quotient. The open item is the Zenodo deposit of the bases (decision 4). The full run from stage 1 reproduces the committed data (order of work,
item 6).

**Stored models decide a curve only when proved.** `FilterByCanonicalQuadrics` reads
`data/models/proved_models.m` and uses a stored model at D > 1 only for a curve listed there with the
argument that proves its equation. Agreement of a model with the trace formula does not prove it. The list
has one entry, X_0(390,1)/{1,6,26,39}: built by Borcherds products at N = 1
(no m = 0 term), every sign of s, 1 − s and y² at the CM points used is forced and the cover equations
are unique up to the twist (`tests/_offline/ProvedModel390_1.m`, ~15 min). The 25 of the 34 that stored
models used to decide (the sixteen of group D, 15_38, 6_385, 714_5, 6_323, and the genus-four models) are
now decided by their certificates; their models stay in the data as equations.

**The genus-0 data of the special-fibre stages.** D = 1 uses Magma's classical modular curves and D = 6
Elkies' hypergeometric description; D = 10 and 22, and the registry of level-M bases (X_0(6,M) for
M = 5, 7, 11, 13, 17, 19 and X_0(14,11)), use our own Borcherds/Schofer values. Measured 2026-11-06: with
those three tests switched off the pipeline ends at the same undecided set, so no verdict depends on them;
and the current code reproduces the stored data of all seven registry bases (values after the Möbius
normalisation, conics up to a square, no ambiguous sign), except 2–8 discriminants per base that the
recomputed table did not reach.

**The 34th curve.** That X_0(141,2)/{1,6,94,141} is not hyperelliptic was first asserted from #128's
enumeration (one cover consistent with the trace formula, a plane quartic;
`paper/heegner-height-models.tex:61-63`, `paper/hyperelliptic-candidates.tex:309-310`). It is now proved by
shape N (section E below), with no enumeration.

## The 39 certified curves by the shape of their proof

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
| 6_341 | 2046 | Borcherds forms, branch pair up to 2-torsion; shape B2: the two curves of the isogeny class (#145, with CHIMP's certificate for C_A) both exclude K + t | #145 | fast (`tests/HyperellipticCertificatesB2.m`, ~2 min) |
| 3390_1 | 3390 | dual test at level 6780 | #155 | fast + offline (basis 43 h) |

**C. Double cover of an elliptic star: the branch pairing is torsion, then zero.** The dual test shows the
pairing is torsion; the reduction of the CM points to the Čerednik–Drinfeld fibre (Molina) and the
component groups Φ_p show it is zero; the Prym fixes the twist (`paper/component-groups.tex`, #158).

| curve | D·N | torsion step | zero step | PR | test |
|---|---|---|---|---|---|
| 15_38 | 570 | dual test, level 1140 | Φ₃ = Z/8, Φ₅ = Z/4 | #151 | fast + offline (~10 min) |
| 6_385 | 2310 | dual test, level 4620 | Φ₂ = Z/20 | #152 | fast + offline (basis 19 h) |
| 714_5 | 3570 | dual test, level 7140 | Φ₃ × Φ₇ × Φ₁₇ | #153 | fast + offline (basis 99 h) |
| 6_323 | 1938 | Gross–Zagier (Yuan–Zhang–Zhang) | reduction at 2 | #154 | fast + offline (shape C4, ~10 min, no basis) |

**D. The sixteen curves modelled from Heegner heights (#128, note #141).** Each is a double cover of an
elliptic curve E branched at CM points whose heights the Gross–Zagier formula of Cai–Shu–Tian gives; #128
enumerates the branch divisors, builds every double cover exactly and lets the trace formula select one.
That gives equations, not proofs; the certificates below prove the curves hyperelliptic.

1326_1, 210_13, 39_11, 55_7, 6_301, 6_455, 6_247, 6_253, 6_185, 6_187, 22_51, 141_2 ({1,3,94,282}),
21_34, 57_10, 15_46, 690_7. (#128 also models 15_38, 6_385, 714_5, 6_323, with caveats; group C proves
those four independently.) Tests: `tests/HeegnerHeights.m` reruns the route end to end only for
1326_1 and 210_13; `tests/ModelChecks.m` compares every stored model with the trace formula.

Certified (2026-11-08): ten by shape C (39_11, 6_185, 141_2, 21_34, 55_7, 57_10, 15_46, 6_253, 6_247,
6_301), the other six by shape D, #128's route with certified L-values (`CertifiedLValueAtOne`,
`CertifiedLDerivativeAtOne`), k² recognised within a proved error and b with denominator dividing
φ(D)ψ(N) (Cai–Shu–Tian's Hodge class): 6_187, 22_51, 6_455, 690_7, and 1326_1, 210_13 over an intermediate
quotient X/W′ (`cert["Wp"]`), identified by its reduction types I_n. `tests/_offline/HyperellipticCertificatesD.m`
reruns 6_187 with a refused control (~3 min).

**E. The one that is not hyperelliptic.** X_0(141,2)/{1,6,94,141}: decided by shape N (one component group,
`paper/classification-certificates.tex`, Section 2), with no enumeration. The only Galois-stable pairing of
its branch places, {P_-4, P_-12} | Z_-1128, has class 2 in Φ₃ = Z/7.

**F. The three quotients with a special fibre at 3** (outside the 34; `paper/classification-certificates.tex`,
Section 6). Each has good reduction at 3 with X mod 3 ≅ Z mod 3 for Z = X_0(D, N/3)/W′. Z mod 3 is
hyperelliptic, computed for 14_51 and 38_15 (`tools/special-fibre-3/zmod3.m`) and a consequence of the
certificate for 158_3, so the special fibre gives no obstruction.

| curve | star | certified by | test |
|---|---|---|---|
| 14_51 | 238a1 | shape C: dual test at 1428, Φ₂ = Z/14, Φ₇ = Z/2, simple Prym | offline (`HyperellipticCertificatesSpecialFibre.m`, ~3 min with 38_15) |
| 38_15 | 190a1 | shape C: dual test at 1140, Φ₂, Φ₁₉ | offline (same file) |
| 158_3 | 158a1 | shape O: one degree-4 branch place (d = −632, D₄ with k = Q(√−79)); the form is old at 3, so the χ-component is multiplied by 1 + χ(𝔭) = 0 at the split prime 3; E^(−79)(Q) is torsion-free and E(Q)[2] = 0 | fast (`HyperellipticCertificatesO.m`) |

**G. The two quotients of X_0(6,35), not hyperelliptic** (outside the 34; `paper/classification-certificates.tex`,
Section 7; `data/certificates/6_35.m`). Shape V (`tests/HyperellipticCertificatesV.m`, a few seconds): each
curve Y has W of index four in the Atkin–Lehner group, the star has genus zero and the three intermediate
quotients genus one; over E = X/{1,3,5,7,15,21,35,105} (Y → E by w_7, resp. w_3) the branch places Z_-420 and
Z_-84 have fields Q(√−5) and Q(√7). So Z_-420 | Z_-84 is the only Galois-stable pairing, each place is a fibre of
E over the star, a hyperelliptic Y would be E(√(cv)) with v a function on the star, and K(star)(√(cv)) would be a
genus-zero quadratic subfield of the Klein-four extension K(Y)/K(star), whose three quadratic subfields have
genus one. The hypothesis on the fields is needed: X_0(210,13)/{1,2,21,42,65,130,1365,2730} and
X_0(1326,1)/{1,17,78,1326}, hyperelliptic by shape D, meet every other hypothesis over one of their elliptic
quotients. Over the committed data the hypotheses hold for 144 genus-3 quotients: 142 ruled out by earlier
stages, and this pair.

Shape T, a second proof for X_0(6,35)/{1,3,5,15} (`tests/HyperellipticCertificatesT.m`, a few seconds). Over
E = X/⟨w_7⟩ (Jacobian 30a6) the branch places Z_-420 and Z_-84 are defined over different quadratic fields, so
they are the only Galois-stable pairing; both are fibres of E → E/w_2 (genus 0), so their difference is the
divisor of a function x and the hyperelliptic cover would be E(√(cx)). At 3 the vertical part of div x is zero,
so that cover is étale over the singular point of thickness two, where the dual graph of the curve has a single
singular point (`CoverIsNotSquareRootAtPrime`, `ComponentGroups.m`).

## The building blocks

| block | where |
|---|---|
| CM points: which d occur on X_0(D,N), their weights | `CarriesCMPoint`, `CMPointWeight` (`DualTest.m`); `tests/CMPointCounts.m` |
| dual test: a CM divisor certified torsion against W | `DualTestSpace`, `ReadWeightThreeHalvesBasis` (`DualTest.m`), with the sharp, 2 ∥ N and odd variants; `tests/_offline/DualTestReproduces.m` reproduces #146–#157 |
| weight-3/2 bases | the directory `DUALW_CACHE`, else `~/dualw` (`DualTestBasisDirectory`); checksums in `data/dual/bases.md5`, download script `tools/fetch_dual_bases.sh` (its Zenodo address is still a placeholder) |
| torsion bound | `JacobianTorsionBound` (`JacobianTorsion.m`) |
| component groups, Molina | `ComponentGroupData`, `CMSpecialisation`, `ClassInComponentGroup`, `TorsionOrderBoundOfCMDivisor` (`ComponentGroups.m`) |
| Borcherds certificates | `WriteBorcherdsForms`, `BorcherdsFormsFromFile` (`BorcherdsForms.m`), files in `data/borcherds/`, `data/dual/` |
| exact covers, heights | `DoubleCoversWithBranchDivisor`, `MatchCoverToCounts` (`ExactDoubleCover.m`), `ModelByHeegnerHeights` (`HeegnerHeights.m`, with `Certified` and `OldformPrime`) |
| certified L-values | `CertifiedLValueAtOne`, `CertifiedLDerivativeAtOne`; `tests/CertifiedLValues.m` |
| Noether's count | `IsHyperellipticFunctionField` (`HyperellipticDetection.m`) |

Measured basis costs on lovelace (single thread, 100–125 GB): level 4620 19 h, 6780 43 h, 7788 50 h,
8260 61 h, 7140 99 h, roughly dim^4 in the dimension of S_{3/2}; precision is cheap and must reach every
target index (B ≥ |d|/4).

## The stage

**Verify, do not search.** Each certified curve has a *certificate file*, `data/certificates/<D>_<N>.m`
(one file can hold the certificates of several quotients of the same curve), read the way banked Borcherds
forms are read: it names the shape (A, B, B2, C, C4, D, N, O) and holds what the proof needs — the CM
divisors or fibres, the banked forms, the levels, precisions and checksums of the dual tests, the
component groups, the star's label. The stage `FilterByHyperellipticCertificates`, last among the deciding
stages, loads the certificate of each undecided curve, rechecks it (`VerifyHyperellipticCertificate`), and
records the verdict with a label naming the shape and the certificate. The searches that
*found* the divisors (the Nelson screen, the scans in `tools/proof-77-10/`, the Borcherds searches) stay as
tools that write certificate files; the pipeline never depends on them.

**The bases.** A dual test needs the weight-3/2 basis of its level. The stage reads it from the cache
directory (`DUALW_CACHE`, default `~/dualw`), refuses a file whose checksum differs from the one in the
certificate, and leaves the curve undecided, with a message naming the missing basis, when it is absent.
Recomputing is opt-in (`Compute`). The bases go to Zenodo (decision 4); the deposit is prepared on lovelace
(`~/zenodo_dualw`, 749 MB, checksums verified) and waits to be uploaded from Eran's or Sachi's account.

**Cost.** With the bases present the stage runs in parallel (`run_par`) in about 20 minutes on lovelace
with 32 jobs (1,277 s on 2026-11-08); `tests/HyperellipticCertificates.m` takes about 6 minutes. Without
the bases, the large levels cost a few hundred CPU-hours (below).

**Registering a stage** touches: `FILTER_STAGES` (`workingcode.m`), `run_pipeline.sh` (a `run_seq` line and
an `UpdateCurvesAfter…` step), `run_sequential_stage.m`, `tests/PipelineStages.m`, `analysis_stages.m`,
`reconstruct_attribution.m`, `make_latex_tables.m`, `make_latex_figures.m` (a token group) and
`docs/RUNNING_PIPELINE.md`; data files are written with `SetColumns(0)`.

## Order of work

1. **Land 61 → 34 in the data.** Done 2026-11-05 up to six quotients: the rerun leaves 40 undecided, the
   34 and the six of `open_gaps` (see "The data" above); `tests/OpenThirtyFour.m` checks both lists.
   390_1 done 2026-11-06; the three special fibres certified 2026-11-08; the 6_35 pair by shape V
   2026-11-09 (and 6_35/{1,3,5,15} also by shape T). Nothing remains.
2. **The dual test as an intrinsic.** Done 2026-11-07: `DualTest.m`; `tests/_offline/DualTestReproduces.m`
   reproduces all nine cases from the saved bases (`data/dual/bases.md5`). One CM dictionary and one dual test, with the 2 ∥ N and odd-D
   variants; the basis cache with checksums. Acceptance: it reproduces, from the saved bases, every
   certified relation and every refused control of #146–#157.
3. **The certificate format and the stage, one shape at a time:** A (46_31 first: no basis needed), then
   B, then C, then D. Started 2026-11-06: `HyperellipticCertificates.m`, `data/certificates/46_31.m`,
   `tests/HyperellipticCertificates.m`; the stage `FilterByHyperellipticCertificates` runs last (after
   `FilterByNonALInvolutions`, before `UpdateCurves8`), so it receives exactly what no filter decides. Shape A
   certifies all seven of its curves: 46_31 (Borcherds forms), 4130_1 and 6_379 (sharp dual test), 21_58,
   77_10, 39_34 (2 || N dual tests, with the check that the Jacobian is old at 2), 6_649 (sharp dual test,
   torsion bound 6 lowered to 2 by Φ₃). Shape B (branch divisor recomputed from the fixed points, a
   torsion-free genus-two star) certifies 2030_1, 2418_1, 6_403 (Borcherds forms) and 3390_1 (dual test).
   2002_1 is shape B through the sharp test at level 4004, precision 2200. Shape C (genus three over an
   elliptic star; Eran's choice 2026-11-08: a general verifier) checks the halves of the branch divisor,
   identifies E exactly (traces to the Sturm bound of level DN, Tamagawa numbers against the component
   groups, unique in its isogeny class), makes the pairing torsion by the dual test (on E's eigenspace
   where given), zero by the component groups and the reduction of E(Q)_tors, and the cover untwisted by
   a simple Prym or by the dual graph of Y (CoverIsUntwistedAtPrime): 15_38, 6_385, 714_5. The pipeline
   ends with 24 undecided at that point. Then 6_341 (shape B2), 6_323 (shape C4, decision 2), the sixteen
   of group D (shape C for ten, shape D for six, decision 8), the 34th (shape N), and the three special
   fibres (shapes C and O): done 2026-11-08. Each shape is accepted when the stage decides its curves from
   the certificates alone
   and the per-curve tests still pass.
4. **#128 made exact** (decisions 1 and 3). Done 2026-11-08 as shape D for the six curves of group D that
   shape C does not decide; the 34th went through shape N instead (decision 8).
5. **Labels of stored models.** Done (decision 6): only models on `data/models/proved_models.m` decide a
   curve; 390_1 is on the list (2026-11-06).
6. **The full run** on lovelace from this branch, from stage 1. Done 2026-11-08 at 05f4886 (an export of
   the commit, every stage into an empty directory; 32 jobs, 1 h 26 min): it ends with the 6_35 pair
   undecided and no other. `stage_counts.txt` and every summary and table are identical to the committed
   ones, and so is every curve file except for the witness of `FilterByAutomorphismGroup` on two of its 46
   curves, X_0(72) and X_0(120)/{1,15}: both are ruled out in both runs, by a different failing subgroup of
   the same order (the stage reports the first of the smallest, and their order is not fixed between
   runs). The 8 fast and 6 offline tests of the certificate stage pass on the fresh data. With that order
   fixed (#160, merged here as 322b2cb), a second run from stage 1 at 322b2cb (2026-11-09, 1 h 26 min) differs
   from the committed data only in the witnesses of `FilterByAutomorphismGroup` on the nine curves #160 lists,
   with the same verdicts; the committed data now comes from that run, with stages 41–42 rerun for shape V.

## Decisions (Eran, 2026-11-05; 7 and 8 later)

1. **The 34th is proved by making #128's enumeration exact for that curve.** L-values enter as exact
   ratios from modular symbols instead of zero-tests at 1e-8, heights and coordinates are recognised
   exactly instead of by `BestApproximation`, and the certificate shows that no branch place needs the
   window |n| ≤ 12. Then "one cover consistent with the trace formula, a plane quartic" is a proof, rerun
   by a test. Superseded by decision 8: the 34th is shape N.
2. **6_323: the Gross–Zagier step is checked exactly.** The certificate computes the root numbers of
   E^(−8) and E^(969) and the rational number L(E^(969), 1)/Ω by modular symbols, and the stage checks it
   is 0; with the signs this gives order of vanishing ≥ 3. The theorem of Yuan–Zhang–Zhang is cited as an
   input, as Molina's theorems are for the component groups.
3. **The sixteen curves of group D: #128 made exact, the same way as for the 34th.** Point counts make the
   Jacobians isogenous, not the curves isomorphic, so the positive verdicts also rest on the enumeration
   missing nothing. Each certificate records its run's inputs (the curves of the isogeny class tried, the
   Mordell–Weil proofs, Manin constant 1) and the oldform multiplier is derived, or shown unneeded, at every
   place where it is used (15_46 at d = −15, p = 2 first). Superseded by decision 8 for ten of them;
   the other six are shape D, which is this decision carried out.
4. **The bases go to Zenodo**, with checksums in the repository and a script that downloads them into the
   cache; the deposit is created from Eran's or Sachi's account.
5. **The list of 34 is a data file**, `data/open_34.m`, with a test that the curves left undecided before the
   new stage are exactly these; certificates are keyed by the same labels.
6. **A stored model decides a curve only when the model is proved**: `FilterByCanonicalQuadrics` uses a model
   at D > 1 only for the curves on `data/models/proved_models.m`, each with the argument that proves its
   equation.
7. **The three special-fibre quotients get a stage that also checks Z mod p.** Superseded 2026-11-08: Z mod 3
   is hyperelliptic for 14_51 and 38_15 (and must be for 158_3), so no check of Z decides them; all three are
   hyperelliptic and the certificate stage proves it (shapes C and O).
8. **Group D and the 34th go through the certificate stage, not through #128 (2026-11-08; supersedes 1 and
   3 for the classification).** The sixteen are genus-three double covers of an elliptic star, so shape C
   proves them from the branch divisor alone (dual test, component groups, twist), with no enumeration;
   the 34th is shape N: the only Galois-stable pairing of its branch places, {P_-4, P_-12} | Z_-1128, has
   class 2 in Phi_3 = Z/7, so it is not hyperelliptic. #128's models stay as equations; making #128 exact
   is needed only if those equations are to be theorems. Done 2026-11-08 by shape C: 39_11, 6_185, 141_2,
   21_34, 55_7 (a_3(E)-eigenspace), 57_10, 15_46 (the odd variant), 6_253, 6_247, 6_301. Shape C did not
   reach the other six — 6_455 and 690_7 (bases at levels 5460 and 9660 still computing), 22_51 (no twist
   test applies), 6_187 (a mixed pairing over Q(√−2)), 1326_1 and 210_13 (W of index four) — and they are
   shape D: #128's route made exact, with certified L-values (decisions 1 and 3), done 2026-11-08.
   Two counting fixes made on the way: the Heegner weight of P_-4 doubles when w_2 is in the group
   (CMPointWeight), and the number of points of X_0^D(N) over a star point is #W over the stabiliser in W,
   without the units of the CM order (cg_stab; it had been wrong at P_-3 and P_-4 off the fixed points).

## Branch

`classification-pipeline` merges the open PRs as they stand on 2026-11-05. Conflicts were resolved by keeping
both sides in the provenance log, the paper README, `.gitignore` and the spec file; in `SchoferFormula.m`
the two fixes to the field of definition (no Artin map when the field is K; fall back to the slow routine
when the Artin map refuses the twist ideal) are combined, as #66 already does. When a PR here changes, merge
its branch again rather than editing it on this branch.

Checked on the merged branch when it was made (2026-11-05; not rechecked since): 13 of 14 fast tests pass (CM fields of definition, the divisor
point at 21_2, dual graph, Kodaira, component groups, plane quartics, 15_38, 46_31, 2030_1, 6_341, 4130_1).
`tests/HyperellipticCRV.m` (#108) failed on a curve that #112 adds: X_0(38,7)/{1,7}, genus 7, stored as a
pair over a conic that has a rational point, so it has a model y^2 = f(x) and should be stored as one, as
#108 does at 10_53. Fixed on #112's branch 2026-11-09 (31b661f, the entry converted, trace formula at eight
primes) and merged here; the test passes on this branch.
