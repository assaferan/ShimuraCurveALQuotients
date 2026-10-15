# Hyperelliptic detection by Noether's count (2026-10-14)

A curve of genus g >= 3 is hyperelliptic iff the pairwise products of a basis of its holomorphic
differentials satisfy (g-1)(g-2)/2 linear relations (canonical image a rational normal curve);
otherwise there are exactly (g-2)(g-3)/2 (Max Noether). Two ways to get the differentials:

* a function field (`IsHyperellipticFunctionField`, HyperellipticDetection.m, PR #108): every curve
  stored in the CRV shape, and any cover whose double covers are stored (fibre product);
* for D = 1, the weight-2 cusp forms on Gamma_0(N) fixed by W, as q-expansions past the Sturm bound
  of weight 4 (`decide_d1.m`); controls `d1_controls.m`: Ogg's hyperelliptic X_0(40), (48), (47),
  (59) and the non-hyperelliptic X_0(34), (45), (64), (53), (55), all right.

Results on the 299 stored CRV curves (`hypsweep.m`, `hypclass.log`): 60 hyperelliptic, all by an
Atkin-Lehner involution; 58 over a conic without a rational point; the two covers of X_0(10,53) by
w_530 and w_53 over the rational conic y^2 = 5x^2 - 8x + 16 got models y^2 = f(x) (`hyp1053c.m`).

The 173 undecided curves of the 2026-09-30 rerun (`census_hyp.m`, `census_und*.m`): 7 at bases
with model files (all seven NOT hyperelliptic, `decide7.m`: 390_1 x2, 462_1, 510_1 x2, 22_15,
10_51); 7 at level 420 with the composite-level branch's files (all NOT hyperelliptic,
`decide420.log`: 6_35 x2, 10_21 x2, 14_15 x3); 22 at D = 1 (`decide_d1.log`: 19 NOT hyperelliptic,
3 hyperelliptic -- X_0(114)/{1,2,19,38}, X_0(130)/{1,2,13,26}, X_0(195)/{1,5,39,195}, the three
curves PR #51 found hyperelliptic with explicit models). The remaining 137: 132 at 114 bases with
D > 1 whose star curve has genus > 0 (no genus-0 quotient, so no Borcherds model and no fibre
product can exist), and 5 at genus-0 bases not yet run (15_14; 210_13 at M = 5460; 798_1, 870_1,
1326_1 at M = 1596, 1740, 2652).

## 2026-10-14, evening: the dual graph at non-squarefree N (PR #114) and the 61 open cases

`decided83.txt` lists the 83 undecided curves the extended dual-graph filter rules out (built with
the Atkin-Lehner ideal from the reflection of the lattice chain, `alideal3.m`); `open61.txt` is what
remains after #113 and #114 (`open61.m`): 55 of genus 3 and 6 of genus 4; 12 at bases whose star curve
has genus 0 (the model route: 7 already known not hyperelliptic from the composite-level files, 15_14
and the three next-wall bases to run), 38 of index 2 over a star of genus 1 or 2, 11 star curves of
genus 3; only two with N not squarefree (22_63, 34_25, both star curves).

## voight-willis/ — MOVED

The Voight-Willis expansion code is a general model-production method, comparable with the
Guo-Yang/Borcherds route, not a hyperellipticity probe; it now lives at
`vvdata/weyl-campaign/voight-willis/`.  It reached this directory only because the motive for
building it was Noether's count on the quotients still undecided here.

## star-elliptic/ — the fixed-point criterion for star curves (2026-10-15)

For a star curve W is the whole Atkin-Lehner group, so a hyperelliptic involution would be an
EXCEPTIONAL automorphism; Rotger (Compositio 132 (2002) 229-241) proves Aut = W at level 1 only
under hypotheses that fail here, and his V_291/W is an explicit star curve whose hyperelliptic
involution IS exceptional, so no general theorem is available.  What survives his Proposition 1 is
that Aut(X*) consists of RATIONAL automorphisms.  The criterion: label each elliptic point by
(stabiliser order, field of definition), both preserved by Aut; an involution fixes at least one
point of every odd-size label class, and being rational it then fixes that point's whole Galois
orbit; so the sum of the field degrees over odd-size classes is a lower bound for its fixed points,
against 2g+2 for a hyperelliptic involution.

`ell.m` lists the elliptic points, `ell2.m` applies the criterion, `ctrl.m` runs it over the
known-hyperelliptic star curves (where it must never fire) and `nonhyp.m` over the known
non-hyperelliptic ones (to see whether the statistic separates them at all).
⚠ The control is WEAK as it stands: over 54 known-hyperelliptic star curves the criterion never
fires, but nearly every sum is 0 and the largest is 3, so it would pass even if the argument were
wrong; several runs are also incomplete because the field routine supplies no field at some
discriminants, mostly at non-squarefree levels.  Treat the two verdicts it produces
(X_0^6(379)* and X_0^4130(1)* not hyperelliptic) as unconfirmed until the separation is measured.
Corroboration so far: it also fires on X_0^26(85)* and X_0^86(15)*, which the metric dual graph of
PR #116 independently rules out.

## open-cases/ — triage of the quotients still undecided (2026-10-15)

`zcomp.m` the special-fibre component at p || N (rules out 3), `cq635.m` the canonical quadrics on
the composite-level models of X_0^6(35), `ecurves.m` identifies the elliptic Atkin-Lehner quotients
in Cremona's tables with rank and torsion, `fixgal.m` the Galois structure of the fixed points
(yield 0), `mod2.m` the mod-2 Galois image test (yield 0; ⚠ p = 2 must be EXCLUDED, since J[2] is
not etale there and the reduction of the Weil polynomial is not the characteristic polynomial of
Frobenius on J[2] — the only apparent hits were at p = 2), `s8.m` the allowed characteristic
polynomials, `stars.m` which open quotients are star curves, `bielliptic.m` the reach of the
bielliptic fibre-product shortcut, `have.m` model coverage.
⚠ `have.m` must be run against a MERGED model directory assembled from every branch: the model
files are scattered, and a single worktree reports no model for any open base.
