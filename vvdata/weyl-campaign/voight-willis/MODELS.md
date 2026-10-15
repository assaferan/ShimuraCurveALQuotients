# Equations built from the Voight-Willis expansions, 2026-10-15/16

Every curve below was produced by `model.m` followed by `finish.m` and nothing else: a basis of
Atkin-Lehner invariant forms from the power series expansions, the branch points of the double
cover over the star curve, the invariants of that branch set recognised as rationals, Mestre's
construction over Q, and the quadratic twist solved for against the trace formula.  Nothing was
taken from the repository's own models or from a published equation except where the line says so.

## What each was checked against

* **Point counts.**  `#C(F_p)` and `#C(F_{p^2})` against `ComputePointsViaTrace` at every good odd
  prime below 150 -- 64 comparisons per curve.  This is also what fixes the quadratic twist, which
  a complex-analytic construction cannot see: the twisting number is read off from the Legendre
  symbols the counts give and solved for by linear algebra over F_2, and the solution is then
  confirmed by the counts it was not fitted to.
* **Genus**, from `GenusShimuraCurveQuotient`, against the genus of the equation.
* **The number of branch points**, 2g+2, which the fit has to produce on its own.
* **Conductor**, for the genus-1 cases: Jacquet-Langlands puts the Jacobian in the D-new part of
  J_0(DN), so the conductor must divide D*N and be divisible by every prime of D.  All four at
  X_0^10(9) pass, and this is independent of the point counts.

## Validation, against two outside sources

**X_0^35(1)/w_5**, genus 2 (not one of the open bases; this is the control).

    y^2 + (x^2 + 1) y = -7x^5 + x^4 + x^3 + x^2 + 5x - 4

* isomorphic over Q to the model this repository already stores for that quotient;
* isomorphic over Q to Molina's published equation (arXiv:1004.3675, Table 2);
* the six branch points computed from the expansions are carried onto the six roots of each of
  those by ONE Mobius change of the coordinate on the star, total residual 1.2e-31 over six
  points, of which three are free after the three that fix the map;
* 64 point counts.

The absolute Igusa-Clebsch invariants came out as 1498959544523423744/102942875,
-139237160005402624/102942875 and -10517118550948864/102942875, exactly the published values, to
a relative 2e-29.  The twisting number between Mestre's representative and the curve over Q was
3000216155931018232719609781369.

## New equations: X_0^10(9), all four index-2 Atkin-Lehner subgroups of genus 1

The repository holds no model for any quotient of this base.  Arithmetic volume 8, the smallest of
the sixteen open bases.

| subgroup W | genus | equation | conductor |
|---|---|---|---|
| {1, 9, 10, 90} | 1 | y^2 + xy + y = x^3 - 19x + 26 | 30 (30a2) |
| {1, 2, 9, 18} | 1 | y^2 + xy + y = x^3 - x^2 + 52x - 53 | 90 (90b3) |
| {1, 5, 9, 45} | 1 | y^2 + xy = x^3 - x^2 - 69x - 235 | 90 (90a3) |
| {1, 2, 5, 10} | 1 | y^2 + xy + y = x^3 - x^2 + 13x - 61 | 90 (90c1) |

Each passed all 64 point counts.  The twisting numbers were 21347557, 1457, -12259 and 391281, in
every case the unique solution of the F_2 system (kernel dimension 0).

X_0^10(9) has three further index-2 subgroups: {1, 5, 18, 90} and {1, 2, 45, 90} have genus 0, and
the remaining quotients have index 4 and are fibre products of three of the above, which is how
the rest of the repository builds such covers.

## Not finished

**X_0^21(5)/{1, 5, 7, 35}**, genus 2.  The geometry is done -- six branch points, fit residual
3.5e-28, the gaps in the two decompositions 5.8e23 and 1.5e25 -- but at a truncation target of
1e-28 the Igusa-Clebsch invariants are not determined to enough digits for the recognition to be
certain, and the recognition refuses rather than guessing.  Genus 2 wants a target near 1e-38
because the invariants have seventeen-digit numerators over nine-digit denominators; genus 1 is
satisfied with far less, which is why the four above came out at 1e-38 with room to spare.  The
branch points are in `out/branch_21_5_1-5-7-35-.m` and a rerun at a deeper target is all that is
needed.
