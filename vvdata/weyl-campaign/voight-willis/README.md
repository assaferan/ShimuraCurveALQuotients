# Power series expansions of modular forms on Shimura curves (Voight–Willis)

A second, independent way to produce **equations for Shimura curves and their Atkin–Lehner
quotients**, alongside the Guo–Yang / Borcherds route the rest of this repository uses.

Reference: J. Voight and J. Willis, *Computing power series expansions of modular forms*, in
Computations with Modular Forms (Springer, 2014); author's copy at
`jvoight.github.io/articles/ModFormPaper-090317.pdf`. The implementation follows its Section 3.

Where the Borcherds route builds a model from an eta-quotient basis obtained by enumerating lattice
points of a polytope, this one expands modular forms numerically in a disc coordinate on the
Fuchsian group's fundamental domain, reads off a basis of forms of each weight, and builds the model
from their ratios. The two are limited by completely different things, which is the point of having
both:

| | Borcherds / Guo–Yang | Voight–Willis |
|---|---|---|
| what blocks it | the number of divisors of M = 2DN (the Normaliz wall at 24), and the m = 0 multiplier when N > 1 | the size of the fundamental domain, which grows with the **arithmetic volume**, not with D·N |
| exactness | exact equations | numerical, so a verdict needs an error statement or an exact post-check |
| conditionality at N > 1 | rests on the m = 0 term | none |

So a base can be out of reach of one and inside the other, and measurement (below) shows the first
bases behind the divisor wall are all well inside this method's reach.

## Status

**Validated at X_0^6(1) and X_0^10(1), base curve and all Atkin–Lehner quotients.** For every
Atkin–Lehner subgroup of each base and weights 2, 4 and 6, the computed dimension matches the one
predicted from the quotient's own signature: **15 of 15 at D = 6 and 15 of 15 at D = 10**, with
gaps of 10^4 to 10^12. At weight 2 that prediction is the genus of the quotient, so this is also an
independent analytic confirmation of `GenusShimuraCurveQuotient` and `NumberOfEllipticPointsByCMOrder`.
D = 14 and D = 15 do **not** yet work; the reason is identified below and is not the one that
blocked D = 10.

**Validated at X_0^6(1)** (signature (0; 2,2,3,3), domain radius 0.5435):

* `dim M_k` for k = 2, 4, …, 14 reproduces the orbifold dimension formula (0, 1, 1, 1, 1, 3, 1),
  with the vanishing singular values separated from the rest by 1e6 to 1e9;
* the Hauptmodul reproduces Elkies (*Shimura curve computations*, arXiv:math/0005160, Sec. 3.1 and
  Table 1, transcribed with its corrections in `tests/CMFieldsElkies.m` on `main`): the cross-ratio
  of the elliptic points is −1 to 7.7e-15, and his t = 0 at discriminant −24 comes out as 1e-34;
* a check **not** forced by symmetry: the valence formula gives deg div(f_4) = 2/3, carriable only
  by the two order-3 points at 1/3 each, and the computed weight-4 form vanishes at exactly those
  two, at relative magnitude 8.4e-21 against 0.503 and 1.000 elsewhere.

## What has to happen next, in order

1. ~~The Atkin–Lehner layer.~~ **Done and validated at D = 6 and D = 10**, see below.
   D = 14 and D = 15 need the side-pairing sampling described below first.
2. **Agreement with the models we already have.** Still open, and it is the step that matters; see
   "What is still NOT done" below for why it is more than a thin wrapper and where to start.
3. **Beyond the Borcherds reach.** The timings below say this is now worth attempting: the first
   bases behind the divisor wall compute their domains in seconds to minutes.

## Two things that had to be fixed, and how they were found

Both were found by asking first whether the **base** library (Gamma relations only, no
Atkin–Lehner) resolves at the base in question — `basediag.m` does exactly that. At D = 10 it did
not, which ruled out the Atkin–Lehner layer as the cause before any time was spent on it.

**1. The disc was centred in the wrong place.** Magma's `UnitDisc` defaults to
`Center := (9/10)*i`, a point of the upper half plane fixed in advance and unrelated to the group.
The series is expanded about that point, so the circumradius rho of the fundamental domain — which
sets both the number of terms (rho^N) and the conditioning of the linear system — is an accident of
that default. `RecentredDisc` in `vwal.m` iterates: compute the domain, take the Euclidean centroid
of its vertices, map that back to the upper half plane and rebuild the disc there, keeping the
better of the two.

| base | rho, Magma's default centre | rho, recentred | vertex moduli before | after |
|---|---|---|---|---|
| X_0^6(1) | 0.5435 | **0.4472** | 0.41–0.54 | 0.41–0.45 |
| X_0^10(1) | 0.8212 | **0.6482** | 0.05–0.82 | 0.49–0.65 |

At D = 10 the default centre gave a domain spanning |w| = 0.05 to 0.82. The symptom was a smooth
geometric ramp of singular values with no gap anywhere — and, decisively, **the ramp moved down with
the truncation floor as terms were added but never separated** (234 terms at 160 digits put the
floor at 1e-20 and the ramp simply followed it). That is what rules out "not enough terms" and
points at conditioning instead. Recentred, the same computation gives gaps of 10^4 to 10^12.

**2. The zero cutoff was a thousand times too loose.** `ZeroLevel` returned `1e6 * rho^N`. Measured
at D = 6 and D = 10 with the disc recentred, a genuine kernel vector sits within about 100 times the
truncation floor rho^N while the smallest non-kernel singular value sits 1e4 or more above it, so
`1000 * rho^N` separates them cleanly. The old value swallowed the first non-kernel value and
reported one dimension too many at weight 2.

## What limits the reach: how OVER-DETERMINED the system is

Measured 2026-10-15, and it supersedes the earlier reading that D = 14 fails for the same reason as
D = 15.  **D = 14 resolves.**  What it needed was more sample points, nothing else.

`ReducedSamples` places Q points on the sampling circle and keeps only those that reduce to a
non-identity element.  Retention varies enormously between bases -- about 15 percent at D = 15, 27
percent at D = 26, 100 percent at D = 22 -- so setting `Q = 3*NT`, as the drivers did, controls the
number of points ATTEMPTED and leaves the number of rows actually obtained to chance.  When that
lands near or below NT the system is square or underdetermined and the singular values come out as
a smooth geometric ramp with no gap anywhere, which is easy to misread as the method failing.

| base | genus | rho | NT | rows at Q = 3 NT | result | rows raised | result |
|---|---|---|---|---|---|---|---|
| X_0^6(1) | 0 | 0.447 | 58 | 73 | resolves | | |
| X_0^10(1) | 0 | 0.648 | 107 | 189 | resolves | | |
| X_0^14(1) | 1 | 0.693 | 126 | 148 (1.2x) | ramp | **984 (7.8x)** | **resolves, gap 67** |
| X_0^15(1) | 1 | 0.736 | 151 | 66 (0.4x) | ramp | 878 (5.8x) | still a ramp |
| X_0^22(1) | 0 | 0.839 | 263 | 494 (1.9x) | ramp |  4117 (15.6x), not finished | |
| X_0^26(1) | 2 | 0.774 | 180 | 147 (0.8x) | ramp | 977 (5.4x) | **still a ramp** |

So the first rule is to size Q by the rows RETAINED, aiming for several times NT -- D = 14 needed
about eight times -- and not by NT itself.  Genus is NOT the obstruction: D = 22 has genus 0 and
fails where D = 14 has genus 1 and succeeds.

But over-determination is necessary and **not sufficient**: D = 26 still gives a ramp at 5.4 times.
What the two bases that resolve easily have in common is a small domain radius, and the table is
monotone in it -- 0.447 and 0.648 resolve at under 2x, 0.693 resolves only at 7.8x, and 0.736,
0.774 and 0.839 have not resolved at any row count tried.  So the radius governs the conditioning
and the row count governs whether that conditioning can be exploited; both have to be right.
Shrinking the radius by sampling on a smaller circle does not by itself help (D = 22 at 70 percent
of the maximal vertex modulus is still a ramp), because the terms needed fall with it and the
relations get weaker in step.

**D = 15 is the genuine exception** and its diagnosis stands: raising the rows to 878, nearly six
times NT, still gives a ramp.  Its reductions barely move the point -- sources at 0.7358 land in
the narrow band 0.666 to 0.736 -- so the relation f(delta w) = j^k f(w) ties two nearly equal points
and carries almost no information, however many of them there are.  Compare D = 10, where sources at
0.648 spread over 0.188 to 0.645.  The fix there is still to pair points across the boundary with
the side-pairing generators from `ShimFDGens` rather than relying on reduction.

## The wall, measured — and it is the volume, not D·N

`FundamentalDomain` is the prerequisite and the bottleneck, but **D·N is a bad proxy for its cost**;
`ArithmeticVolume` is a good one. For a fixed D·N the split between D and N changes the volume by a
factor of several, because each ramified prime p contributes p − 1.

| base | D·N | volume | domain | sides |
|---|---|---|---|---|
| X_0^6(1) | 6 | 0.33 | 0.19 s | 6 |
| X_0^10(1) | 10 | 0.67 | 0.21 s | 6 |
| X_0^6(35) | 210 | 16 | 248 s | 102 |
| **X_0^330(1)** | 330 | 13.3 | **9.3 s** | 70 |
| **X_0^510(1)** | 510 | 21.3 | **120 s** | 116 |
| **X_0^462(1)** | 462 | 20 | **233 s** | 98 |
| **X_0^6(55)** | 330 | 24 | **317 s** | 150 |
| X_0^77(10) | 770 | **180** | not finished in 30 min | — |

The earlier reading, that the method stops somewhere below D·N = 770, was an artefact of the base
chosen to test it: X_0^77(10) has volume 180 because D = 77 contributes (7−1)(11−1) = 60, eleven
times the volume of X_0^6(35), which does finish. **The bases behind the Borcherds divisor wall are
not like that.** X_0^330(1), X_0^462(1) and X_0^510(1) have volumes of 13 to 21 and their domains
come out in seconds to minutes, so the reach of this method and of the polytope route genuinely do
differ, in the direction that matters.

Volumes of the other bases at those levels, for planning: X_0^10(33) 32, X_0^22(15) 40,
X_0^6(77) 32, X_0^22(21) 53, X_0^10(51) 48, X_0^6(85) 36.

## The Atkin-Lehner layer (`vwal.m`)

Implemented and validated at X_0^6(1): for all five Atkin-Lehner subgroups and weights 2, 4 and 6,
**15 of 15** computed dimensions match the dimension predicted from the quotient's own signature,
with gaps of 10^9 to 10^10. At weight 2 that prediction is the genus of the quotient, so this is
also an independent analytic confirmation of `GenusShimuraCurveQuotient` and of
`NumberOfEllipticPointsByCMOrder`, which supply the signature.

Three things the work turned up, all recorded in the code:

* The Atkin-Lehner element of reduced norm m is found in the **two-sided ideal** of that norm, where
  every element of reduced norm exactly m is a generator and so normalises the order. Searching the
  order itself with a box of 4 misses several; searching the ideal with a box of 3 succeeds at every
  base tried. Magma's `Normalizer` is implemented for definite orders only and cannot be used.
* The reduced norm **must come out positive**. A generator of negative norm sends the upper half
  plane to the lower one and does not act on the curve, and `IsPrincipal` returns exactly such a
  generator at D = 6, m = 3.
* `ShimuraReduceUnit` **fails on 20 to 30 percent** of Atkin-Lehner images. Applying mu throws the
  point to modulus up to 0.997, and the reduction then returns a point still outside the domain,
  with modulus up to 0.91 against a domain radius of 0.54, where a 50-term series is worthless.
  Those rows are dropped. On the rest the relation holds to 3e-10, against 3e-11 for the group's own
  relation, which is how we know the relation is right and only the reduction is at fault.
  Recovering those samples would need a better reduction and is not done.

**The number of terms and the vanishing tolerance must follow the domain radius.** At D = 6 the
radius is 0.54 and 50 terms suffice; at D = 10 it is 0.82, where 50 terms leave an error of 1e-5 and
a fixed tolerance of 1e-7 reports every space as empty. `TermsFor` and `ZeroLevel` in `vwal.m` set
both from the radius. The cost of this is real: D = 10 needs 164 terms where D = 6 needs 53.

## Agreement with published equations: DONE, and it holds

The step that makes this a second source rather than a second guess.  `starcmp.m` builds the
Hauptmodul of a **star** curve from the expansions -- which exercises the Atkin-Lehner layer, since
the star is the quotient by the full group -- and compares it with Elkies' published coordinates
(*Shimura curve computations*, arXiv:math/0005160, Sec. 3.1 and Table 1 for D = 6, Sec. 4 and
Table 3 for D = 10, transcribed with their corrections in `tests/CMFieldsElkies.m` on `main`).

A computed Hauptmodul is defined only up to a Mobius change, so the test is whether **one** Mobius
map carries the published value onto the computed one at **every** CM point at once.  Writing
s = (a t + b)/(c t + d) gives one linear equation per point in four unknowns defined up to scale, so
three points fix the coordinate and every further point is a prediction.  The smallest singular
value of that system is the residual of the fit and the next one its conditioning; nothing is
fitted point by point.

| base | weight | rho | terms | CM points | predictions | residual | next sv |
|---|---|---|---|---|---|---|---|
| X_0^6(1)* | 24 | 0.4472 | 58 | **14 of 14** | 11 | 5.4e-13 | 1.30 |
| X_0^10(1)* | 12 | 0.6482 | 107 | **12 of 12** | 9 | 4.1e-16 | 1.51 |

Per-point relative residuals are 1e-13 at D = 6 and 1e-15 to 1e-16 at D = 10.  **D = 10 is the case
with real content**: its star curve has four elliptic points (orders 2, 2, 2, 3) at t = 0, oo, 2, 27
for discriminants -3, -8, -20, -40, so their cross-ratio is an invariant no change of coordinate can
fake, and the twelve points fitted include all four.  At D = 6 the star has only three elliptic
points, which any Mobius map matches, but the eleven further CM values carry the content there.

Convergence is clean, which is what rules out a coincidence (D = 6, weight 24):

| terms | 46 | 69 | 92 |
|---|---|---|---|
| kernel | 7.5e-10 | 1.0e-15 | 3.9e-22 |
| fit residual | 1.8e-3 | 1.1e-6 | 1.6e-9 |

## Why the comparison against OUR OWN stored models is vacuous at D = 6 and D = 10

Worth stating so nobody builds it twice.  For these bases every quotient has genus 0, and a stored
model of an index-2 quotient is a double cover y^2 = f_W(x) of the star whose branch points are the
roots of f_W.  Those branch points are exactly the elliptic points of the star not fixed by W, so
the whole complex-geometric content of the stored file is the **position of three points** on a
genus-0 curve -- and any three points are equivalent under a Mobius change.  At D = 10 the models
give the three order-2 points (x = 0, 1, 2/27 for the points fixed by w_2, w_5, w_10 respectively)
but say nothing about where the order-3 point sits, so there is no fourth point to over-determine
anything.  The comparison therefore has **zero independent predictions** at both bases, which is why
the Elkies CM tables, not the model files, are what the expansions are tested against above.

A non-vacuous model comparison needs a base whose stored quotients have **positive genus**, where
the j-invariant or the Igusa invariants are genuine moduli rather than a coordinate choice.
Separately, a complex-analytic method cannot see the **quadratic twist** of y^2 = f_W(x) at all:
that is arithmetic, not geometry, and nothing here tests it.

## Traps, each of which cost a day

* **Magma's `SingularValueDecomposition(M)`** on an m × n matrix returns (m × n, m × m, n × n), and
  the singular values are the diagonal of the **first** return, descending; the second is the left
  orthogonal matrix. Reading its diagonal instead gives numbers bounded by 1 and of arbitrary sign,
  and produced an entirely spurious earlier conclusion that the method does not work.
* **Equilibrate the columns.** They scale like rho^n, so without it every singular value sits at
  rho^N and no gap is visible at all.
* **Magma's own `CMPoints` is broken at D = 6**, failing with a `RingClassGroup` coercion error
  inside `package/Geometry/GrpPSL2/GrpPSL2Shim/reciprocity.m` line 29. Take instead the fixed point
  of a trace-zero element of reduced norm m.
* The undocumented `ShimFDGens` and `ShimFDDisc` hold the side-pairing generators and the domain
  vertices; the Cayley conjugate of the half-plane matrix lands exactly in SU(1,1), so the
  automorphy factor is unambiguous at even weight.
* **The absolute zero cutoff is the wrong instrument once Atkin-Lehner rows are present.**
  `ZeroLevel` assumes the kernel sits near the truncation floor rho^N, which holds for Gamma rows
  alone.  Atkin-Lehner rows are limited by the quality of the reduction instead, so the kernel
  floats far above rho^N: at D = 6 weight 24 it sat at 7.5e-10 and 1.4e-8 against a floor of 1e-16,
  while the first non-kernel value was 0.0195 -- a ratio of 1.4e6, unmistakable, yet every kernel
  value was above the cutoff and the dimension came out 0.  Use `GapDimRel`, which takes the largest
  ratio jump among values below a RELATIVE threshold and uses no absolute cutoff at all.
* **Use `Embed(R, O)` to place CM points, never a box search.**  It returns an OPTIMAL embedding, so
  the fixed point really has discriminant disc(R) and not that of some larger order containing it,
  which a search for an element of given trace and norm does not guarantee.  It also succeeds where
  a box search simply fails: at D = 10 the generator of discriminant -40 lies outside a box of 6 on
  the order basis, and brute force over that box returns nothing at all, silently losing 9 of the 12
  CM points.  Note the reduced norm form is INDEFINITE here, so basis elements have negative norms
  and "small norm" is no guide to small coefficients.
* **`AttachSpec` with a relative path breaks when a script is run from its own directory**, and the
  failure is quiet: the intrinsics are simply undefined, so every check that uses them is skipped
  while the rest of the script runs on and prints plausible-looking output.  All scripts here now
  use the absolute path.

## Files

`vw3.m` is the working library and `vwal.m` the Atkin-Lehner layer (which also holds
`RecentredDisc` and `ZeroLevel`); `aldims.m` is the validation sweep over all Atkin-Lehner
subgroups of a base, `basediag.m` tests the base library alone at one base and is the first thing to
run when something does not resolve, `recentre.m` reports how much recentring gains, `aldiag.m`
instruments the reduction quality; `dims.m` computes dimensions of M_k, `haupt.m` builds a Hauptmodul
as a ratio of weight-12 forms, `cmchk.m` evaluates it at a CM point, `predict.m` is the
valence-formula prediction. `vw2.m`, `run2.m`, `shape*.m` and `dbg_*.m` are development history,
kept only because they record what was tried. `vwlib.m` and `NOTES-firstattempt.md` are from the
first attempt and are superseded, except that `time2.m` still serves for timing the prerequisite.
`NOTES.md` is the working log.
