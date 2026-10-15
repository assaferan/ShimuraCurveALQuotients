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

1. ~~The Atkin–Lehner layer.~~ **Done**, see below.
2. **Agreement with the models we already have.** Still open, and it is the step that matters; see
   "What is still NOT done" below for why it is more than a thin wrapper and where to start.
3. **Beyond the Borcherds reach.** The timings below say this is now worth attempting: the first
   bases behind the divisor wall compute their domains in seconds to minutes.

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

## What is still NOT done: agreement with the stored models

This is the step that would make the method a second source, and it is not built. What exists is a
dimension match, which is not a model match.

The missing piece is model reconstruction: functions on the quotient as ratios of forms of equal
weight, their degrees, and the algebraic relation between two of them fitted from values. Two
things make it more than a thin wrapper:

* **At D = 6 the comparison has no content.** The star curve has three elliptic points, and any
  three points are equivalent under a Mobius change of coordinate, so matching the stored model's
  branch points is free. The branch assignment is already determined by the signatures, which come
  from the repository rather than from the expansions. The first base where a coordinate-free
  comparison bites is D = 10, whose star curve has four elliptic points (orders 2, 2, 2, 3) and
  hence a genuine cross-ratio to compare.
* The stored models are written over a Hauptmodul of the **star** curve, so a comparison needs that
  Hauptmodul reconstructed and normalised the same way, and at D = 10 the star curve first reaches
  dimension 2 at **weight 12**, which at 164 terms is a much heavier computation than weight 4.

Coordinate-free invariants are the way to do the comparison when it is built: the j-invariant for a
genus-1 quotient, Igusa invariants for genus 2, the cross-ratio of the elliptic points for genus 0.
Comparing raw coefficients would require recognising numerical values as rationals and matching a
coordinate, which is avoidable.

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

## Files

`vw3.m` is the working library; `dims.m` computes dimensions of M_k, `haupt.m` builds a Hauptmodul
as a ratio of weight-12 forms, `cmchk.m` evaluates it at a CM point, `predict.m` is the
valence-formula prediction. `vw2.m`, `run2.m`, `shape*.m` and `dbg_*.m` are development history,
kept only because they record what was tried. `vwlib.m` and `NOTES-firstattempt.md` are from the
first attempt and are superseded, except that `time2.m` still serves for timing the prerequisite.
`NOTES.md` is the working log.
