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

**D = 15, 22 and 26 now work too** (2026-10-15), at weight 4 and the correct dimensions 2, 3 and 3,
with gaps of 10^9 at a truncation target of 1e-12 and 10^17 at 1e-20. The cause of their failure
was a defective reduction, not the size of the domain; see "What actually blocked the larger bases"
below, which supersedes the retention argument in the section before it. Hecke relations
(eq. 18) are now implemented and validated as well. What is left is not resolution but cost, and
that is quantified at the end of this file.

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
   ~~D = 14 and D = 15 need the side-pairing sampling described below first.~~ Done: D = 15, 22
   and 26 resolve with the side-pairing reduction.
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

⇒ **Superseded in its conclusions by "What actually blocked the larger bases" below.**  The
retention figures quoted here were measuring a broken reduction; with it fixed, retention is 100
percent at every base and two rows per term suffice everywhere.  The section is kept because the
rule it states -- size the sampling by the rows RETAINED, not by the points attempted -- is still
right, and because it records what the numbers looked like before.

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

But over-determination is necessary and **not sufficient**, and the reading that followed from
this table -- that the domain radius is a hard wall somewhere near 0.693 -- was **WRONG**.  It came
from testing every base at the same truncation target, 1e-12, which is far too loose once rho is
large.  See the next section.

## ⇒ CORRECTION (2026-10-15): the radius is not a wall.  The term count was too low

**D = 22, at rho = 0.839, resolves.**  Holding everything else fixed and raising only the
truncation target (and the working precision with it), `walldiag.m` gives, for weight 4 where the
expected dimension is 3:

| target | NT | prec | rows | floor rho^NT | dimension found | gap |
|---|---|---|---|---|---|---|
| 1e-12 | 158 | 120 | 297 | 9.1e-13 | 4, **wrong** | 14.8 |
| 1e-20 | 263 | 150 | 494 | 9.1e-21 | **3, right** | 63.1 |
| 1e-28 | 368 | 170 | 692 | 9.1e-29 | **3, right** | 620 |

The three kernel singular values sit at 0.9 to 8 times the floor at every term count, while the
first non-kernel value sits at 62, then 503, then 2176 times it.  That is the signature of
TRUNCATION, not of conditioning: the kernel tracks the floor downward while the gap above it widens
by an order of magnitude per step.  Conditioning would instead show the whole spectrum sliding down
together with no gap ever opening, which is exactly what the mis-centred D = 10 computation did
before `RecentredDisc` (see above), and the contrast between the two is what settles it.

**So the rule is not a radius cutoff but a relation between the two.**  The term count must be set
so that the first non-kernel singular value stands clear of rho^NT, and since the kernel sits AT the
floor, the way to check a dimension is to raise NT and confirm the gap grows.  A single fixed
target, as the drivers used, silently under-resolves every base with a large radius.

**There are TWO distinct failure modes, and the radius does not separate them.**  D = 15 and
D = 26 do not resolve however many terms they are given, and both have a SMALLER radius than D = 22,
which does.  So the monotone-in-rho reading of the previous table is wrong twice over.

| base | rho | NT | rows | floor | smallest sv | verdict |
|---|---|---|---|---|---|---|
| D = 22 | **0.839** | 368 | 692 | 9.1e-29 | 8.2e-29, i.e. **at** the floor | resolves |
| D = 26 | 0.774 | 252 | 1026 | 1.0e-28 | 3.6e-38, **10 orders below** | fails, got 11 for 3 |
| D = 15 | 0.736 | 211 | 828 | 7.7e-29 | 1.2e-74, **45 orders below** | fails, got 1 for 2 |

Where the smallest singular values sit AT the floor, the limit is truncation and more terms fix it.
Where they sit far BELOW the floor, the matrix is rank-deficient for reasons that have nothing to do
with the series: the condition numbers are 9.0e37 at D = 26 and 2.8e74 at D = 15, against 1.4e12 at
D = 22.  That one measurement -- the smallest singular value against rho^NT -- tells the two apart,
and it is what `walldiag.m` exists to report.

**What correlates with it is how much of the sampling circle lies OUTSIDE the fundamental domain**,
as the retention of `ReducedSamples` measures: 63 percent at D = 22, 27 percent at D = 26, 14
percent at D = 15, in exactly the order of the table.  A point inside the domain reduces to itself
and carries no relation, so when the domain covers most of the circle the relations come from a few
short arcs, and recovering two hundred coefficients from short arcs is exponentially
ill-conditioned however many samples are taken from them.

⇒ **The mechanism is right and the numbers are an artefact.  See the next section**: those
retention figures were measuring a broken reduction, not the geometry, and the remedy is not a
larger sampling radius.

**Remark 3.3 is NOT the explanation, contrary to a natural first guess.**  The paper recommends
working with f_norm = f (Im z)^{k/2}, whose automorphy factor has absolute value 1, "if rho is
large".  Measured (`jfac.m`), |j|^k is bounded by 1 at every base and its spread across samples is
only 1.5 at D = 6 rising to 11 at D = 22 -- nothing like the seven orders of magnitude that would
be needed to explain the failures.

## ⇒ What actually blocked the larger bases (2026-10-15): the reduction, in one line of `vw3.m`

A circle of radius rho touches the fundamental domain only at its outermost vertices, so **the
true retention at radius rho is essentially 100 percent at every base**, and on a circle outside
the domain it is 100 percent exactly.  The 14 and 27 percent above were not geometry.  They were
`ShimuraReduceUnit` failing.

`ReducedSamples` called it with `Generators(G)`, which is a **minimal generating set** -- three
elements where the domain at D = 15 has fourteen sides.  Given too few generators it returns the
identity on points it cannot reduce, with no indication of failure, and `ReducedSamples` discarded
those points on the correct-looking test that a trivial reduction carries no relation.  So the
points that survived were the few the three generators happened to handle, and they cluster on
short arcs.  The ill-conditioning the previous section describes was real; its cause was this.

The fix is to descend over the **side pairing** of the domain, which Magma records in
`G`ShimFDGens` (fourteen elements at D = 15, one per side), repeatedly applying the pairing that
most decreases |w|.  The domain Magma returns here behaves as the Dirichlet domain centred at the disc centre, so a
point should be reduced exactly when no side pairing decreases its modulus; that the descent does
terminate inside was checked directly, on the Hecke images at D = 6, 15, 22 and 26, where the
worst modulus after descent is below rho at all four and was 0.91 to 0.99 before.  `vw3.m` now
does this; `ShimuraReduceUnit` is no longer used for sampling.

Measured at weight 4, truncation target 1e-20, two rows per term, recentred disc:

| base | volume | rho | retention before | retention after | condition before | condition after | dimension | gap |
|---|---|---|---|---|---|---|---|---|
| X_0^6(1) | 0.33 | 0.447 | 53% | **100%** | 1.3e21 | 1.5e18 | 1 of 1 | 6.9e17 |
| X_0^10(1) | 0.67 | 0.648 | — | **100%** | — | 3.1e18 | 1 of 1 | 8.9e17 |
| X_0^15(1) | 1.33 | 0.736 | 15% | **100%** | 2.9e53 | 4.5e18 | **2 of 2** | 2.3e17 |
| X_0^22(1) | 1.67 | 0.839 | 63% | **100%** | — | 4.5e18 | 3 of 3 | 7.5e16 |
| X_0^26(1) | 2.00 | 0.774 | 27% | **100%** | — | 3.2e18 | **3 of 3** | 1.6e17 |

D = 15 is the clearest case: the condition number falls from 2.9e53 to 4.5e18 and the gap rises
from 59 to 2.3e17.  The Elkies comparison is unchanged by the fix -- 14 of 14 CM values at D = 6
weight 24 and 12 of 12 at D = 10 weight 12, exactly as before.

**Consequences for the two knobs that were expected to matter.**

* *More terms.*  Worth nothing on a base that is rank-deficient.  Before the fix, D = 15 at 91,
  151 and 211 terms gave gaps of 22.6, 58.8 and 40.4 while the smallest singular value ran from 20
  to 47 orders of magnitude BELOW the truncation floor.  After the fix the same base gives 5.5e9,
  2.3e17 and 1.5e25: the gap grows by eight orders per step, which is what a sound system does.
  So raising the term count is how one confirms a dimension, not how one rescues a base.
* *Sampling outside the domain.*  Now buys nothing, because retention is already 100 percent at
  radius rho, and it costs a great deal: the series is evaluated out to the sampling radius, so
  the term count follows that rather than rho, and at D = 15 a contour at 1.25 rho needs 331 terms
  where rho needs 91.  Useful only in the Cauchy formulation, where the contour has a different
  job.

**A reporting defect found along the way.**  `GapDimRel` returns the index of the largest relative
gap and therefore always names a dimension, even on a spectrum that is a smooth continuum with no
gap at all.  Every "wrong dimension" recorded for D = 15 and D = 26 was really "no resolution":
the largest ratio in those spectra was about 50.  `GapDimStrict` refuses to answer unless the gap
reaches a stated size and the first value above it is a real singular value rather than more
noise.  Use it; the paper calls a gap of 4 unsatisfactory and a gap of 1e19 convincing.

## Hecke relations (eq. 18): implemented and validated

`vwhecke.m`, with `heckecal.m` to calibrate and `heckedim.m` to apply.

The paper writes (18) with the upper half-plane automorphy factor, and everything else here works
in the disc.  Writing w = (z-p)/(z-pbar), so that 1-w = (p-pbar)/(z-pbar), and J(g,w) for the
factor of the SU(1,1) matrix, one has the identity

    J(g,w) = j(g,z) (1-w)/(1-w')

and the three half-plane factors in (18) collapse: the (1-w) weights cancel and

    K_n^h = sum_i [ J(pi_i, 0) J(g_i, w_{p,i}) ]^-k (w'_{p,i})^n ,   sum_n K_n^h b_n = a_p b_0 .

The same identity is why the rows in `vw3.m` are the bare `w'^n - J^k w^n` with no (1-w)^k weights:
that is eq. (11) divided by -(1-w')^k.

**The coset representatives need no group cohomology.**  The p+1 left ideals of reduced norm p are
`O*x + O*p` for x in O of reduced norm p, and such an x generates its own ideal, so a box search
over a basis of O with deduplication by ideal equality gives the cosets.  A box of 5 gives all p+1 at
D = 6 and D = 15 for p up to 13, but only 8 of 14 at D = 22 for p = 13 and 2 of 6 at D = 26 for
p = 5.  A cube in the coordinates of a basis is the wrong shape for a skew order, so this wants
replacing by an enumeration of short vectors of the positive definite form given by the Frobenius
norm of the matrix representation, whose determinant is the reduced norm.  Until then, check the
count against p+1 and skip the prime if it is short -- `heckedim.m` does.

**The normalisation, fixed by measurement and confirmed by derivation.**  With the coset matrices
scaled to determinant 1, as `DiscMatrix` does, the operator of (17) has eigenvalue
p^(1-k/2) a_p, where a_p is the classical one.  The check is at D = 6 weight 4, where the space is
one-dimensional and so the computed expansion is itself an eigenform: the measured ratio comes out
as a_p/p at p = 5, 7, 11 and 13, to the full 149 digits of the computation and with an imaginary
part of the same size.  Four primes agreeing exactly is what distinguishes a normalisation from a
fitted constant.

**It does what the paper says it does.**  At D = 15 weight 4 the space is two-dimensional; adding
one Hecke row per prime for p = 2 and 7, with the eigenvalues of each of the two newforms of
S_4^{15-new}(Gamma_0(15)), splits the kernel into two one-dimensional ones with gaps of 6.4e9 and
6.7e9.  **With a deliberately wrong eigenvalue the kernel disappears**, the smallest singular value
rising by a factor of 1e9; that negative control is what shows the rows are constraining something.

Two things worth knowing before relying on it:

* The accuracy of a Hecke row is set by truncation at the WORST reduced Hecke image, which can sit
  much further out than the samples do.  Before the reduction was fixed the images at D = 15 landed
  at modulus 0.98 against rho = 0.736 and the rows were pure noise, changing nothing.
* A Hecke relation cannot rescue a rank-deficient system.  Applied to the pre-fix D = 15 matrix it
  moved the smallest singular value by a factor of 1.0.  It reduces the gap a base must achieve,
  not the conditioning it must have.

## One thing the paper prescribes that is still not implemented

It is in Section 3 and bears directly on large rho.

1. **The Cauchy-integral formulation** (eq. 13-14 with Simpson's rule at eq. 16).  The paper is
   explicit that it prefers this: relations are assembled "by taking appropriate linear
   combinations of these relations, or (what seems better in practice) using the Cauchy integral
   formula".  All three of its Shimura-curve examples use it; the automorphy relations of eq.
   (9)-(11), which `vw3.m` implements, are presented first and then superseded.  The matrix is
   SQUARE, (N+1) x (N+1), where ours is overdetermined, and Remark 3.4's two rescalings -- b'_n =
   b_n rho^n, then dividing row n by (K_nn - 1) -- are what make it nearly diagonally dominant.
   `vwcauchy.m` is a prototype.  It reproduces the correct dimension at D = 10 with fewer terms (64
   against 107) and reaches 1e-28 where the automorphy system reaches 1e-12, but at D = 15, 22 and
   26 it returns large numbers of spurious near-zero singular values and does not yet match the
   paper's conditioning.  It is not debugged and should not be trusted as it stands.

   ⇒ Its spurious near-zero singular values at D = 15, 22 and 26 were measured BEFORE the
   reduction was fixed, and that fix is the likeliest cause; it has not been re-measured since.
   Re-measure before treating it as broken.

**For the record, the radii the paper itself reaches**: Example 2 (D = 6) rho = 0.447213, which
matches our recentred value exactly; Example 3 rho = 0.71807 with N = 150, one singular value below
1e-20 and the next at 0.4994; Example 4 rho = 0.885611 with N = 600.  Its Example 1 is the
non-cocompact case at rho = 0.85, where the paper says plainly that "our method does not work well
in this non-cocompact situation" and reports a ramp, the next singular value being "only 4 times as
large".

## An alternative worth knowing about

Nelson, *Evaluating modular forms on Shimura curves*, Math. Comp. 84 (2015) 2471-2503
(arXiv:1210.1243), evaluates a newform by computing the Shimizu lift, a realisation of
Jacquet-Langlands, as "an explicit, rapidly-convergent formula that is well-suited for numerical
computation".  The Voight-Willis paper cites it and says "It will be interesting to compare the two
techniques."  Two things to weigh before reaching for it: its input is a newform specified by its
Hecke eigenvalues and it evaluates THAT form, so it gives one eigenform at a time rather than a
basis of a space; and implementing the theta correspondence effectively is a substantial piece of
work that Magma does not provide.  For our purpose -- a basis of weight-2 forms on a quotient, to
feed a canonical model -- the eigenform-at-a-time shape is workable, since the repository can
supply the eigenvalues, but it is a bigger build than finishing Section 3.

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

## ⇒ How far this reaches, now that resolution is not the obstacle (2026-10-15)

⇒ **THE COST CONSTANT IN THIS SECTION IS WRONG BY A FACTOR OF ABOUT SEVEN HUNDRED, and the table
of timings and the verdict that follow from it are wrong with it.  See "What the decomposition
actually costs" and "Expand on the normaliser's domain" below, both 2026-10-15 (later).**  The
shape of the argument -- that the cost is set by rho alone, through NT = T ln10 / (-ln rho), and
grows as NT^3 -- is right, and so is the observation that the working precision need not grow
with NT.  Only the constant, and therefore every absolute time and every conclusion about what is
affordable, is wrong.

With the reduction fixed, nothing in the method fails at the bases tested.  What stops it is cost,
and the cost is governed by **one number per base: the circumradius rho of the recentred domain**.

Everything follows from it.  The term count is NT = T ln10 / (-ln rho) for a truncation target of
1e-T.  The singular value decomposition is then the whole expense and it goes as NT^3, measured at
3.1e-8 seconds per NT^3 at working precision 200 for a system of 2 NT rows -- consistent to within
15 percent across D = 15, 22 and 26.  The working precision does **not** have to grow with NT,
because the dynamic range of the matrix is rho^NT, which the target pins to 1e-T by construction.
Memory for the real form of the matrix is 8 NT^2 complex entries, about 332 bytes each.

Measured by `rhoprobe.m`, at a target of 1e-20:

| base | volume | domain | sides | rho | NT | decomposition | memory |
|---|---|---|---|---|---|---|---|
| X_0^6(1) | 0.33 | 0.5 s | 6 | 0.4472 | 58 | < 1 s | 0.01 GB |
| X_0^10(1) | 0.67 | 0.3 s | 6 | 0.6483 | 107 | < 1 s | 0.03 GB |
| X_0^15(1) | 1.33 | 0.3 s | 14 | 0.7358 | 151 | < 1 s | 0.06 GB |
| X_0^22(1) | 1.67 | 1.6 s | 10 | 0.8391 | 263 | 1 s | 0.17 GB |
| X_0^26(1) | 2.00 | 1.8 s | 14 | 0.7742 | 180 | < 1 s | 0.08 GB |
| X_0^35(1) | 4.00 | 5.4 s | 26 | 0.8788 | 357 | 1 s | 0.32 GB |
| X_0^51(1) | 5.33 | 23 s | 38 | 0.9484 | 870 | 20 s | 1.9 GB |
| X_0^58(1) | 4.67 | 14 s | 22 | 0.9692 | 1474 | 99 s | 5.4 GB |
| X_0^330(1) | 13.3 | 29 s | 70 | 0.9852 | 3100 | 15 min | 24 GB |
| X_0^6(35) | 16.0 | 267 s | 102 | 0.9883 | 3938 | 32 min | 38 GB |
| X_0^462(1) | 20.0 | 886 s | 100 | 0.9943 | 8117 | 4.6 h | 163 GB |
| X_0^510(1) | 21.3 | 236 s | 116 | **0.9977** | **19944** | **68 h** | **984 GB** |

**Volume is only a rough proxy and rho must be measured per base**: X_0^58(1) has a smaller volume
than X_0^51(1) and a larger rho, and X_0^510(1) needs two and a half times the terms of
X_0^462(1) on six percent more volume.  The geometric lower bound rho >= sqrt(V/(V+2)), from the fact
that a domain of hyperbolic area 2 pi V cannot fit inside |w| < r unless r^2/(1-r^2) >= V/2, is
weak -- it gives 0.957 at X_0^510(1) where the truth is 0.9977, and the term count depends on the
difference through 1/(-ln rho), a factor of 19 here.

**So the answer for the 34 open quotients is no, and by a wide margin.**  Their bases have
arithmetic volumes of 46 to 240.  At volume 21 the decomposition already needs about a terabyte,
and the memory grows as NT^2 while rho is still climbing towards 1.  Lowering the target to 1e-12
scales NT by 0.6 and the memory by 0.36, which does not change the verdict.

Neither of the remaining levers touches this.  The Cauchy formulation makes the matrix square,
(N+1) x (N+1) instead of 2N x N, which saves a factor of two in rows and nothing in the N^2
memory.  Hecke relations reduce the gap a base has to achieve, not the number of coefficients it
has to carry.  **The wall is NT ~ 1/(-ln rho), and it is structural**: a base whose domain nearly
fills the disc needs that many coefficients to describe a form on it, whatever one does with them
afterwards.  Reaching volume 46 and above needs a method whose cost is not a dense decomposition
of a matrix of that size -- Nelson's Shimizu lift (below) evaluates one eigenform at a time and
has no such matrix, which is the shape worth weighing.

What this does open, since it was not reachable before: everything up to about volume 16 is now
routine on one machine, which includes **X_0^330(1)** at 15 minutes and 24 GB and **X_0^6(35)** at
32 minutes and 38 GB, both behind the Borcherds divisor wall.  Volume 20 and 21 are a large
server, and the domain itself starts to dominate there -- 886 seconds at X_0^462(1).

## ⇒ What the decomposition actually costs (2026-10-15, later)

The constant above was never measured; it was an estimate written into `rhoprobe.m`, and the
"measured" claim beside it is an error.  Timed directly, Magma's `SingularValueDecomposition` on a
real m x n matrix with m = 2n takes

| n | 80 | 160 | 320 |
|---|---|---|---|
| precision 40 | 1.4 s | 8.4 s | 56 s |
| precision 70 | 1.5 s | 9.2 s | 60 s |
| precision 120 | 1.9 s | 11.2 s | 71 s |

so the cost grows as n^2.8 and is **almost independent of the working precision** -- going from 40
to 120 digits costs a third more, not three times more.  The system here has NT complex columns
and about 2 NT complex rows, whose real form is 4 NT x 2 NT, so n = 2 NT, and the timing agrees
with a run of `heckedim.m` at D = 15: NT = 151, 302 rows, **79 seconds**.  The rule is therefore

    one decomposition  ~  79 s x (NT/151)^2.8

against the 0.11 s the old constant predicts there.  **Everything in the reach table above is
about seven hundred times too fast**: X_0^330(1) at NT = 3100 is not 15 minutes but something over
a hundred hours, and a quotient of the cheapest of the open bases would have been nineteen hours
rather than a minute.  Precision is not a lever, and neither is the normal-equations form of the
decomposition (forming M^t M and decomposing that is 38 s against 60 s at n = 320, a factor of
1.6).  The lever is the next section.

## ⇒ Expand on the NORMALISER's domain, not on Gamma's (2026-10-15, later)

Every model this project stores is of a quotient by an Atkin-Lehner subgroup, so the forms wanted
are invariant under the group generated by Gamma and the Atkin-Lehner elements.  That group's
Dirichlet domain is 1/#W* of the area and has a smaller circumradius, and since NT goes as
1/(-ln rho) with rho close to 1, a modest reduction in rho is a large reduction in NT and a cubic
one in the cost.  Measured by `rhostar.m`, which reduces several thousand points of the disc over
the larger group and takes the largest modulus reached:

| base | volume | #W* | rho(Gamma) | rho(normaliser) | NT | NT* | cost ratio |
|---|---|---|---|---|---|---|---|
| X_0^6(1) | 0.33 | 4 | 0.4472 | 0.4185 | 58 | 53 | 1.3 |
| X_0^35(1) | 4.0 | 4 | 0.8788 | 0.7584 | 357 | 180 | 6.9 |
| X_0^10(9) | 8.0 | 8 | 0.9515 | 0.7300 | 927 | 147 | 180 |
| X_0^6(25) | 10.0 | 8 | 0.9615 | 0.8094 | 1174 | 218 | 115 |
| X_0^15(7) | 10.67 | 8 | 0.9810 | 0.9020 | 2406 | 447 | 115 |
| X_0^85(1) | 10.67 | 4 | 0.9656 | 0.9442 | 1316 | 803 | 4.0 |
| X_0^142(1) | 11.67 | 4 | 0.9825 | 0.9194 | 2608 | 548 | 81 |
| X_0^14(9) | 12.0 | 8 | 0.9820 | 0.8489 | 2536 | 282 | 490 |
| X_0^21(4) | 12.0 | 8 | 0.9761 | 0.8652 | 1905 | 319 | 154 |
| X_0^21(5) | 12.0 | 8 | 0.9861 | 0.8627 | 3286 | 312 | 765 |
| X_0^91(1) | 12.0 | 4 | 0.9847 | 0.9630 | 2991 | 1220 | 13 |
| X_0^158(1) | 13.0 | 4 | 0.9877 | 0.8826 | 3731 | 369 | 682 |
| X_0^123(1) | 13.33 | 4 | 0.9751 | 0.9535 | 1826 | 968 | 6.0 |
| X_0^166(1) | 13.67 | 4 | 0.9926 | 0.9706 | 6170 | 1544 | 50 |
| X_0^141(1) | 15.33 | 4 | 0.9903 | 0.9647 | 4737 | 1282 | 40 |
| X_0^15(8) | 16.0 | 8 | 0.9846 | 0.9160 | 2954 | 525 | 131 |

X_0^15(11) and X_0^35(3), both of volume 16, are not in the table: their fundamental domains had
not finished after an hour and a half, and the domain is the prerequisite for everything else.

The gain is not predictable from #W* -- X_0^21(5) and X_0^158(1) gain a factor of several hundred
while X_0^85(1) with the same volume gains four -- so it has to be measured per base, which takes
a minute or two and needs no decomposition.  Together with the corrected constant it is the
difference between these bases being unaffordable and being an hour's work each.

⚠ **At non-squarefree level the figure depends on the Atkin-Lehner elements being right.**  The
rows for X_0^10(9), X_0^14(9), X_0^21(4) and X_0^15(8) above are the ones measured AFTER the fix
described next; before it, the search returned a scalar for the Atkin-Lehner element of a square
Hall divisor, which acts trivially, so the group was too small and the radius too large -- 0.837
instead of 0.730 at X_0^10(9), a factor of five in the cost and a weight-4 star space of
dimension 4 where the signature says 3.  **A scalar is not an Atkin-Lehner element**: 3 has
reduced norm 9 and normalises the order at level 9, and nothing but an explicit test rules it
out.  `alchk.m` is that test, over every Hall divisor of a base.

The reduction descends over the side pairing of Gamma and over the Atkin-Lehner translates
alternately, and records which coset of the normaliser it ended in as a Hall divisor.  That label
is what lets ONE set of samples serve three spaces at once: a form invariant under W with
eigenvalue eps on the other coset satisfies f|g = eps f, so the invariant space, the star space
and the anti-invariant space differ only in the table of signs.

## ⇒ From the expansions to an equation over Q (`model.m`, `finish.m`)

For W of index 2 in the full group W*, with the star curve of genus 0, the model is a double cover
y^2 = f(t) of the star written over a hauptmodul t -- which is the form this project stores.  Two
functions come out of the expansions:

* **t = f1/f2**, a ratio of two forms on the star.  A space of weight-k forms on a genus-0
  quotient is the sections of a line bundle of degree dim - 1, so a ratio of two of them has
  degree dim - 1 and a hauptmodul needs dim - 2 common zeros imposed; each extra zero is one more
  row asking the series to vanish at a point of the disc.
* **y = g/h**, with g changing sign under the deck transformation and h a form on the star of the
  same weight.  y^2 is invariant under all of W*, hence a rational function of t, and fitting
  P(t) - y^2 Q(t) = 0 at points of the disc recovers it.  The branch points are the roots of P*Q
  of **odd** multiplicity, together with infinity when the degrees differ by an odd number, and
  there must be exactly 2g+2 of them -- an internal check the fit has to pass.

The coordinate the computation chose for itself is not defined over Q, so the equation cannot be
read off in it.  `finish.m` goes through the **invariants**, which do not depend on the coordinate
at all and are rational for a curve over Q: the absolute Igusa-Clebsch invariants in genus 2, the
j-invariant in genus 1.  Recognise them as rationals, rebuild the curve over Q from them (Mestre's
construction, which Magma provides), and what is left is exactly the **quadratic twist**, which is
arithmetic and which a complex-analytic method cannot see.

**The twist is solved for, not searched for.**  Mestre's construction returns an arbitrary member
of the isomorphism class over Qbar, and the number relating it to the curve over Q is not small:
at X_0^85(1) it is 288862917341177298307577.  But twisting a genus-2 curve by c replaces a_p by
-a_p, so comparing the point count of the reconstructed curve with the trace formula at a prime p
reads off the Legendre symbol (c|p); and c can only be divisible by primes dividing the
discriminant of the reconstructed curve or 2DN.  Writing c as a product of those with exponents 0
or 1 turns each prime into one linear equation over F_2, and the system is heavily
overdetermined.  At X_0^85(1): 11 candidate primes, 25 usable symbols, a unique solution, and the
resulting curve then matches 64 point counts over F_p and F_{p^2} that it was not fitted to.

**Validated end to end on the arithmetic stage.**  Starting from the six branch points of Molina's
X_0^85(1)/w_5 carried through a random Mobius change and perturbed at the 1e-28 level -- which is
what a truncation target of 1e-28 leaves -- `finish.m` recovers the three absolute invariants
exactly (to a relative 1.3e-24), solves the twist exactly, and returns
y^2 + (x^3+1) y = -x^6 - 4x^5 + 5x^4 + 12x^3 - 9x^2 - 5x - 2, isomorphic over Q to the published
equation, with 64 point counts against the trace formula.  **That is the right control for this
stage**: it tests the step from numerical branch points to a verified rational model, with the
geometry supplied by a published answer rather than by the expansions.

**How many digits are needed.**  Two different numbers, and conflating them wastes a run.  The
WORKING precision must be comfortably above 80 or `IgusaClebschInvariants` over the complex
numbers is meaningless -- at precision 40 it returns I10 wrong in every digit, which is
cancellation inside the invariant and nothing to do with the input.  What a recognition may rely
on is the TRUNCATION target, less about four digits for the degree-10 arithmetic.  The invariants
of these curves have numerators of seventeen digits, so a target of 1e-28 is about the minimum.

**Five of the sixteen open bases are out of scope for this construction**, because their star
curve has genus 1 and there is no hauptmodul to write over: X_0^142(1), X_0^91(1), X_0^158(1),
X_0^123(1) and X_0^141(1).  A quotient of index 4 is likewise not a double cover of the star; it
is the fibre product of the three index-2 quotients between, which is how the rest of the
repository builds them.

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

* **A CUBE IN THE COORDINATES OF A BASIS IS THE WRONG SHAPE FOR AN ORDER.**  The reduced norm
  form is indefinite of signature (2,2), so its level sets are hyperboloids and a cube cuts them
  along directions unrelated to the lattice.  Three searches here did that and all three missed
  elements silently: the Hecke cosets (8 of 14 at D = 22 for p = 13, 2 of 6 at D = 26 for p = 5),
  the Atkin-Lehner elements (nothing at all at D = 35 for m = 7), and the CM points.  Enumerate
  short vectors of the positive definite majorant `FrobeniusGram` instead, raising the bound until
  the expected number appears -- and check the number against what theory says it is, which for
  the Hecke cosets is p+1.
* **CAPPING A ROW BUDGET WITH A BREAK KEEPS A SHORT ARC.**  Samples come in order of angle around
  the sampling circle, so taking the first 2 NT of several thousand takes a sixth of the circle,
  and fitting a degree-NT polynomial on a short arc is exponentially ill-conditioned however many
  points it has.  Take them with a STRIDE.  The symptom is worth recognising because it is the
  same one the broken reduction produced: a smooth geometric ramp of singular values, running
  far BELOW the truncation floor, with no gap anywhere.  Separately, each row is normalised to
  unit length so that every relation carries the same relative truncation error whatever its
  automorphy factor; that was not the cause here (the factors span less than one order of
  magnitude at D = 35) but it costs nothing.
* **`by` IS A RESERVED WORD IN MAGMA** (`for i in [1..n by 2]`), so it cannot be a variable name;
  using it gives "bad syntax" at the line that assigns it and nowhere else.
* **MAGMA READS A SCRIPT AS IT RUNS IT.**  Editing a driver while a run of it is in flight makes
  that run fail with a syntax error in the middle of a finished computation -- it happened twice
  here, each time throwing away a decomposition that had already succeeded.  Same hazard as
  editing `tests/` during a suite.  Write progress to a separate file that is closed after every
  line (Magma buffers stdout to a file, so the log tells you nothing), and leave the script alone
  until the run ends.
* **`Generators(G)` is NOT the side pairing of the fundamental domain.** It is a minimal
  generating set: three elements at D = 15, where the domain has fourteen sides. The side pairing
  is `G`ShimFDGens`, for the domain Magma most recently computed, so recompute the domain for the
  disc in hand before reading it. Its entries are order elements (`AlgQuatOrdElt`), not group
  elements, so `Quaternion` fails on them; coerce into the algebra instead.
* **`ShimuraReduceUnit` returns the identity on a point it cannot reduce**, with no indication of
  failure, so "the reduction was trivial" and "the reduction failed" are indistinguishable at the
  call site. Given too few generators it fails often: at D = 15 on 85 percent of a circle of
  radius rho. Check the modulus of what comes back against rho rather than trusting it.
* **`GapDimRel` always names a dimension**, being the index of the largest relative gap, so it
  reports a number on a spectrum with no gap at all. Every "wrong dimension" at D = 15 and D = 26
  was really "no resolution". Use `GapDimStrict`.

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

`vw3.m` is the working library (holding the sampling, the side-pairing reduction, the automorphy
system and the positive definite majorant `FrobeniusGram` that the three element searches share),
`vwal.m` the Atkin-Lehner layer (which also holds `RecentredDisc`, `ZeroLevel`, `GapDimStrict`,
and the reduction over the normaliser with its coset labels) and `vwhecke.m` the Hecke relations
of eq. (18).

For building an equation: **`model.m`** takes a base and one Atkin-Lehner subgroup of index 2 and
writes the branch points of the double cover over the star to `out/`, with `vwmodel.m` holding the
labelled system, the imposed-zero rows and the rational-function fit; **`finish.m`** reads that
file and returns the curve over Q with its twist settled against the trace formula.  `rhostar.m`
measures how much the normaliser's domain saves at a base, in a minute and with no decomposition,
and is the first thing to run at a new base.  `heckecos.m` checks the Hecke coset search finds all
p+1.  For measurement: `reach.m` asks
whether a base resolves and at what cost, sweeping truncation targets and sampling radius, and is
the first thing to run at a new base; `rhoprobe.m` gives the circumradius and the predicted cost
without computing anything else; `heckecal.m` calibrates the Hecke normalisation against a known
eigenvalue and `heckedim.m` applies the relations, with a wrong-eigenvalue control.
`aldims.m` is the validation sweep over all Atkin-Lehner
subgroups of a base, `basediag.m` tests the base library alone at one base, `recentre.m` reports
how much recentring gains, `aldiag.m`
instruments the reduction quality; `dims.m` computes dimensions of M_k, `haupt.m` builds a Hauptmodul
as a ratio of weight-12 forms, `cmchk.m` evaluates it at a CM point, `predict.m` is the
valence-formula prediction. `vw2.m`, `run2.m`, `shape*.m` and `dbg_*.m` are development history,
kept only because they record what was tried. `vwlib.m` and `NOTES-firstattempt.md` are from the
first attempt and are superseded, except that `time2.m` still serves for timing the prerequisite.
`NOTES.md` is the working log.
