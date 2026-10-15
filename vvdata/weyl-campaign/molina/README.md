# Equations of hyperelliptic Shimura curves by specialisation of Heegner points (Molina)

A third way to produce **equations for Shimura curves and their Atkin–Lehner quotients**, next to
the Borcherds/Guo–Yang route the rest of the repository uses and the Voight–Willis power series in
the neighbouring directory.

Reference: S. Molina, *Equations of hyperelliptic Shimura curves*, arXiv:1004.3675v1 (2010);
PDF in Dropbox at `/Eran Assaf/Research/HyperellipticQuotients/Molina.pdf`. Sections 7 and 10 are
the algorithm, Sections 8 and 9 are worked examples at discriminants 39 and 55, Section 10.3 is the
worked quotient at 35, and Section 10.4 Table 2 lists the ten equations he obtained.

Where the Borcherds route builds a model from an eta-quotient basis, this one never leaves the
arithmetic of quaternion orders. A hyperelliptic quotient has its Weierstrass points at Heegner
points, so `y^2 = R(x)` factors as `prod p_{R_i}(x)`, one factor per Galois orbit of Heegner points.
Cerednik–Drinfeld at a prime `p | D` says where those points go in the special fibre, and the
thickness of the singular point two of them meet at is the valuation of the difference of the
corresponding roots of `R`. That pins the discriminants and resultants of the factors; the Gross–
Zagier formula for endomorphism rings over Artinian rings pins the leading coefficients; class
field theory pins the number field each factor defines. What is left is to find integer polynomials
with a prescribed field, leading coefficient and discriminant, which is an index form equation.

The three methods are limited by entirely different things:

| | Borcherds / Guo–Yang | Voight–Willis | Molina |
|---|---|---|---|
| what blocks it | the number of divisors of M = 2DN (the Normaliz wall at 24), and the m = 0 multiplier when N > 1 | the size of the fundamental domain | **the degree of the number field a factor defines** — `IndexFormEquation` is implemented for degree 3 and 4 only |
| needs D odd | no | no | **yes** (good reduction at 2 is used throughout) |
| reaches non-hyperelliptic curves | yes | yes | **no** |
| exactness | exact | numerical | exact |

So it complements the Borcherds route exactly where that route is worst — odd discriminant, where
the eta-quotient pool is sparse — and it is useless where the Borcherds route is strongest.

## Status

**All ten equations of Molina's Table 2 are reproduced**, each checked isomorphic over Q to the
published one and each checked against the trace formula. Two of the ten —
`X_0^85(1)/<w_5>` and `X_0^85(1)/<w_85>` — are quotients this repository holds no model for, so they
are new to us, and **two further ones are new outright**: `X_0^115(1)/<w_115>` and
`X_0^161(1)/<w_161>`, both of genus 2.

Every equation in `equations.m` passes, in `verify.m`:

* the genus from `GenusShimuraCurveQuotient`, and the genus of the curve itself;
* the number of points over `F_p` and over `F_{p^2}` against `ComputePointsViaTrace`, at every odd
  prime of good reduction below 60 — **336 counts, no disagreement**.

Neither check knows anything about the construction.

## How far it reaches, measured

Of the 43 hyperelliptic Atkin–Lehner quotients with odd discriminant, level one and genus at least
2 that have no model anywhere in this repository (`census.m` builds that list from a merged copy of
every branch's `data/models`), **exactly 5 are inside the method** (`scope.m` decides this from
Step 1 alone, in a minute, without building a single ring class field):

| why it is out | how many |
|---|---|
| `W` is bigger than `<w_m>`, so every Atkin–Lehner involution is already trivial on the quotient and none of them can be the hyperelliptic one — Molina's standing hypothesis in Section 10.2 | 23 |
| no divisor `n` of `D` has `2g+2` fixed points on the quotient, so again the hyperelliptic involution is not Atkin–Lehner | 5 |
| the hyperelliptic involution is `w_n` but `D` is not one of `m, n, n'`, so the two points at infinity are not exchanged by complex conjugation and Molina's choice of them in Sections 3 and 10.2 is unavailable — this is what rules out `X_0^57(1)` and `X_0^93(1)`, whose hyperelliptic involutions are `w_19` and `w_31` | 2 |
| inside the method but a factor has degree 5 or more, beyond `IndexFormEquation` | 8 |
| **inside and reachable** | **5** |

Of those 5, four are done: `85/<w_5>`, `85/<w_85>`, `115/<w_115>`, `161/<w_161>`, all of genus 2.
The fifth, `115/<w_5>` of genus 4, is the ceiling and is **not** solved — see below.

**The ceiling falls at genus 4.** All four reachable cases of genus 2 resolve; at genus 4 the
factor degrees are already `[1,3,3,3]` (D = 115) or `[5,5]` (D = 143 and 237), and from genus 6
upward every one of them has a factor of degree 5, 6, 7 or 8. This is exactly what Molina predicts
at the start of Section 10: the degrees grow with the genus until the index form equation is out of
reach. Note that it is the degree of one factor, not the genus itself, that decides: the degree-6
factor at `237/<w_79>` puts a *genus 2* curve out of reach.

## What is implemented, and what is not

Implemented, following the paper:

* **Step 1, the decomposition.** Which quadratic orders carry the Weierstrass points
  (Prop. 3.1 and Section 10.2) and the degree of each factor (Prop. 10.1, which is
  Gonzalez–Rotger Thm 5.12 — the repository's `DegreeOfFieldOfDefinitionOfCMPoint`).
* **Step 1, the discriminant of R.** At `p | D` the valuation of `disc R` is the sum of the
  thicknesses over the singular points of the special fibre of the quotient (Thm 2.3 summed over
  the colliding pairs, Lemma 4.1 for the thickness, Section 10.1 for the effect of `w_m`). The
  fibre and the Atkin–Lehner action come from `DualGraphData`. At 2 the valuation is `4(g+1)`,
  since `X_0^D` has good reduction at 2 for odd `D` (equation 7.30). Nowhere else.
* **Step 2, the leading coefficient.** The norm equation (6.26) and the intersection index (6.25),
  summed over the two involutions by the projection formula (10.32). The sign and square class come
  from the field of definition of `pi(P_inf)`.
* **Step 3, the fields.** `FieldOfDefinitionOfCMPoint`, put into an optimised presentation.
* **Step 4.** `IndexFormEquation`, then the resultants to fix the relative translations.

**Not implemented: the bimodule computation of Section 5.** Theorem 5.4 says exactly which singular
point each Heegner point specialises to, and hence how `nu_p(disc R)` splits between the individual
`disc p_{R_i}` and the resultants `Res(p_{R_i}, p_{R_j})`. In place of it the implementation
*enumerates* the possible splittings. The enumeration is small because three things constrain it
hard: each `disc p_{R_i}` is `J_i^2 disc(M_i)` for an integer `J_i`, so its valuation at `p` is
congruent to `v_p(disc M_i)` mod 2 and at least it; the product of all the discriminants and the
squares of all the resultants is `disc R`, which is known exactly; and a factor of degree `n` can
account for at most `floor(n/2)` singular points, so its discriminant is bounded by the sum of the
`floor(n/2)` largest thicknesses. What the enumeration costs is a few hundred to a few thousand
index form equations per curve instead of two or three; what it does not cost is correctness,
because every candidate is checked against the trace formula at the end.

**Not implemented: the blow-down.** Molina obtains the Weierstrass model `W` from Morita's model by
blowing down exceptional components, which changes the thicknesses by (4.15). The implementation
uses the fibre as `DualGraphData` gives it. That is a real gap, and it is *not* harmless in general
— it happens to be harmless on everything computed here because the blow-downs that occur merge
singular points along a single chain, and (4.15) then leaves the *sum* of the thicknesses, which is
all that `disc R` needs, unchanged. `diag.m` reports where a blow-down has occurred (the number of
singular points of the fibre is not half the number of Weierstrass points with singular
specialisation): it has, at 8 of the 24 prime-and-curve pairs computed here, and the answers were
right at all of them. **If a case is ever found where the sum is not preserved, this is the first
thing to suspect.**

**`115/<w_5>` is the one reachable case that did not come out**, with either choice of the field at
infinity (`a_R = -3` and `a_R = -8`), with the 2-adic part of each `J_i` allowed up to 2^6 and
indices up to 3 000 000. It is also one of the cases where `diag.m` says a blow-down has occurred,
at `p = 5`, and where the sum might therefore be wrong; that is the first thing to check. The other
possibilities are that the correct field at infinity is a third one, or that the 2-adic bound is
still too small.

## Files

| file | what it does |
|---|---|
| `molina.m` | the library: Steps 1 to 4 as functions. Load after `AttachSpec`. |
| `driver.m` | one curve end to end, with verification. `magma -b D:=39 m:=1 g:=3 logfile:=out.txt driver.m < /dev/null` |
| `equations.m` | the twelve equations obtained |
| `verify.m` | genus and point counts against the trace formula, for everything in `equations.m` |
| `compare.m` | isomorphism over Q against Molina's Table 2, transcribed |
| `census.m` | the 43 missing models, from a merged copy of every branch's `data/models` |
| `scope.m` | which of the 43 the method can reach, and why not for the rest |
| `diag.m` | where a blow-down has occurred |
| `probe2.m`, `probe3.m` | the two measurements the implementation was built on: that the thickness sum reproduces `nu_p(disc R)` on all ten published equations, and that `IndexFormEquation` finds the published generator at every index that occurs |

**These scripts need `DualGraph.m`, which is on `main` and not yet on this branch** — the special
fibre and the Atkin–Lehner action on it come from `DualGraphData`. Merge `main` down before running
anything here; without it every script stops at `Identifier 'DualGraphData' has not been declared`.
Everything was developed and measured against the tree at `main`.

Everything is run from the repository root, since `load` and `AttachSpec` resolve against the
working directory, and everything takes a `logfile:=` because Magma buffers stdout to a file and a
long run would otherwise look dead:

    magma -b D:=85 m:=85 g:=2 logfile:=out.txt vvdata/weyl-campaign/molina/driver.m < /dev/null
    magma -b logfile:=out.txt vvdata/weyl-campaign/molina/verify.m < /dev/null
    magma -b logfile:=out.txt vvdata/weyl-campaign/molina/scope.m < /dev/null

`driver.m` takes `maxj2` (how large the 2-adic part of each `J_i` may be, default 2), `maxidx` (skip
an index form equation above this index), `dinf` (force a particular field at infinity) and
`onlyfirst` (default 1: use only the smallest `|a_R|`). `census.m` takes `modeldir:=`, a directory
holding a copy of every branch's `data/models`; `compare.m` takes `res:=`, a file holding a Magma
list of `<D, m, [polynomials]>`.

Timings on one core, for the ten published cases: seconds to a couple of minutes each, except
`D = 57` with `w_3` (20 minutes) and `D = 35` with `w_5` (4 minutes), where a factor of degree 1 or
2 leaves its discriminant nearly unconstrained and the enumeration is correspondingly wide.

## Two traps

**The field of definition comes back with enormous coefficients.** `FieldOfDefinitionOfCMPoint`
returns, for instance, `x^4 - 4x^3 + 28831568902x^2 - 57663137796x - 440334637816243891199` for a
field that `OptimizedRepresentation` rewrites as `x^4 + x^3 - 2x^2 + x + 1`. `IndexFormEquation` on
the first form does not return in ten minutes; on the second it takes 0.09 seconds.

**`IndexFormEquation` returns its generators only up to translation by an integer**, and
`p(x) -> p(x+r)` moves the root of the monic `q(x) = a^{n-1} p(x/a)` by `a*r`, not by `r`. So the
translates have to be put back by hand, one residue class modulo `a` at a time; only one of them
gives a `p` with integer coefficients. Without this the search returns nothing at every leading
coefficient other than ±1, which looks exactly like the method failing.
