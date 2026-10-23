# A genus-two curve from a modular abelian surface, by its periods

The Jacobian of a genus-two Atkin–Lehner star is a modular abelian surface, so its periods come
from a classical newform — where precision is cheap. That is the point of this route. On the
Nelson–Shimizu side every digit costs an evaluation of Theorem 3.1, minutes apiece, and the
recognition of invariants stalls around twelve digits; here 198 real digits take under a second at
small level. Both routes aim at the same object, the exact model of a genus-two star.

The chain is: periods → principal polarization → Riemann matrix → Rosenhain invariants → Igusa
invariants → recognition over **Q** → Mestre's construction → minimal model → the quadratic twist,
determined from point counts against the newform.

`periodcurve.m` runs it at a level, taking the first two-dimensional factor of `JZero(N)`.
`periodstar.m` is the same with the factor selected by Atkin–Lehner signs: the star's differentials
are the newforms whose classical sign is −1 at every prime dividing D, which after the
Jacquet–Langlands flip is the quaternionic +1 defining the star.

## The step that is the actual content: the polarization

Magma's `Periods` and its `IntersectionPairing` are **not in the same homology basis** —
symplectifying the pairing and applying that transformation to the periods leaves the Riemann
matrix asymmetric by 0.94 — so they cannot simply be composed. The polarization is therefore
recovered *from* the periods. For a period matrix Π (2 × 4, columns the lattice generators) and an
integral alternating Q, Riemann's first relation Π Q Πᵗ = 0 is one complex equation in the six
unknowns q_ij; its integral solutions are the Néron–Severi group, of rank 2 for a surface with real
multiplication. Among them the principal polarizations are those of Pfaffian ±1, cut down by
Riemann's second relation (positivity). This is the computational face of the classification of
principal polarizations on GL₂-type surfaces (González, Guàrdia, Rotger, *Abelian surfaces of
GL₂-type as Jacobians of curves*, Acta Arith. **116** (2005), 263–287).

## Controls: two pass, and the third does not

| level | factor | curve recovered | outcome |
|---|---|---|---|
| 23 | `JZero(23)` entire | `y² + (x³+x+1)y = −2x⁵ − 3x² + 2x − 2` | equals `SmallModularCurve(23)` |
| 67 | a proper quotient of a 5-dimensional J₀ | `y² + (x³+x²+1)y = −x³ − x` | counts agree at 41 primes |
| **35** | the dimension-2 factor of `JZero(35)` | `y² + (x³+x)y = −x⁵−x⁴+3x³−3x²−5x` | **⚠ NOT the Guo–Yang model of X_0^35(1)/w_5** |

Level 23 is the easy case: the polarization is principal for free, because J_0(23) *is* a
genus-two Jacobian. Level 67 is a proper quotient, which is the situation at every star, and the
route returns a curve its polarizations all agree on. **Level 35 is the one that matters most**,
because it is the only control with an independent model: there the route returns a curve that is
not the right one, and the point counts cannot tell. See entry 8 below and `models.md`.

## ⚠ Seven ways this computation returns a confident wrong answer

Every one of these produced plausible output rather than an error, and all but one were caught only
by the known-answer control. They are recorded because the same shapes will recur.

1. **`Periods(A, n)`: `n` is the NUMBER OF TERMS, not digits**, and the accuracy it delivers is not
   what the complex field claims. At level 23 the real accuracy runs at about n/7 digits and then
   saturates against the ambient field, which defaults to 30 regardless of what is requested.
   `SetDefaultRealFieldPrecision(p)` lifts the ceiling and then gives p − 2 real digits.
2. **The terms needed scale with the LEVEL.** The expansion converges with |q| ~ exp(−2π/N), so the
   digits obtained run roughly as terms/(0.37·N). At N = 2002, 1600 terms give about *two* digits —
   which is why the first run there found no polarization at all.
3. **The LLL multiplier must be MODERATE, not maximal.** Set to the full precision it skews the
   lattice so badly that LLL returns only its long vectors (heads of 10⁴⁷ and more) and the genuine
   short relations never surface. Sweep it and keep what recurs.
4. **A relation is believable because it is SHORT, not because its residual is small.** LLL's long
   vectors achieve the same *relative* residual as genuine ones; comparing residual to coefficient
   size passes both.
5. **`[ f(i,j) : i in I, j in J ]` varies i FASTEST in Magma.** Used to build a 2 × 4 conjugate
   matrix row-major, it silently scrambles the entries — which suppressed every principal
   polarization while reporting none exists.
6. **Mestre's output is far from minimal** — 128-digit coefficients at level 23 — and reducing such
   a model mod small p gives bad or singular reductions, so point counts taken on it are
   meaningless. `ReducedMinimalWeierstrassModel` first.
7. **`Trace(HeckeOperator(A, p))` is 2·Tr_{K/Q}(a_p)**, not the trace of a_p: T_p acts on the
   four-dimensional homology as a_p on K². Verified against X_0(23), where the ratio to
   p + 1 − #C(F_p) is exactly 2 at every prime.

Two smaller ones: the real-multiplication field must be read at a prime **not dividing** the level
(at 39 the characteristic polynomial of T_3 is (x−1)⁴ because 3 | 39, so the test is vacuous), and
primes with a_p = 0 carry no information about the twist's character and must be excluded from the
sign-pattern fit.

## ⚠⚠ 8. THE POINT COUNTS DO NOT IDENTIFY THE CURVE

This is the most important entry and it invalidates the obvious reading of the section below. All
the principal polarizations on the surface give curves with isomorphic *unpolarized* Jacobians, so
they share an L-function and agree in #C(F_p) at **every** prime. The counts therefore confirm the
isogeny class and cannot select among the curves. Demonstrated at level 35 against the independent
Guo–Yang model of X_0^35(1)/w_5: the period route returns a curve with different G2 invariants —
so not even geometrically isomorphic — whose point counts agree at all eight primes tested. See
`models.md`. Until a selector exists, a curve produced here is a candidate, not a model.

## The certification, and what it does and does not establish

Mestre's construction gives the curve only up to quadratic twist, and the right twist is not
generally small — at level 23 it is −16108971167423174820, so scanning |d| ≤ 2000 finds nothing.
The twist is therefore *determined* rather than guessed: at each good prime the count is
p + 1 ∓ a_p, the sign pattern is the quadratic character of the twisting discriminant, and that
discriminant can only ramify where the model has bad reduction, so it is found among signed
squarefree products of those primes. A wrong *twist* has no single character fitting all the signs, so this does pin the twist. It does
not pin the curve: see the warning above.
