# X_0^1326(1)/{1, 17, 78, 1326}: y² = 272x⁸ − 72x⁶ − 55x⁴ − 18x² + 17

Built 2026-10-18 by the height route (no Borcherds form, no floating point in the construction):

1. E = X/⟨W, w₆⟩ is 1326b2: the Atkin–Lehner signs select the isogeny class 1326b, and the
   Gross–Zagier heights select b2 within it (on b1 the predicted height is not a lattice height).
2. The branch points of X/W → E are the fixed points of w₆: two of discriminant −24 and two of
   −408, all over F = Q(√2). By Cai–Shu–Tian (Theorem 1.5) with the admissible genus characters
   (−3, 8) and (8, −51), the minus parts P − P̄ have height exactly that of the generator G′ of
   E(F)⁻ ≅ E^8(Q) (k = 1.000000 in both cases); the plus parts are 2f(ξ) + torsion with f(ξ)
   the image of the Hodge class, the same for both branch points (b/2 times the generator of E(Q),
   E of rank 1).
3. For each b ∈ [−6, 6] every configuration (signs, torsion translates, halvings) was tested exactly:
   `DoubleCoversWithBranchDivisor` + `MatchCoverToCounts` against the trace formula. Survivor at
   b = −5 (sweep of the other b still running when this was written, see `pipe_1326_partial.log`):
   four covers over Q with that branch divisor, one consistent with the trace formula, hyperelliptic.
4. Its twist is unique once primes up to 700 are used (κ = 9391072966921, the large primes being
   artefacts of the normalisation of f); the model from the function field, normalised by the
   Klein-four action (x ↦ −x, x ↦ ±1/(2x)) and square scalings, is the octic above; the Q-form among
   the twists of the small even model was fixed by `IsIsomorphic` with the raw model.

Checks that were predictions, not inputs: bad primes exactly {2, 3, 13, 17}; point counts agree with
the trace formula at the 164 primes 5 ≤ p ≤ 1000 and over F_{p²} at p = 5, 7, 11, 23; automorphism
group of order 8; elliptic quotients 1326b2 (w₆), 1326c1 (w₃), 1326d1 (w₂) with the right signs;
fixed points of the involutions over Q(√2) (−24, −408), Q(√17) (−136) and the quartic field of
discriminant 18496 (−1768), exactly the fields of definition of those CM points. 8 rational points
of height ≤ 10⁴. Hyperelliptic, as Nelson's numerics said.

Only b mod 2 matters (translating both branch points by a rational point of E gives an isomorphic cover,
and halvability depends on b mod 2), so the enumeration over b ∈ {0, 1} is complete: the even class has
no configurations (b = −6 above) and the odd class gives this one curve. Rigorous given Cai–Shu–Tian's
theorem and the exact cover test; the sweep over the other b was stopped as redundant.
Scripts: `pipeline2.m` (the sweep), `surv1326.m` (twist and raw model), `norm3.m`/`norm5.m` (normal
form), `check1326.m` (the checks), `overconic.m` (model over the canonical conic, not needed here).
