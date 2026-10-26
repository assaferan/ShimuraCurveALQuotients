# The local Whittaker factor at a split level prime

`level-prime-kappa.tex` — a self-contained account of the level-prime obstruction found while
computing CM values of Borcherds forms on `X_0^D(N)`.

Build with `pdflatex level-prime-kappa.tex` (twice, for cross-references). No external `.bib`;
the bibliography is inline.

What it records, and the status of each claim:

* **Theorem 3.1** — closed form for the local Whittaker factor of the N-scaled binary lattice at a
  split level prime, all valuations. Verified in exact arithmetic (180 checks, N = 2, 3, 5,
  m ≤ 60); pinned by `tests/EisensteinLocalFactors.m`.
* **Corollary 3.2** — the empirical "level support rule" is exactly the vanishing of that factor.
* **Proposition 4.1** — re-running Kudla–Rapoport–Yang's case analysis with the level section adds
  no contributing configuration.
* **Section 5** — the correction the CM values require, measured exactly, and the table any
  corrected formula must reproduce.
* **Section 6** — eight eliminated explanations, with the reason each fails.
* **Section 7** — the two possibilities we cannot presently distinguish.

The computations live in `VectorValuedForm.m` and `EisensteinLocalFactors.m`; `vvdata/` banks the
runs the numbers come from.

# X_0^21(58)* is hyperelliptic

`hyperelliptic-21-58.tex` — the proof, with no equation for the curve, and the lemma it needs. Four
CM divisors of degree two are pairwise torsion-related against a torsion bound of 2, so two coincide
and the genus-three curve carries a `g^1_2`. The torsion relations are certified by the dual test
against weight-3/2 forms; at a level with `2 || N` the plain test is blind to half of the relevant
forms, and the note proves (Section 3) that the Heegner components of `X_0^D(N)*` are
`(U(4) + 2)` applied to Kohnen plus-space forms attached to `X_0^D(N/2)*`: a fibre count along the
degeneracy map, the injectivity of the correspondence on a 2-old Jacobian, and Borcherds' modularity
of the vector-valued Heegner series restricted to the 2-part of the discriminant form. Tests
`tests/Hyperelliptic21_58.m` (minutes) and `tests/_offline/Hyperelliptic21_58.m` (half an hour);
solver in `tools/proof-21-58/`.
