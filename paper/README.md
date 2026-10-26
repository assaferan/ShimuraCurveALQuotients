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

# X_0^77(10)* is hyperelliptic

`hyperelliptic-77-10.tex` — the proof, with no equation for the curve: five CM divisors of degree
two (four pairs of rational points, two of them outside the class-number-four screen, and one
quadratic point) are pairwise torsion-related against a torsion bound of 2, so two coincide and the
genus-three curve carries a `g^1_2`. The relations come from the dual test for `2 || N` proved in
`hyperelliptic-21-58.tex`, run as a search over all CM divisors of degree at most two. Tests
`tests/Hyperelliptic77_10.m` (a minute) and `tests/_offline/Hyperelliptic77_10.m` (twenty-five
minutes); solver and scan in `tools/proof-77-10/`.
