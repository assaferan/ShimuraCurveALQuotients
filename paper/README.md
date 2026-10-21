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

# Hyperellipticity from several Borcherds divisors

`torsion-counting-hyperelliptic.tex` — why a search that returns several degree-two Borcherds
divisors proves hyperellipticity outright, where one divisor proves it only under Guo-Yang's
Assumption 40.

Build with `pdflatex torsion-counting-hyperelliptic.tex` (twice). No external `.bib`.

* **Lemma 1.1** — a Borcherds form with integral principal part has a divisor class of finite order
  in J(Q). Borcherds' theorem, with the correction that restores the finite-order condition.
* **Theorem 1.2** — more effective degree-two divisors, pairwise differing by such divisors, than
  the rational torsion has elements, forces two of them to be linearly equivalent, hence a function
  of degree two. The assumption becomes a count.
* **Section 3** — why both sides of the count are computable with no model: the torsion bound from
  the trace formula, and the forms from their recorded combinations of eta quotients.
* **Example 4.1** — X_0^46(31)*, genus three with no Atkin-Lehner involution left: four effective
  divisors against a torsion bound of 2, so it is hyperelliptic unconditionally.

Checked by `tests/Hyperelliptic46_31.m`, which rebuilds the six forms from
`data/borcherds/forms_46_31.m` and recomputes the torsion bound.
