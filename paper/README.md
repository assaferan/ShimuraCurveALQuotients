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

# The CM value formula, stated once for every coset and index

`schofer-formula.tex` — the formula `SchoferFormula` evaluates (Theorem 2.2 with Theorem 2.3 and
Corollary 2.4), with a proof that treats composite squarefree level and non-fundamental
discriminant on the same footing as the prime, fundamental case: the level, the conductor and the
ramification of the primes of DN enter only through the local lattices, and the statement about
those is one formula (Theorem 2.3(c)); the case list is Part II, carried over from
`kappa0-proof-standalone.tex`, which this note is meant to supersede. Draft; not yet agreed.
