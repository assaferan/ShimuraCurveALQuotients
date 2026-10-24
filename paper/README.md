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

# Twenty genus-three quotients from the heights of their branch points

`heegner-height-models.tex` — the argument behind the twenty models of `data/models` added on this
branch (PR #128), written as a proof: the branch locus of `X -> E` from Ogg's fixed-point counts
(Section 2); the explicit Gross–Zagier formula of Cai–Shu–Tian in the normalisation the code uses,
the admissibility rule, the rank-one shift by the Hodge class, and the oldform multiplier (Section
3); why the candidate branch divisors form a finite, complete list, and the two places searched over
a window (Section 4); the classification of the double covers along a given divisor and the
proposition that the point counts select the curve (Section 5); what was predicted rather than
input, and what the theorem rests on (Section 6). Table 1 lists all twenty with `W` by generators,
`E`, the branch places and the model; the two candidates and the two window-conditional cases are
marked. Code: `HeegnerHeights.m`, `ExactDoubleCover.m`, `tools/heegner_model.m`; the per-curve run
record is the PROVENANCE section of the same name.
