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

# The character of a Borcherds form in positive genus, and reduction types

`character-and-reduction-types.tex` — what Guo–Yang's Assumption 40 amounts to and how to do
without it, and the Kodaira type of a genus-one quotient from the Čerednik–Drinfeld dual graph.
Build with `pdflatex character-and-reduction-types.tex` (twice). Inline bibliography.

* **Proposition 3.1** — given Lemma 7's parity condition, the character is trivial exactly when the
  divisor is principal (not merely torsion); its order is the order of the class in the Jacobian.
* **Corollary 3.2** — no rational torsion among CM classes, no assumption needed; so 142a1, 158a1,
  190a1 are unconditional.
* **Example 3.4** — the assumption appears to fail at X_0^14(17)*: one run of the genus-one star
  driver reported a form whose value at the rational CM point of discriminant −4 is 2·7^{3/2}·17^{1/2},
  which Remark 3.3 reads as a character of order 2. The driver and its log are not in this branch, so
  the example is an observation and not reproducible from the files beside the note.
* **Section 4** — e-th powers of the forms remove the assumption.
* **Proposition 5.1** and the table — the Kodaira index is the sum of the thicknesses of the nodes on
  the cycle, that is half the length of the metric cycle (`KodairaIndexFromDualGraph`,
  `tests/KodairaFromDualGraph.m`): 238a2 and 282b1. The method is González–Rotger's;
  `tests/KodairaGonzalezRotger.m` runs it on the 29 genus-one curves and quotients of their two
  tables, taking each isogeny class from the trace formula, and reproduces every identification.
