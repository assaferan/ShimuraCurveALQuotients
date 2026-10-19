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

# Evaluating forms on Shimura curves by the Shimizu lift

`nelson-method.tex` — how the differentials of an Atkin–Lehner quotient are evaluated at points
of the upper half-plane by Nelson's formula for the Shimizu lift, the extension of that formula
from the newforms of level exactly `DN` to the oldform part of the space, and the two uses of
the values: Noether's count of canonical quadrics (hyperellipticity, numerically) and the
hyperelliptic model of a genus-3 quotient with an Atkin–Lehner involution of genus-1 quotient.

Build with `pdflatex nelson-method.tex` (twice). Bibliography inline.

What it records, and the status of each claim:

* **Section 1** — Nelson's Theorem 3.1 and its hypothesis (newforms of level exactly `DN`).
* **Section 2** — the oldform extension: the oldclass of a newform of level `DM` as the span of
  the pullbacks `g|w_S`, the action of the Atkin–Lehner group on it (Proposition 2.1), and the
  evaluation of a pullback as the level-`M` sum at the moved point times a phase. Checked by the
  dimension counts and by direct measurement of every column's eigenvalues; a proposition with a
  proof, not a numerical claim.
* **Section 3** — the sample points and the cost law (linear in `DN`), measured.
* **Section 4** — Noether's count from numerical values, with the reading rule and what a
  numerical kernel dimension does and does not prove.
* **Section 5** — the model step through the elliptic quotient: the marked invariants, the
  certainty criterion, the twist constants from point counts, and where the fit is ill
  conditioned.
* **Section 6** — the four checks and their status: numerical plus point counts; one agreement
  with an exact model (`X_0^1326(1)/<w_17, w_78>`).
* **Section 7** — the two traps: precision must rise with the cutoff; the Eichler order is not
  deterministic across processes.

The code is `NelsonShimizu.m`, tested by `tests/NelsonShimizu.m`; `data/nelson/` banks the stored
values the model test reads.
