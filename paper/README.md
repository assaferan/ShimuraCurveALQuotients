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
  evaluation of a pullback as the level-`M` sum at the moved point times a phase. Stated with a
  proof, the sign convention at the ramified primes resting on the Jacquet–Langlands
  correspondence; the dimension counts and the measured eigenvalues of each column are
  consistency checks. Not to be called proved until both maintainers have read it.
* **Section 3** — the sample points and the cost law (linear in `DN`), measured.
* **Section 4** — Noether's count from numerical values, with the reading rule and what a
  numerical kernel dimension does and does not prove.
* **Section 5** — the model step through the elliptic quotient: the marked invariants, the
  certainty criterion, the twist constants from point counts, and where the fit is ill
  conditioned.
* **Section 6** — the four checks and their status: numerical plus point counts; agreement with
  the exact model of `X_0^1326(1)/<w_17, w_78>`, and at `X_0^6(185)` agreement of the computed
  invariants with an exact model to twelve digits although the recognition step refuses them.
* **Section 7** — the two traps: precision must rise with the cutoff; the Eichler order is not
  deterministic across processes.

The code is `NelsonShimizu.m`, tested by `tests/NelsonShimizu.m`; `data/nelson/` banks the stored
values the model test reads.
