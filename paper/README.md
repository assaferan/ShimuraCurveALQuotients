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

# X_0^6(379)* is hyperelliptic

`hyperelliptic-6-379.tex` — the proof, with no equation for the curve: a Borcherds form on the
star, given only by its principal part at infinity (a weight-1/2 form of level 4548 whose existence
is decided by linear algebra against the 377 weight-3/2 cusp forms of that level), has divisor a
positive multiple of `P_-123 + P_-163 - P_-52 - P_-168`; the Jacobian has no rational torsion, so two
distinct effective divisors of degree two are linearly equivalent and the genus-three curve carries
a `g^1_2`. The certificate is `data/dual/cert_6_379.m`, checked by `tests/DualCertificate6_379_divisor.m`
(12 s) and `tests/_offline/DualCertificate6_379.m` (an hour); the solver is `tools/proof-6-379/`.
