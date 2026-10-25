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

# X_0^6(341)/<w_3, w_11, w_62> is hyperelliptic

`hyperelliptic-6-341.tex` — the proof: the Borcherds certificate of PR #142 (the branch pair is
canonical up to the star's torsion), the star's Jacobian identified with that of the curve C of
PR #143 up to isogeny (CHIMP endomorphisms and conductor, Ribet, the orbit at level 2046), the
isogeny class shown to contain exactly two principally polarised surfaces (a chain of length two in
the Bruhat–Tits tree, since the star has a rational 2-torsion point), the second one computed as
the Jacobian of an explicit curve C_A, and C_A excluded because no double cover over Q can be
branched along any class the certificate allows. Scripts and CHIMP logs in `tools/proof-6_341/`.
