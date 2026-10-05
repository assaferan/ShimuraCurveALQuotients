# The m = 0 multiplier without a transcendental step (2026-10-05)

`exactm0.m` computes the m = 0 multipliers (1/2) c_eta(0) of Schofer's formula, per support class
of the isotropic coset eta, as EXACT rational numbers: every intermediate quantity is a root of
unity, a rational, or the positive square root of a rational, carried in the group ring
Q[Z/n] (polynomials modulo x^n - 1, n = 24 M) and reduced modulo the n-th cyclotomic polynomial at
the end, where the result must be a constant.

## The two closed forms

The per-coset product rho*(w^-1) e_0 [eta] * a_0(f | w) is constant on the class g = gcd(c, M)
(thetag-derivation.md), so one representative gamma = [a b; c d], c > 0, per class suffices, with
the canonical lift (gamma, sqrt(c tau + d)) on BOTH factors:

* **Weil representation.**  rho*(gamma^-1) e_0 [eta] = rho_B(gamma)_{0,eta}
  = e(1/8) c^{-3/2} |L^v/L|^{-1/2} e(d Q(eta)/c) G(c, a; eta), with
  G(c, a; eta) = sum_{nu in L/cL} e((a Q(nu) + (eta, nu))/c), a 3-dimensional Gauss sum that splits
  over the prime powers of c (the quadratic term picks up the cofactor c/p^k, the linear term does
  not).  Derived from the theta transformation formula (gamma tau = a/c - 1/(c(c tau + d)), Poisson
  summation on the inversion); the factor e(d Q(eta)/c) carries d, the Gauss sum carries a, and the
  swapped version is wrong by O(1) (rhotest2.m).
* **Slash constant.**  With [d 0; 0 1] gamma = g_d [a_d b_d; 0 e_d] (g_d in SL_2(Z) with positive
  lower-left entry), Apostol Thm 3.4 gives
  a_0(prod_d eta(d tau)^{r_d} | gamma) = e(-1/8) prod_d [eps(g_d) e(b_d/(24 e_d)) e_d^{-1/2}]^{r_d}
  times the constant term of prod_d [prod_n (1 - e(n b_d/e_d) q^{n a_d/e_d})]^{r_d} q^{sum r_d a_d/(24 e_d)},
  eps the Dedekind-eta multiplier (a 24th root of unity from a Dedekind sum).  For the polytope
  forms every monomial is holomorphic at the cusps other than oo and 0, so at the middle classes the
  constant term is 1 or 0 and no series is needed; at the cusp 0 the representative S gives a
  series with RATIONAL coefficients (eta(tau/d)); the general q-series path over Q(zeta) exists
  for monomials with poles elsewhere (the LLL probe forms of compmult5.m).

The e(1/8) and e(-1/8) cancel.  Class 1 (cusp 0) uses S itself; class M (the identity) enters only
c_0(0), which must vanish.

## Validation

* `rhotest.m`, `rhotest2.m`: the rho row from the FFT (VVRhoInvE0FFT), from the square of the c = 1
  matrix, and from the Gauss-sum formula agree to 1e-58 for c = 1, 2, 3, 4, 5, 10, including the
  representatives [3 1; 5 2], [1 2; 4 9], [7 2; 10 3] with a != d; the slash constant from the closed
  form agrees with the numerically pinned one on 140/140 (monomial, point) pairs.
* `exactm0_15_2.log`: all nine multipliers of X_0^15(2) equal the measured ground truth
  (tests/M0MultiplierExact.m), c_0(0) = 0, every class agrees between two representatives, 3 s;
  M0MultipliersBySupport agrees (9.4 s for both).
* Further logs: `exactm0_21_2.log`, `exactm0_10_3.log` (with the numeric routine), `exactm0_6_35_lll.log`
  (the four compmult5 forms: expected {5}, {7}, {5,7} = 1/2, 1/2, 1/4; -23/4, -39/4, 0; 78, 64, 0;
  -1/4, 9/8, 0).

## Traps met

* A RATIO test of two vectors that mostly vanish is noise: compare by difference with the lift sign
  fixed on a nonzero component (the first RHOCHECK "failed" on every class except c = 1 for this
  reason while the formula was right).
* `gt` is a reserved word in Magma.
* The CRT split of the Gauss sum: only the quadratic term gets the cofactor.  With the linear term
  multiplied by a too, every representative with a = 1 is right and the two-representative check
  catches it at the first a != d (class 5 of M = 60).
