
## 2026-10-05 (late evening) — exactm0/: the m = 0 multiplier by exact algebra

exactm0/exactm0.m computes (1/2) c_eta(0) with no transcendental step: one coset per class g = gcd(c, M),
rho row = e(1/8) c^{-3/2} |L^v/L|^{-1/2} e(d Q(eta)/c) sum_{nu in L/cL} e((a Q(nu) + (eta,nu))/c) (theta
transformation), slash constant = Dedekind-eta multiplier (Apostol 3.4), everything in Q[x]/(x^n - 1) with
n = 24 M and reduced mod the cyclotomic polynomial. 15_2 = ground truth 9/9 in 3 s; 21_2, 10_3, 22_3 =
M0MultipliersBySupport 9/9; the four compmult5 monomials at 6_35 = compmult5's twelve values; the production
forms of 6_35/10_21/14_15 (17 each) = the genmodels logs as multisets (cmpprod.py); 34_11 in 148 s. Now the
production route (branch exact-m0, M0MultipliersAlgebraic). Traps: ratio tests of mostly-zero vectors; the
CRT split of the Gauss sum (only the quadratic term takes the cofactor); `gt`/`LLL` are reserved. README there.
