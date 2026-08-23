# The per-class functional sitting (2026-08-23 evening): a0_g is a cusp constant term, and the assembly is a weight-3/2 residue pairing

Scripts: classfunc.py (Ramanujan-modulus scan), classfunc2.py (exact per-class
rref over the pp channels). Data: the six cusp3_*_full dumps + pp3/pp2 dumps.

## Facts measured tonight (exact, panel-level)

1. **Inter-class relations of the S_g row vectors** (S_g = class sum of
   contrib = (M/g) rho_g a0_g, snapped to rationals from the 10-digit dumps):
   * 15_2: rank 3/4 on the 9-form panel: **3 S_1 = 2 (S_3 + S_15)**.
   * 35_2: rank 2/4 on its 8-form panel: **S_5 = -S_1** and
     **S_35 = 2 S_1 + S_7**.
   * 10_7, 21_2, 33_2, 6_11: live classes independent on their panels.
   (Panel-span statements, not yet intrinsic.)

2. **The per-class functionals are NOT principal-part functionals.** On 21_2
   every live class (1, 3, 7, 21) is exactly INCONSISTENT as a linear
   functional of the dumped oo- and 0-side pp channels, and adding a
   square-class column sum_k c(-k^2) does not restore consistency. (On the
   other five bases the systems are consistent but underdetermined -- no
   contradiction, their panels are smaller than their channel lists.)

3. **The class g=1 sums are the 0-cusp constant terms.** S_1 on 21_2 =
   (6, -4886, 10, -4886, 6, 12, -2, 2, 8) over forms (-2,-1,9..15): the
   -4886 is the long-recorded 21_2 0-cusp constant (p2-local-density memory),
   appearing on exactly the forms (-1, 10) where that handle lives. Together
   with the exact S_1(f) = Coefficient(f0, 0) on all nine 15_2 forms (banked
   3391b33), the identification is: **a0_g(f) is literally the constant term
   of f at the class-g cusp of Gamma_0(M)** (as the assembly derivation
   states), so no per-class pp law exists to find. The Ramanujan-modulus
   ansatz per class (classfunc.py) was ill-posed for this reason -- recorded
   so nobody refits it.

## The route this forces (the reading, to derive next)

The assembly c_{eta*}(0) = sum_g (M/g) rho_g a0_g(f) is the constant-term
(residue) pairing of f against a FIXED weight-3/2 object:

* Let E be any weight-3/2 form on Gamma_0(M) with multiplier conjugate to
  f's eta multiplier (so f*E is a weight-2 form with trivial multiplier).
  Sum over all cusps of the constant terms of f*E vanishes (residue theorem).
  Expanding: sum_m c_f(-m) a_E(m) [one block per cusp where f has principal
  part -- oo AND 0: this is exactly the JOINT (oo,0)-channel structure the
  cusp5/cusp6 sitting discovered] + sum_cusps a0_cusp(f) * a0_cusp(E) = 0.
* If E's cusp-constant vector matches (M/g) rho_g -- now fully explicit from
  the rho-entry theorem (support, g sqrt2/M and g/M magnitudes, E0 phase,
  mu(g/4)) -- then the m=0 multiplier functional IS m |-> a_E(m):
  **W(m) = Fourier coefficients of the weight-3/2 Eisenstein series with
  prescribed cusp constants rho**. DERIVE the ledger, don't fit it.
* Well-definedness: E is unique up to cusp forms, and cusp forms pair to
  zero against every principal part (Bruinier-Funke) -- reproducing the
  "weights defined modulo the relation ideal" phenomenon EXACTLY, and the
  per-base rationality (Eisenstein coefficients at level M are rational).
* Expected dividends: the square-class channel = the Zagier E_{3/2}
  component of E_rho (w_sq = its coefficient; conjecture 1/(N-1) becomes an
  Eisenstein computation); the psi(k) = lambda(k)(k/N) twist = the character
  of the relevant Eisenstein component; the N|m support = rho's support law
  pushed through the cusp-constant matching; base-dependence of the kappas =
  dependence of the Eisenstein basis on M -- no universal symbol law needed,
  explaining the three refuted conjecture families.

## Next computational step (a sitting of its own)

Build the weight-3/2 Eisenstein space at level M with the conjugate eta
multiplier (per base), impose the rho cusp constants, and read off W(m).
First targets with full acceptance data: 15_2 (level weights -1, +1, 0 and
w_sq = 0), 22_3 (w_sq = 1/2, kappa_11 = 3/2), 21_2 (w_sq = 1), then the
forced singles (10_3, 6_5, 10_7, 6_13, ...). The multiplier bookkeeping
(eta character of each panel's forms) is in the cusp1/cusp4 machinery.
