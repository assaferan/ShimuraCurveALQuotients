# The level-prime local Whittaker at GENERAL m, both cosets (2026-09-10)

Computed, not fitted. `wlocal.py` in this scratchpad. Setup is the paper's own recipe:
  alpha_k = p^{-k(n-1)} #{x in (mu+L)/p^k L : Q(x) = m mod p^k},  n=2
  W_{m,p}(s,phi) = (1-X) G(X),  G = sum alpha_k X^k,  X = p^{-s}
At p = N the lattice is the level-N hyperbolic plane, L = Ze+Zf, Q(ue+vf) = N u v; then
L^v/L = (Z/N)^2 with Q(a/N,b/N) = ab/N, so isotropic <=> ab = 0 mod N: exactly 2N-1 cosets,
which is MEASURED FACT 1 of [[isotropic-cosets-live-at-level]] (2N-1 over 12 bases).

## VALIDATION FIRST (reproduce a known value before trusting a new one)
mu = 0, m = 0 gives coefficients 1, N-1, N-1, ... at N = 2,3,5
  = (1 + (N-2)X)/(1-X), which is the paper's recorded W_{0,N}(s,phi_0) (Lem 11.1 / Prop 11.2).
As s->0, X->1: numerator -> N-1, denominator 1-X ~ s log N, so W ~ (N-1)/(s log N) --
the exact POLE recorded in [[m0-kappa0-solved]]. Both the value and the pole are reproduced.

## RESULT 1 -- ZERO coset, general m
    W_{m,N}(X, phi_0) = 1 + (N-1)(X + ... + X^v) - X^{v+1},   v = ord_N(m)
Verified N = 2,3,5 and m = 0..6. Specialises to the validated m=0 case as v -> infinity.
    ==> W_{m,N}(s=0, phi_0) = (N-1) * ord_N(m)      for m != 0
    ==> IT VANISHES EXACTLY WHEN N does not divide m.

## RESULT 2 -- NONZERO isotropic coset, general m
    W_{m,N}(s, phi_mu) = 1   identically, for EVERY m and every nonzero isotropic mu.
By hand: for mu = (a/N, 0) with a a unit, Q(x) = t*(a + N s), so for each s there is exactly
one t mod N^k; alpha_k = 1 for all k, G = 1/(1-X), W = 1. Verified N=2,3,5, m=0..3, all mu.
Extends the recorded m=0 fact ("identically 1 at every nonzero isotropic coset") to all m.

## WHAT THIS BUYS
* The level prime is now CLOSED for general m at every isotropic coset -- the piece the open
  theorem was phrased around ("general m at a nonzero isotropic coset") is the EASY factor.
* Since a nonzero isotropic mu lives ENTIRELY at N and contributes 1 there, b_mu(m) carries no
  level-prime data at all: all 2N-2 nonzero isotropic cosets must give the SAME coefficient.
  (Consistent with m=0, where all 2N-2 gave kappa = -log N/(N-1).)
* ==> THE DIFFICULTY IS NOT AT N. It is at p | D (anisotropic/ramified) and at infinity --
  exactly where [[b-eisenstein-coefficients-solved]] localised the defect ("P VANISHES at
  supported indices ... because the local factor at p = 3, a prime dividing D, is 0").

## ⚠ TENSION TO RESOLVE BEFORE QUOTING ANY OF THIS
A_r = -b^{eta*}_0(r)/4 has SUBSCRIPT 0 -- the ZERO-coset component. If that is right, RESULT 1
says the level factor vanishes iff N does not divide r, which would DERIVE the support rule
N | r. But [[support-rule-is-a-gauge]] records that rule as a GAUGE, "undecidable by
derivation" (Conjecture 8.2). Either (i) the gauge entry concerns a different support statement,
(ii) the subscript is eta not 0 in the case that matters, or (iii) one of the two is wrong.
DO NOT claim the support rule is derived until this is settled -- read the paper's Conj 8.2 and
the b^{eta*}_eta(r) indexing at [[b-eisenstein-coefficients-solved]] first.
Regression set for any claim: the 9 exact b at 15_2, 6_5, 10_3.

---

## ⚠⚠ RETRACTION (same day, before anyone builds on this): NEITHER RESULT IS NEW.

Both are already in `paper/level-prime-kappa.tex`, and I should have read it before computing:

* **Result 1 IS Theorem `thm:closed`**, verbatim:
  `W_{m,N}(X) = 1 + (N-1)(X + ... + X^j) - X^{j+1}`, `W_{m,N}(1) = (N-1) ord_N(m)`,
  with **`cor:support`**: it vanishes iff `N` does not divide `m`. Verified there for
  `N = 2,3,5` and all `1 <= m <= 60` on the `Lm` of 15_2, 6_5, 10_3 -- **180 checks**, against
  my 3 values of `N` and `m <= 6`.
* **Result 2 is in `sec:open`**, which carries the same `alpha_k`/`G(X)` recipe AND the counts:
  "at `mu != 0` ... `alpha_k = 1`; at `mu = 0` one has `xy = 0 mod p^{k-1}`, giving
  `alpha_k = (k-1)(p-1)+p`". That is exactly this computation.
* The paper also already records the consequence I drew ("the difficulty is at `p | D`"), via
  `prop:closedcoef` (coefficients supported by an EMBEDDING condition at the primes of `D`).
* And the gauge worry was answered there too: `rem:gauge` explains why the `N | m` support rule
  (a statement about a REPRESENTATIVE) and `prop:closedcoef`'s embedding-supported coefficients
  are compatible, and why no panel can separate them.

**What actually remains of this work:** an INDEPENDENT brute-force reproduction. The paper itself
asks for one -- "That test alone would establish only that the code is self-consistent, since the
closed form was read off the same implementation" -- though `sec:open` then supplies one, so this
is at best a second. Modest, and not what I claimed.

**HOW I GOT IT WRONG, since the pattern is the point.** I read the memory entries, saw "the next
theorem is general `m` at a nonzero isotropic coset", and inferred the level prime was open at
general `m`. It is not: `thm:closed` IS general `m` at the level prime. The memory sentence means
the *intersection* -- general `m` AND the `D`-part -- and I resolved its ambiguity in the direction
that made my computation look new. Every number I produced was right; I was wrong about **which
object was already known**, which is the failure mode `CLAUDE.md` opens with, in a fresh disguise:
validating arithmetic cannot detect redundancy. **Read the paper before deriving, not after.**
