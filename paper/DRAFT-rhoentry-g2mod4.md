# DRAFT — the missing case of `thm:rhoentry` (ii): `g ≡ 2 (mod 4)`

**Status: a measured OBSERVATION plus an unresolved discrepancy with the proof's own `|A[c]|`.**
Not proposed text yet — see §3, which is why. Found 2026-09-05 while writing
`tests/RhoEntryClosedForm.m`, the first check of the code against our own §10.

## 1. The gap

`thm:rhoentry` (ii) reads:

> has absolute value `g/√|A|` for odd `g` and `g/M` for `4 | g`

`g ≡ 2 (mod 4)` is not covered. Those cases are not vacuous — at `X_0^6(1)` and `X_0^10(1)` they
are `g = 2, 6, 10`.

## 2. What is measured

Building `ρ(γ)` from the `S`/`T` matrices for `γ = [[1,0],[c,1]]` and reading the `e_0` column,
**every** nonzero entry at `g ≡ 2 (mod 4)` has

    |entry| = 2g/M          i.e.  |entry|^2 = 4g^2/M^2

8 of 8 such `γ`, across both bases — `6_1` at `c = 2, 6, 10` and `10_1` at `c = 2, 6, 10, 14, 18`.
That is exactly **4×** the `4|g` value, uniformly.

## 3. RESOLVED — the proof's `|A[c]|` has the wrong condition

*(settled 2026-09-05 by computing `|A[c]|` directly as the `c`-torsion of the discriminant group,
`c = 1..20` on both bases.)*

The proof of `thm:rhoentry` states

    |A[c]| = prod_{p | g odd} p^2,  times 8 if 4 | c

Measured against the actual `c`-torsion, the ratio (measured / formula) is

    1   whenever c is ODD or 4 | c
    8   exactly when c = 2 (mod 4)

⇒ **The condition should be `2 | c`, not `4 | c`.** The 2-part of `|A[c]|` is 8 whenever `c` is
even. Sample: `c=2 -> 8`, `c=6 -> 72 = 8*9`, `c=10 -> 200 = 8*25`, against a formula predicting
`1`, `9`, `25`. Odd `c` and `4 | c` are unaffected, which is why the error never showed:
`c=3 -> 9`, `c=4 -> 8`, `c=12 -> 72 = |A|`, `c=20 -> 200 = |A|`, all matching.

## 4. With that fixed, the missing case FOLLOWS

It is no longer an observation needing separate justification. For `c = 2 (mod 4)` write `g = 2k`
with `k` odd; then `prod_{p|g odd} p^2 = k^2 = g^2/4` and the 2-part is 8, so

    |A[c]| = 8 * g^2/4 = 2g^2      and     |entry|^2 = |A[c]|/|A| = 2g^2/(M^2/2) = 4g^2/M^2

giving `|entry| = 2g/M`, **exactly the measured value**. And the two stated cases fall out of the
same corrected rule:

    odd c   : 2-part 1,  |A[c]| = g^2,      |entry|^2 = 2g^2/M^2  ->  |entry| = g/sqrt|A|   (as stated)
    4 | c   : 2-part 8,  |A[c]| = 8(g/4)^2
                                = g^2/2,    |entry|^2 = g^2/M^2   ->  |entry| = g/M        (as stated)
    c = 2(4): 2-part 8,  |A[c]| = 2g^2,     |entry|^2 = 4g^2/M^2  ->  |entry| = 2g/M        (NEW)

All three from one rule, verified on 32 values of `c` across `X_0^6(1)` and `X_0^10(1)`.

## 5. ⚠⚠ HOSTILE REVIEW, 2026-09-05: edit (1) SURVIVES, edit (2) is BROKEN

An adversarial review was commissioned specifically to break this. Verdict: SPLIT.

### Edit (1) — `4|c -> 2|c` — SURVIVES, with stronger support than §3 had

* `A[c]` IS the c-torsion, from Stromberg's own definition (arXiv:1108.0202): `D_c = ker(phi_c)`,
  "all elements of `D` of order dividing `c`", and the coefficient is `xi(A) sqrt(|D_c|/|D|)`.
  That EXCLUDES the `cA` and `A/cA` readings — the attack most likely to kill it.
* Verified at **14 bases**, `c = 1..3M`, computed exactly from invariant factors and
  brute-force cross-checked. Paper formula fails at exactly `c = 2 (mod 4)`, ratio exactly 8, at
  all 13 even-`DN` bases; corrected formula: 0 failures anywhere.
* **The origin of the error**: `A_2 = (Z/2)^3` has **exponent 2** but **level 4**. The proof's
  "(iv) for `4|c` the level-4 component is annihilated" is about the level; the annihilation
  actually happens already at `2|c`.
* **No collateral damage.** Purely 2-adic, so odd-`p` factors, `c_p^Eich/c_p^ram`, `g_p = f_p/2`
  and the telescoping in `thm:noresplit`/`thm:Ngeneral` are untouched. `c_2` enters
  `thm:noresplit` only as a factor common to both lattices and cancels; the one quoted numeric
  2-adic value (`loc_2 = (1+i)/4` at `X_0^{22}(3)`) is at odd `c`, where `|A_2[c]| = 1` either
  way. Production is unaffected regardless — the pipeline computes rho numerically
  (`VVRhoInvE0`), never from this closed form.
* ⚠ Also observed: `|A| = M^2/2` is an **even-`DN`** statement. At `15_1` (odd `DN`)
  `|A| = 450`, not `M^2/2 = 1800`. Harmless for the theorem's stated scope, but worth knowing.

### Edit (2) — the third case — **BROKEN. Do not add it.**

`thm:rhoentry` (ii) is about `(rho_Abar(gamma~)e_0)_{eta*}` with **`eta*` a NONZERO ISOTROPIC
element**. §2's measurement never touched such an entry: `X_0^6(1)` and `X_0^{10}(1)` have **no
nonzero isotropic cosets at all** (1 of 72, 1 of 200 are isotropic — the zero one). The
`|entry| = 2g/M` measured in §2 is the modulus at the **non-isotropic** cosets that carry the
support.

**At the actual `eta*`, the entry is EXACTLY ZERO for `c = 2 (mod 4)`** — `<= 4e-61` at 60 digits,
at every odd-`N` base with genuine `eta*` (`6_5`, `10_3`, `14_3`), and unchanged over 318
representatives with varying `d` and both rho conventions. The reason is Stromberg's Lemma 2.4:
`x_c != 0` exactly when `2 || c` (the odd Jordan component `A_2^1` has scale 2), and then `x_c`
has nonzero 2-part while `cA` has trivial 2-part, so nothing in the support has trivial 2-part —
while every nonzero isotropic `eta` DOES have trivial 2-part when `2 | D`.

⇒ **This is the repo's headline failure mode, committed by me, one step from the paper: the right
number about the wrong object.** §2 measured a real quantity at cosets the theorem is not about.

### A THIRD error the review found, which I was not looking for

**`thm:rhoentry` (i) is also wrong as stated.** "vanishes iff `N | g`" fails at `X_0^6(5)`,
`c = 2` (`g = 2`, `N = 5` does not divide 2), and at `c = 6, 14, 18, 22, 26, 34, 38`; likewise
`10_3` and `14_3`. The correct statement is

    vanishes iff  N | g   OR   c = 2 (mod 4)          (equivalently: N | g or x_c != 0)

and that also explains why (ii) has only two cases — but not for the reason the current proof
gives.

### Adjacent risk, not caused by the correction

§10's `ct_L = glob * prod_p loc_p` with `loc_p = e(-sigma_p/4) G_p sqrt(|A_p[c]|/|A_p|)` is never
zero, and `thm:Ngeneral`'s coset bookkeeping tracks only the odd-part indicator
`prod_{p in S}[p nmid c]`, implicitly assuming `x_c = 0`. At `c = 2 (mod 4)` the truth is 0 (e.g.
`6_1`, `c = 2`: formula `sqrt(8/72) = 1/3`, truth `0`). `rho_mu = kappa rho_0` survives there only
because both sides vanish. The factorisation is valid for `c != 2 (mod 4)`. Deployed numerics are
unaffected. **Worth a look before `sec:determined` is submitted.**

## 6. Proposed edits to `level-prime-kappa.tex`, as amended by the review

1. **KEEP**: in the proof of `thm:rhoentry`, `times $8$ if $4\mid c$` -> `times $8$ if $2\mid c$`,
   justified by `A_2 = (Z/2)^3` having exponent 2 (level 4), and `A[c]` being Stromberg's
   `D_c = ker(\cdot c)`.
2. **DROP** the third case in (ii). Replace with an amendment to **(i)**: vanishes iff `N | g` OR
   `c = 2 (mod 4)`.
3. In the proof of (iv), separate the two facts: `x_c = 0` for `4 | c` is correct, but the
   annihilation of the 2-part (the factor 8) already happens at `2 | c`.
4. `tests/RhoEntryClosedForm.m`: the `g = 2 (mod 4)` skip should become an assertion that the
   **isotropic** entries VANISH there — not that `|entry| = 2g/M`, which holds only at
   non-isotropic cosets. The current comment states the observation in a way that reads as a
   claim about the theorem's entry; it must say which cosets it measured.

⚠ None applied. Edit (1) touches a published-track proof; edits (2)/(3) change a theorem
statement. All four want a second reading.
