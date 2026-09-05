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

## 5. Proposed edits to `level-prime-kappa.tex`

Both are small and independent of everything else in the paper:

1. **In the proof of `thm:rhoentry`**, change `times $8$ if $4\mid c$` to `times $8$ if $2\mid c$`.
2. **In `thm:rhoentry` (ii)**, add the third case: `and $2g/M$ for $g\equiv2\pmod 4$`.

⚠ Neither is applied — this is a draft, per the repo's derive-and-draft habit. Edit (1) is a
correction to a published-track proof, so it deserves a second reading before it lands; note the
error is invisible in the two cases the theorem states, which is presumably why it survived.

## 6. Provenance

`tests/RhoEntryClosedForm.m` currently SKIPS `g = 2 (mod 4)` rather than asserting there. Once
edit (2) lands, that skip should be removed and the case asserted -- the test already computes
everything needed. An earlier version applied the `4|g` formula to these cases and reported 8
"failures"; that was the test's scope error, not a divergence.
