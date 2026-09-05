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

## 3. ⚠ Why this is NOT yet proposed as theorem text

The proof of `thm:rhoentry` computes the absolute value as `√(|A[c]|/|A|)` with

    |A[c]| = prod_{p | g odd} p^2   times 8 if 4 | c,        |A| = M^2/2

Applying that to `g ≡ 2 (mod 4)`: write `g = 2k` with `k` odd, so `prod_{p|g odd} p^2 = k^2 = g^2/4`,
and `4 ∤ c`, so the factor 8 does not apply. Then

    |entry|^2 = (g^2/4) / (M^2/2) = g^2/(2 M^2)

but the **measurement is `4g^2/M^2`** — larger by a factor of **8**. Note 8 is exactly the
`4 | c` factor, so the discrepancy looks like that factor being present when the stated condition
says it should not be (or the condition being about something other than `4 | c`).

⇒ Either the `|A[c]|` formula needs a case for `c ≡ 2 (mod 4)`, or my reading of it is wrong.
**Until that is settled the observation should not go into the theorem** — adding a third case
whose value contradicts the proof's own machinery would make the paper inconsistent with itself.

## 4. What would settle it

Compute `|A[c]|` directly — the `c`-torsion of the discriminant group — for the same `γ` the test
uses, and compare against both the stated formula and `4g^2/M^2 · |A|`. That is a small
computation on top of `tests/RhoEntryClosedForm.m`, which already has the lattice and the
matrices in hand. If the direct `|A[c]|` matches the measurement, the formula in the proof needs
the extra case; if it matches the formula, the discrepancy is elsewhere in the chain and the
absolute-value claim needs re-deriving for this residue class.

## 5. Provenance

`tests/RhoEntryClosedForm.m` currently **skips** `g ≡ 2 (mod 4)` rather than asserting anything
there, with a comment pointing here. An earlier version of that test applied the `4|g` formula to
these cases and reported 8 "failures"; those were the test's scope error, not a divergence between
code and paper.
