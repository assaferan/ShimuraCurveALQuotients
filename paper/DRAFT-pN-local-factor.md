# DRAFT — the `p | N` local factor on the Schofer side

**Status: a derivation TARGET, not a derivation.** This states precisely what has to be proved
before anything is implemented, and records exactly which ambiguities block writing the formula
down today. Drafted 2026-09-05; not yet part of `level-prime-kappa.tex`.

Written because the repo's own rule is to derive first
(memory `feedback-derive-before-implementing`): *"only when we have something written down that
precisely describes the formula we want to implement, should we implement it."* Guessing the
conventions and fitting numerically is exactly what produces a plausible wrong formula.

---

## 1. What is broken

`X_0^{26}(3)`, Guo-Yang's CM table, our computed hauptmodul values:

| disc | ours | correct | relation |
|---|---|---|---|
| −267 | 8/25 | 17/25 | ours = 1 − correct |
| −708 | 11/49 | 38/49 | ours = 1 − correct |
| the other 12 | correct | | |

`s + s̃ = 1` holds at every disc, so the pair is *consistent* and the two hauptmoduls' values are
**exchanged** at exactly these two. Both are non-coprime to `N = 3` — but so are nine discs that
come out right, so coprimality alone does not predict it.

At `X_0^{15}(2)` the same class of defect appears differently: `d = −7, −15, −60` are each short by
exactly `2^4`, i.e. the correct value needs **+4·log 2**.

## 2. Why it is a `p | N` problem

Schofer's Thm 4.1 rests on the assumption stated verbatim on p.51 of the thesis: *"For an
unramified prime `p`, the local lattice `A_p = A ⊗ Z_p` is UNIMODULAR."* That is true for a CM-field
ideal at level 1. It **fails at `p | N`**: the Eichler order of level `N` makes `L^-` non-unimodular
at `p | N` even though `p` is unramified in the quaternion algebra. Verified directly: `L^- = λ^⊥`
at 2 for `d = −7` is unimodular on `15_1` (`v_2(det) = 0`) and **2-modular** on `15_2`
(`v_2(det) = 2`).

## 3. What §10 supplies

`sec:determined`, in the proof of `thm:noresplit`: at each `p | DN` the local discriminant group is
`(Z/p)^2` and the two relevant forms are

| | form | from | Gauss sum | signature | `c_p` (`p ∤ c`) |
|---|---|---|---|---|---|
| Eichler level `p` (unramified, **non-unimodular**) | hyperbolic `Q̄(b,c) = −bc/p` | `Q = −x² − pyz` | `p` | 0 | `+1/p` |
| ramified (`p | D`) | anisotropic `(uc² − b²)/p` | `nrd = −ux² − py² + upz²` | `−p` | 4 | `−1/p` |

(both are `1` when `p | c`), whence `c_p^ram = ε_p c_p^Eich` and, by `thm:Ngeneral`,

```
rho_0 = prod_{p|D} c_p^ram * prod_{p|N} c_p^Eich
```

⇒ **This is the missing local structure.** The `15_2` investigation ended by concluding the missing
`+4 log 2` "must come from a NON-kappa level-`N` bad-place factor that Schofer Thm 4.1 LACKS", and
identified the dead `kappaminuszero` (`SchoferFormula.m:569`, `log_coeffs[p] := 1` for
`p | N/gcd(d,N)`) as its vestige. §10 is what that vestige was waiting for.

## 4. The translation problem, stated

§10 gives a **multiplicative** factor on `ct_L`. Schofer gives an **additive** identity in logs,

```
log |Psi_f(z_d)|^2  =  -(|CM(d)|/4) * sum_m c_f(-m) kappa(-m),
```

where `kappa` picks up `log p` only from the *derivative* of the local Whittaker at the unique
prime of the Diff set. The theorem to prove is the bridge. **Four things are genuinely undetermined
by §10 as written**, and each one is a place where a numeric fit would silently choose wrong:

**(T1) Sign.** `|c_p^Eich| = 1/p`, so the naive reading gives `log|c_p| = −log p`. `15_2` needs
`+4 log 2`. The sign is wrong before multiplicity is even considered.

**(T2) Multiplicity.** `15_2` needs the factor **four** times over. Does it attach once per CM
value, once per `m` in `Σ_m c(−m) κ(−m)`, or once per term of the representation sum inside
`kappa`? §10 has no `m` in it, so it cannot answer this.

**(T3) Normalisation.** Schofer computes `log|Ψ|²` with the `|CM(d)|/4` prefactor and the
`[GY]`/`[Err]` square-and-quarter convention this repo already documents as delicate; `ρ_0` is a
bare constant term. The two must be reconciled explicitly, not by matching magnitudes.

**(T4) Coset bookkeeping — the likely source of the `26_3` EXCHANGE.** §10's coset enters as an
**indicator, not a phase**: `rho_mu = kappa * rho_0`, `kappa = prod_{p in S} [p ∤ c]`, with
`S = {p | N : mu_p ≠ 0}`. Note `c_p^Eich = −c_p^ram` differ **only in sign**, so no Eich↔ram
substitution can produce a magnitude change like `2^4`; and conversely a *sign*/indicator error at
`p | N` is exactly the shape that would exchange the values of two forms while preserving
`s + s̃ = 1`. Which `S` each Borcherds form sees at `p | N` has to be pinned.

## 5. Two independent regression targets, with exact known answers

Any candidate translation must hit **both**, which is what makes it hard to fit a wrong formula:

* `26_3`, `N = 3`: Guo-Yang's CM table, all 14 discriminants, with `−267 → 17/25` and
  `−708 → 38/49` the two currently wrong.
* `15_2`, `N = 2`: `d = −7, −15, −60`, each currently short by exactly `2^4`
  (`1/128` vs `1/8`, `5/128` vs `5/8`, `1/384` vs `1/24`).

Different level primes, different failure shapes (exchange vs scaling). A formula that fixes one
and not the other is wrong.

## 6. What is already ruled out — do not re-derive

* point selection — every disc contributes exactly ONE CM point (`-267` h=2 and `-123` h=2 alike)
* `ScaleForSchofer` — uniformly `−1/4` across the bad/good pairs
* the optimal embedding `λ` — all of correct norm `−2d`, `v_3 = 1`, no structural difference
* `find_signs_hauptmodul` — never reports ambiguity at the bad discs
* seed / quaternion-order representative — stable across 3 seeds; PR #18's isometry fix holds
  (`ElementOfNorm` returns different `λ` per run while the values stay fixed)
* any invariant of `d` — `−267` matches the GOOD `−123`, `−708` matches the GOOD `−132`, on class
  number, fundamentality, every valuation and every Kronecker symbol. **So a `Keep` list keyed on
  `d` cannot work.**

## 7. Recommendation

Resolve **T1–T4 on paper**, as a numbered statement in `level-prime-kappa.tex` (natural home: a
subsection of `sec:determined`, or a short section after it, since it consumes `thm:Ngeneral`'s
local values). Only then implement, treating that statement as the spec — and test against §5's
two targets before believing it.

⚠ Do **not** start by fitting a constant to make `15_2`'s `2^4` come out. That is precisely the
move the derive-first rule exists to prevent, and with four free choices (T1–T4) a fit will
succeed while being wrong.
