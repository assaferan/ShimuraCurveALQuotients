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

⚠⚠ **CORRECTED 2026-09-05 — the `15_2` `2^4` DEFECT NO LONGER EXISTS.** Every earlier version of
this draft (and memory `embedding-selection-root-cause`, 2026-08-16) said `d = −7, −15, −60` at
`X_0^{15}(2)` are each short by `2^4`. **Measured today, they are not.** Recomputed values against
the expected `our_s = their_s / 2`:

    d = -7    ours 1/8    expected 1/8    OK
    d = -15   ours 5/8    expected 5/8    OK
    d = -52   ours 1/2    expected 1/2    OK
    d = -88   ours 2      expected 2      OK
    d = -60   ours -1/24  expected +1/24  WRONG -- SIGN ONLY
    (-12 -> oo, -40 -> 0, -120 -> 1: the three normalising anchors, all correct)

Some intervening fix repaired the magnitude. **The surviving defect at `15_2` is a pure SIGN error
at `d = −60`** — and `gcd(60, 2) = 2`, i.e. non-coprime to the level, exactly like `26_3`'s two bad
discriminants.

⇒ This **invalidates the magnitude premise** these notes were built on, and **strengthens the
conclusion**: the only surviving defects at BOTH bases are sign errors at discriminants
non-coprime to `N`. That is what a sign-valued `p | N` correction predicts, and it is now the
evidence rather than the `2^4`.

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

## 4b. DERIVATION ATTEMPT, 2026-09-05 — what is now settled, and what is not

### (a) The correction is an INDEX factor, not a Whittaker one. [T1/T2, structurally settled]

Decisive, and it rests on a result this repo already proved rather than on a guess. §10 factors the
local contribution as

    loc_p = e(-sigma_p/4) * G_p * sqrt( |A_p[c]| / |A_p| )

and for the Eichler plane `|A_p| = p^2`, with `A_p[c]` trivial when `p ∤ c` — which is precisely
where `c_p^Eich = 1/p` comes from. It is the third factor, the **index**, not the Gauss sum and not
a Whittaker value.

That the Whittaker route is EMPTY is not an assumption:
`tests/KudlaYangLocal.m::test_prop55_nonzero_isotropic_coset` proves (1440 checks, five bases,
including `N = 2`) that at the level prime, for a nonzero isotropic `N`-only-supported coset, KY
Prop 5.5 **collapses to the constant polynomial 1 — no `m`-dependence, no pole, every `m`, every
conductor.** A constant-1 local Whittaker never vanishes, so `p | N` is **never in the Diff set**
and contributes **no `log p` from any derivative**.

⇒ **T1/T2 resolved structurally.** The missing factor is `m`-INDEPENDENT. So it does NOT enter the
`Σ_m c(−m) κ(−m)` sum term by term; it multiplies `Ψ` as a whole, and therefore appears with the
multiplicity the `log|Ψ|²` convention gives it, not with a multiplicity that has to be discovered.

### (b) Where the `d`-dependence comes from — and why only SOME discs are wrong

A factor depending only on `N` would shift every discriminant equally, and it demonstrably does
not: at `15_2` **16 of 19 discs are exact** and only 3 fail. The resolution is in §10's own
statement — `c_p^Eich = 1` when `p | c`, and `1/p` otherwise. The factor depends on **`c`**, the
lower-left entry of the coset representative, hence on the CM point. Discs with `p | c` need no
correction; discs with `p ∤ c` do. That is the mechanism, and it predicts a partition of the discs
of exactly the observed kind.

Consistent with this, `15_2`'s three failing discs (`−7, −15, −60`) are exactly those of **odd
fundamental discriminant**, i.e. those where 2 is split — the case KRY singles out as never
entering the Diff set.

### (c) Magnitude — SUGGESTIVE, NOT DERIVED [T3 still open]

`15_2` is short by exactly `2^4`. With `|A_2| = p^2 = 4`, `|A_2|^2 = 16 = 2^4`, which is what a
factor of `|A_p|` entering `|Ψ|` linearly would give in `|Ψ|^2`.
⚠ **This is numerology until the normalisation is derived.** `sqrt(|A_p[c]|/|A_p|) = 1/2` at
`p = 2` would give `2^{-2}` in `|Ψ|^2`, not `2^{+4}` — wrong sign AND wrong exponent. So the bare
index does not produce the observed factor, and the gap between `1/p` and `p^2` has to be closed by
the actual normalisation (the `|CM(d)|/4` prefactor, the `[GY]`/`[Err]` square-and-quarter
convention, and whether `s` scales as `|Ψ|` or `|Ψ|^2`). **Do not fit this constant.**

### (d) ⚠ The two regression targets may NOT be the same bug

An index factor **scales** a value; it does not **exchange** two values. `15_2`'s defect IS a
scaling (`2^4`, uniform across the three affected discs). `26_3`'s is an **exchange**
(`s ↔ s̃`, with `s + s̃ = 1` preserved exactly). Those are different signatures, and the derivation
above only addresses the first.
⇒ Treat `26_3` as possibly a SECOND defect — most likely T4, the coset bookkeeping: which `S =
{p | N : μ_p ≠ 0}` each Borcherds form is assigned at `p | N`. §10's coset enters as an
**indicator**, `ρ_μ = κ ρ_0` with `κ = ∏_{p∈S}[p ∤ c]`, and an indicator assigned to the wrong form
exchanges two values while preserving their sum. That matches `26_3` exactly and matches `15_2` not
at all.
⚠ Corollary for §5 below: "a formula that fixes one and not the other is wrong" is now itself
**suspect** — if these are two defects, the right fix for one need not touch the other. Keep both
targets, but stop treating agreement on both as the pass condition until (d) is settled.

### (e) THE CONVENTION IS ALREADY PINNED, AND N=1 IS A CONTROL [T3 sharpened]

*(assaferan, 2026-09-05: the Errthum convention was resolved and there is a CI test. Correct.)*
`tests/Schofer.m` is **in CI** and validates `SchoferFormula` against published values —
`[Errthum, p.850]` and `[Errthum, Table 2]` at `D=6, N=1`, `[Err, Table 4]` at `D=10, N=1`.

⇒ **The decisive control: `-267` and `-708` are in that list and PASS at `6_1`.** The very
discriminants that come out exchanged at `26_3` (`N=3`) are correct at `N=1`, against an
independent published source, with the same code. Same `d`, same machinery, only the LEVEL differs.

Two consequences for T3:
* The `N=1` normalisation — the `|CM(d)|/4` prefactor and the square-and-quarter convention — is
  **already validated**. So the missing constant is **purely the `p|N` contribution** and nothing
  else. That is a far sharper target than "reconcile the normalisations".
* **Any fix MUST be inert at `N=1`**, or it breaks a CI test against published values. §10's
  structure satisfies this automatically: `prod_{p|N}` is empty when `N=1`. That is a real
  constraint on the answer, and a free regression test.

⇒ **Third regression target: `6_1` and `10_1` must keep passing** (`tests/Schofer.m`). Unlike the
other two this one is already automated.

### (f) ⚠ A lead: the `10_1`, `d = -68` "error in [Err]" may be OURS

`tests/Schofer.m`, `D=10` block, carries:

    // The next one does still not work at 2! Off by a factor of 4, this is confusing.
    // Probably an error in [Err] ?
    // <-68, { <2, 2>, <5, 1> }>,     <- Errthum's published value
       <-68, { <2, 6>, <5, 1> }>,     <- ours, committed as the expected value

The discrepancy is **4 in the `log 2` coefficient** — the SAME `2^4` signature as `15_2`'s three
failing discs. It is currently attributed to a mistake in the source. Given we are now hunting a
`2^4`-shaped index factor at exactly this prime, **that attribution should be re-examined**: if it
is our defect rather than Errthum's, then it appears at `N = 1` too and the whole "it is a level
effect" framing needs revisiting.
⚠ Note this cuts BOTH ways and is why it matters: `d = -68 = -4*17` has `2 | d` and `2 | D = 10`,
so 2 is ramified — a different local situation from the `p | N` Eichler case. Do not assume the two
`2^4`s have one cause; check.

### (g) ⇒ THE CORRECTION IS A SIGN, NOT A MAGNITUDE. [T1 and T2 RESOLVED]

This inverts §4b(c) above, which chased a magnitude and could not make one work. Two facts from
`thm:rhoentry` and §10 settle it.

**The Weil-representation coefficient.** `thm:rhoentry`'s proof quotes Stromberg Thm 6.4: the
`eta*`-component of `rho_Abar(gamma~)e_0` is the SINGLE term

    xi(a,c) * sqrt( |A[c]| / |A| ) * e( ac*qbar(y_0) + B(x_c,y_0) ),      |A| = M^2/2

**Fact 1 — the index cannot distinguish the cases.** `15_1` and `15_2` have the SAME `M = 60`,
hence the same `|A| = M^2/2`; and `|A_p| = p^2` for the anisotropic plane and the hyperbolic plane
alike. So `sqrt(|A[c]|/|A|)` is identical in the ramified and Eichler cases. **No index factor can
produce a magnitude change**, and §4b(c)'s hunt for one was misconceived.

**Fact 2 — §10 says the difference is exactly a sign.** `c_p^ram = eps_p c_p^Eich`, `eps_p = ±1`:
Gauss sum `+p` vs `-p`, signature `0` vs `4`, `|c_p| = 1/p` in both. The Eichler and ramified local
factors differ ONLY in sign.

**⇒ The consequence, and it explains the `2^4` exactly.** In a log-linear identity, using `-x`
where `+x` belongs shifts the answer by `2x`, not by `x`. So:

⚠ **The `2^4` arithmetic that motivated this is now void** — see §1: `15_2`'s magnitudes are
correct today and only `d = −60`'s SIGN is wrong. The conclusion survives on stronger evidence:
the surviving defects at `15_2` (`−60`) and `26_3` (`−267`, `−708`) are ALL sign errors, and ALL at
discriminants non-coprime to `N`. A sign-valued correction at `p | N` predicts exactly that; a
magnitude-valued one predicts nothing of the kind.

⇒ **T1 RESOLVED**: the correction is a sign, and its direction is fixed by §10 — at `p | N` the
form is the HYPERBOLIC plane, so `c_p^Eich = +1/p`; using the ramified `-1/p` is the error.
⇒ **T2 RESOLVED**: the "multiplicity 4" was never a multiplicity. It is `2 x (a coefficient-2
term)`, the factor 2 being the sign-flip doubling. Nothing has to attach four times.
⇒ **T3 DISSOLVED**: there is no magnitude constant to derive from KY/KRY. The normalisation is
already validated at `N = 1` (§4b(e)), and the `p | N` correction contributes no magnitude at all.

**This also explains `26_3` — so §4b(d) is superseded.** A sign error is not a scaling, so the two
targets need not have different causes after all. At `26_3` the two hauptmoduls are EXCHANGED with
`s + s~ = 1` preserved; a sign error in the Weil-representation phase at `p | N` selects the wrong
branch of exactly that pair. The `15_2` and `26_3` signatures are the same defect seen through
different forms.

### (h) The concrete, checkable prediction

**At `p | N` the code should be using the HYPERBOLIC plane (`Qbar(b,c) = -bc/p`, from
`Q = -x^2 - pyz`, Gauss sum `+p`, signature `0`), not the anisotropic/ramified one.** If it selects
the local form by `p | D*N` rather than distinguishing `p | D` from `p | N`, that is the bug, and it
is a sign error of exactly the derived shape.
⇒ NEXT: read the local-form selection in the Whittaker/Jordan path (`Wpoly_scaled`, `Wpoly2`'s
Jordan classification, and `get_lattice_data`) and check whether `p | N` is distinguished from
`p | D`. This is now a code question with a derived answer to check against, not a search.

⚠ Still to verify before implementing: that a term of coefficient exactly `2 log 2` is present and
sign-flipped at `15_2`'s three discs — i.e. confirm the mechanism at the level of the actual terms,
not just the totals. The totals matching is necessary, not sufficient.

### (i) ⚠ THE PREDICTION OF (g)/(h) IS NOT CONFIRMED BY TERM-LEVEL PROBING [2026-09-05]

The sign derivation of §4b(g) predicts a sign error in the `p | N` local factor. Probed directly at
an invariant-identical pair — `d = -267` (BAD) against `d = -123` (good), both `h=2`, both
non-coprime to `N=3`, one CM point each — **the `p = 3` term is identical in sign and magnitude on
both sides of the formula**:

    oo-side (Kappa0, mu = 0):
      m=6   -6Log3 - 6/7Log13        |  -6Log3 - 6/7Log13        p=3 term IDENTICAL
    cusp-0 side (sum over nonzero gammas in the mM bucket):
      mM=3  -2Log3 - 8Log5           |  -2Log3 - 4Log127         p=3 term IDENTICAL
      mM=4  -4/3Log2                 |  -4/3Log2                 identical

The other differences (`Log5`, `Log7`, `Log127`, ...) are the ordinary arithmetic of different
discriminants. **No sign anomaly at the level prime appears anywhere.**

⚠ Also ruled out by code reading: the local form is NOT selected by a `p | D*N` test. `Wpoly_scaled`
builds `S` from the lattice and hands it to `Wpoly`/`Wpoly2`, which take `eps` from the Jordan block
diagonals and the sign from `LegendreSymbol(-1,p)^(l_mu(k) div 2) * prod legendre(eps[i],p)`. The
hyperbolic-vs-anisotropic distinction is read off the lattice, correctly, by construction. So §4b(h)'s
"check whether it selects on `p | D*N`" is answered: it does not, and that is not the bug.

### (j) A STRUCTURAL CONSTRAINT any explanation of the `26_3` exchange must satisfy

Both sides of the formula compute a **form-independent kernel**:
* `SchoferFormula` (oo-side) calls `Kappa0` only — `mu = 0` for EVERY form;
* `SchoferFormula0` (cusp-0) buckets by `mM mod M` and sums `Kappa(gamma, mM/M, d, ...)` over that
  bucket — the bucket depends on `mM`, NOT on the form.

The form enters only as `Coefficient(f, -mM)`. So every form's value is `sum_m c_f(-m) K(m)` with a
COMMON `K(m)`. Therefore an exchange of two forms' values requires either

  (A) `K(m)` wrong at exactly those discriminants — but the `p=3` part of `K` is measured identical
      to the good twin's, so if this is it, the error is at some other prime and is not a `p | N`
      effect at all; or
  (B) the coefficient vectors `c_1, c_2` swapped — but that is `d`-INDEPENDENT and would corrupt all
      14 discriminants, not 2.

**Neither is consistent with what is measured.** That is the sharpest statement of the puzzle
currently available, and it means the `26_3` exchange is NOT yet explained by this derivation.

⚠ **Do not implement the §4b(g) sign fix on the strength of the derivation alone.** It is a clean
argument whose empirical prediction fails at the only place it can be tested. Either the derivation
applies somewhere the probes did not reach, or it is wrong. Settle that before writing code.

## 5. Two independent regression targets, with exact known answers

Any candidate translation must hit **both**, which is what makes it hard to fit a wrong formula:

* `26_3`, `N = 3`: Guo-Yang's CM table, all 14 discriminants, with `−267 → 17/25` and
  `−708 → 38/49` the two currently wrong.
* `15_2`, `N = 2`: **`d = −60` only**, and it is a SIGN error (`−1/24`, expected `+1/24`).
  ⚠ `−7` and `−15` are NO LONGER failing — do not use them as targets; the `2^4` shortfall recorded
  for them is stale.

Different level primes, and BOTH failure shapes are now sign-valued (an exchange at `26_3`, a sign
flip at `15_2`). A correct sign rule should fix both; a magnitude constant would fix neither.
⚠ Both bases also share the pattern that only SOME non-coprime discriminants fail — 1 of 4 at
`15_2`, 2 of 11 at `26_3` — so the rule must be finer than "non-coprime to `N`".

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
