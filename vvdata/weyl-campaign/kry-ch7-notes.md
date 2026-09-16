# KRY Chapter 7 ("An inner product formula") — read 2026-09-16, and why it's retired as a deficit predictor

Source: S. S. Kudla, M. Rapoport, T. Yang, *Modular Forms and Special Cycles on Shimura Curves*
(Annals of Math. Studies 161, Princeton, 2006).
PDF: https://www.math.uni-bonn.de/people/rapoport/myalggeom/preprints/kry.pdf (1.87 MB, fetched
2026-09-16; `pdftotext -layout` gives a 21303-line text dump — kept locally in scratch, not
committed, re-fetch if needed).

**Scope actually read**: §7.1 ("Statement of the main result") in full, §7.6 ("Contributions for
`p | D(B)`") in full. §7.2–7.5 and §7.7–7.11 were located and their section headers read (so the
overall shape of the chapter's argument is known) but NOT read line by line. Treat this as a
strong first pass, not an exhaustive survey — see `PLAN.md`.

## What Theorem C says

`M` is the integral model over `Z` of a Shimura curve for an indefinite quaternion algebra `B` of
discriminant `D(B)`. `φ̂1(τ) = Σ_t Ẑ(t,v) q^t` is a weight-3/2 generating series valued in the
arithmetic Chow group `ĈH^1(M)` (Ch. 4). `φ̂2(τ)` is a genus-2 analogue valued in `ĈH^2` for
0-cycles (Ch. 6), a Siegel modular form of weight 3/2 in each variable after pulling back along the
diagonal `H × H → H_2`.

> **Theorem C.** `⟨ φ̂1(τ1), φ̂1(τ2) ⟩ = φ̂2( diag(τ1,τ2) )`.

Equivalently, for `t1, t2 ∈ Z`, `v1, v2 > 0`:

> `⟨ Ẑ(t1,v1), Ẑ(t2,v2) ⟩ = Σ_{T: diag(T)=(t1,t2)} Ẑ(T, diag(v1,v2))`,  `T` ranging over
> `Sym2(Z)^∨ = { T ∈ Sym2(Q) : tr(Tb) ∈ Z ∀ b ∈ Sym2(Z) }`.

The left side `⟨ , ⟩` is the **arithmetic height pairing** — a real number coming from Arakelov
intersection theory on `M` (Gillet–Soulé style `ĈH^1 × ĈH^1 → R`), built by decomposing each cycle
into horizontal + vertical parts and summing local contributions at every prime plus an archimedean
(Green-function) term. This is a fundamentally different kind of object from a q-expansion
coefficient: it needs the actual integral model, its special fibers, and metrized line bundles.

Theorem 7.1.1 (a corollary, the `t1=t2=0` boundary case) evaluates `⟨ω̂,ω̂⟩` — the self-height of the
Hodge bundle — as `ζ_{D(B)}(-1) [ 2ζ'(-1)/ζ(-1) + 1 - 2C - (1/2)Σ_{p|D(B)} (p+1)/(p-1) log p ]`,
`2C = log(4π)+γ`. This is the kind of closed form Ch. 7 produces: an explicit REAL NUMBER for a
specific, fixed height, not a rank statement about a family of forms.

## The local formula at ramified primes (§7.6, the part closest to our question)

For `p | D(B)` and `T` with `diag(T) = (t1,t2)`, `det(T) ≠ 0`, the local contribution factors
through intersection numbers of special cycles `Z(j1), Z(j2)` in Drinfeld space, via
`ν_p(T) = 2χ(Z(j), O_{Z(j1)} ⊗^L O_{Z(j2)})` (twice the intersection number). Proposition 7.6.4
gives the explicit closed form for the SINGULAR `T` case (`det T = 0`, i.e. `t1 t2 = m^2`): writing
`t1 = n1 t`, `t2 = n2 t`, `gcd(n1,n2)=1`, `4t = n^2 d` with `-d` the discriminant of
`k = Q(√-t)`, `k := ord_p(n)`, `χ := χ_d(p)`:

* `p` ramified or inert in `k`: `C_p(T) = (1/2) deg(Z(t)_Q) · ν̃_p(T) · log p`, with
  `ν̃_p(T) = ord_p(t1 t2) + 2(1 + χ - ord_p(d/4)) - [ (p+1)(p^k - 1)/(p-1) if inert,
  2(p^{k+1}-1)/(p-1) if ramified ]`.
* `p` split in `k`: `C_p(T) = δ(d; D(B)/p) · H_0(t; D(B)) · ν̃_p(T) · log p`, with
  `ν̃_p(T) = -2(p^k - 1)`.

This is genuinely a "local density at the ramified prime" formula, structurally in the same family
as this project's own `κ_p`/`SchoferFormula.m` (Yang-style local densities), and it is finite and
computable for any fixed `(t1,t2)` — **but it computes an intersection multiplicity for one fixed
pair, not a rank statement over a basis.** `T` here is a single quadratic-form index (one point in
the sum defining the height pairing); it never plays the role our target set `T` plays in
`deficit.m` (a FIXED collection of many discriminants, against which a whole basis of `S_{3/2}` is
tested for degeneracy).

## Why this doesn't transfer to `deficit.m`

Re-read `vvdata/weyl-campaign/deficit.m` (already committed, unchanged) to pin down exactly what it
computes: `Ncols(mat) - Rank(ech_basis * mat)`, where `ech_basis` is an echelonized weakly-
holomorphic basis (truncated q-expansion coefficients up to a pole order `P`) and `mat` is
`coeffs_to_divisor_matrix(-P, D, N, ncols)` — a matrix built purely from q-expansion coefficients
and a fixed list of target discriminants. **No scheme, no integral model, no Green functions, no
heights, no archimedean data anywhere in this computation.** It is finite-dimensional linear algebra
over `Q`.

KRY Ch. 7's closest global quantity (the height pairing itself) is, by the book's own framing
(culminating in Ch. 9, "Central derivatives of L-functions"), the Shimura-curve analogue of the
Gross–Zagier formula: heights = central derivatives of `L`-functions. That is exactly the kind of
"hard arithmetic" (Waldspurger-type nonvanishing) object the 2026-09-16 (earlier) session already
identified as the true nature of `deficit.m`'s rank question, when it found that
`dim M_{3/2}(ρ_L^*)` (Riemann–Roch) computes the wrong thing because the real question is whether
SPECIFIC Fourier coefficients vanish at the SPECIFIC discriminants in the small target `T`.

⇒ **KRY Ch. 7 does not give a shortcut around that difficulty — its own central theorem is proved
BY relating heights to central derivatives of `L`-functions, which is the hard side of the question,
not a local closed form that bypasses it.** Retiring this as a deficit-predictor lead. See
`PLAN.md`'s "⚠ THE KRY LEAD" section for the recorded verdict and what, if anything, might still be
worth trying (a representation-theoretic vanishing argument for the SPECIFIC targets that occur,
not a generic formula).
