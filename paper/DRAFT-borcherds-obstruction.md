# Why a base is obstructed: it is Borcherds' obstruction theorem, and the space is `S_{3/2}`

*Drafted 2026-09-15 from the primary sources. Everything in the first three sections is quoted or
directly derived; §4 is the open implementation and §5 is explicitly conjectural.*

## 1. The theorem

Borcherds, *The Gross–Kohnen–Zagier theorem in higher dimensions*, Duke Math. J. **97** (1999),
219–233, Theorem 3.1, verbatim:

> Suppose that `k ∈ ½Z`, `Γ` is a subgroup of `Mp_2(R)` which is commensurable with `Mp_2(Z)` and
> `ρ` is a finite dimensional complex representation of `Γ` factoring through a finite quotient of
> `Γ` such that `ρ = σ_k` on `Γ ∩ K`. Then the space of obstructions `Obstruct(Γ, 2−k, ρ)` is finite
> dimensional and **dual to the space `HolModForm(Γ, k, ρ*)`**.

This is Serre duality on the modular curve: the obstructions to finding a weakly holomorphic form of
weight `2−k` with prescribed singularities are dual to the *holomorphic* forms of weight `k` in the
**dual** representation. There is a published erratum (Duke **105** (2000), 183–184) — check it
before quoting a numbered statement in the paper.

⇒ **Borcherds' map is not unconditionally surjective.** It is surjective onto exactly the divisors
annihilated by that dual space. This is the answer to "why is this base obstructed".

## 2. It specialises to OUR setting, in Borcherds' own words

Example 5.3 of the same paper is this project:

> Suppose that `K` is a non split 4 dimensional central simple algebra over the rationals which is
> split at infinity and `R` is some order in `K`. We let `M` be the lattice of elements of `R`
> orthogonal to 1, with the inner product given by minus that of `R`. Then `M` is a lattice in
> `R^{2,1}`, and the group of units of `R` acts on `M` by conjugation. The Grassmannian `G(M)` is
> isomorphic to the upper half plane `H` and the group `Γ` acts on it with the quotient being a
> compact Riemann surface (or more precisely an orbifold) called a **Shimura curve**. The points on
> the Shimura curve associated to vectors of `M` are called **CM points**. So theorem 4.5 implies
> that certain divisors associated to CM points are coefficients of **modular forms of weight 3/2**.

Our `L` is exactly that `M` (`QuaternionLatticeData`, signature `(b+, b−) = (1,2)` — see
`WeilRepresentation.m`; Borcherds writes the opposite sign convention, `R^{2,1}`). Our input forms
are weight `1 − n/2 = 1/2` with `n = 1`, which is what `thm:respair` in `level-prime-kappa.tex`
already assumes. So `2 − k = 1/2` gives **`k = 3/2`**.

## 3. Holomorphic or cuspidal? In OUR dimension it is cuspidal.

Theorem 3.1 says *holomorphic* forms, not cusp forms; the two differ by whether the constant term of
the prescribed singularity is constrained. ⚠ This distinction is easy to gloss and it matters.
Borcherds settles it in Example 5.4:

> In the case of the Gross–Kohnen–Zagier theorem the modular form we get is a **cusp form** and
> `y_{00} = 0`. We give an example to show that **in higher dimensions** the modular forms we get are
> not necessarily cusp forms and `y_{00}` can be nonzero.

GKZ is the signature-`(2,1)` case, i.e. ours. So the holomorphic-vs-cuspidal gap is a *higher*
-dimensional phenomenon and here the obstruction space is genuinely the cusp forms.

Corroborated independently by Bruinier, *On the converse theorem for Borcherds products*
(arXiv:1210.4821), which builds the map

> `Λ : S_{1+n/2, L} → H^{1,1}(X_Γ)` from the space of cusp forms of weight `1 + n/2` with
> representation `ρ_L` for the group `Mp_2(Z)`

— with `n = 1` that is `S_{3/2, L}`, the same space `thm:respair` pairs against.

⚠ Bruinier's *converse* theorems (which meromorphic forms are Borcherds lifts) assume `n ≥ 3`, and
he notes that "in the case of lattices with signature `(1,2)` there are orthogonal modular forms
which cannot be obtained as a Borcherds lift". That is a different statement from the obstruction
theorem for *divisors* and does not weaken §1–§2; do not conflate them.

## 4. What this makes of `deficit.m` — and the prediction it yields

`vvdata/weyl-campaign/deficit.m` computes `Ncols − Rank` of `coeffs_to_divisor_matrix` applied to
the weakly holomorphic basis truncated at pole order `P`: the dimension of the **cokernel of
Borcherds' map at finite `P`**, restricted to the tracked divisor classes. By §1–§3 that cokernel is
dual to the part of `S_{3/2}(ρ_L^*)` pairing nontrivially with those classes. So the deficit is not
a proxy for the obstruction — **it is the obstruction, measured numerically**, and:

    deficit(P)  ≤  dim S_{3/2}(ρ_L^*)   for every P,   with equality once P is large enough.

Three consequences, in increasing order of usefulness:

1. the `deficit = 1 / 2 / 3` stratification is the dimension of a cusp space, not an empirical
   curiosity, and the "2-dimensional obstruction space" bases acquire a meaning;
2. **a measured deficit ABOVE the bound is an artifact of the ladder** — which decides, by a closed
   form instead of hours of compute, the open question of whether the deficit-2/3 readings at
   `D ≈ 300–630` are real (the pole list in `deficit.m` is hard-coded and caps at `P = 266`,
   while every base the screen was validated on had a flat ladder from `P = 102`);
3. obstruction becomes computable without running anything.

## 5. ⚠ RESOLVED (2026-09-16): `dim S_{3/2}(ρ_L^*)` is now computable, and it is the WRONG target

This section originally asked for `dim S_{3/2}(ρ_L^*)` in closed form as a predictor for `deficit.m`'s
measured number. Both halves of that plan are now settled — the dimension is computable, and it does
**not** predict the deficit, for a reason that is itself informative.

### 5a. Route B (the naive guess) is refuted by measurement

The first candidate, motivated by Eichler/Shimizu/Shimura, was `deficit = genus(X_0^D(N))` (the
weight-2 space on the full Shimura curve itself, `GenusShimuraCurve` in `ShimuraQuotients.m`):

    38_5  : genus  9   deficit 1   match=false
    146_1 : genus  7   deficit 0   match=false
    194_1 : genus  9   deficit 0   match=false   <- SAME genus as 38_5, different deficit

`146_1` and `194_1` share a genus but not a deficit, so deficit cannot be a function of the curve's
genus alone. Refuted immediately, no further work needed on this route.

### 5b. The real Riemann–Roch formula, sourced and validated

The correct general statement is Borcherds' own (*GKZ in higher dimensions*, Duke 97 (1999), p. 9,
right after Lemma 4.4 — not the Bruinier/Kuss citation this section originally guessed at): for a
`d`-dimensional representation `ρ` of `Mp₂(Z)` on which the metaplectic central element `Z = S²` acts
as `e^{−iπk}·Id`,

    dim HolModForm(ρ, k) = d + dk/12 − α(e^{iπk/2}S) − α((e^{iπk/3}ST)^{−1}) − α(T)

where `α(X)` sums the fractional eigenvalue-phases of `X`. Applying it to `ρ = ρ_L^*` at `k = 3/2`
needs eigenvalues of `ρ_L^*(S)` and `ρ_L^*(ST)` restricted to the unique subspace where `Z` is
compatible with `k = 3/2` — the "antisymmetric under `γ ↦ −γ`" eigenspace of the dual, confirmed by
direct computation of `Z = ρ_L^*(S)²` on both eigenspaces, not assumed.

**The obstacle was scale, not theory**: this repo's own `WeilRepresentationST` builds the honest
`|L'/L| × |L'/L|` matrix, and `|L'/L|` runs from 72,200 to 1,299,272 on the calibration bases —
far beyond what can be diagonalized. The fix: every quantity Borcherds' formula needs reduces to a
handful of `O(n)` Gauss sums `Σ_γ e(c·Q(γ))` (never the full matrix), via two algebraic tricks —
`tr(S²·X)` isolates one eigenspace of the negation involution by a projector trick, and `tr((ST)²)`
factors into a *product* of two single Gauss sums because shifting `γ ↦ γ − 2δ` is a bijection of the
group for fixed `δ`. **All six resulting trace identities were checked EXACTLY (bit-for-bit, in the
same cyclotomic field) against the real matrices** on `6_1` (n=72) and `10_1` (n=200) before being
trusted on anything larger — this is the "reproduce a known value first" habit, applied to a formula
rather than a number.

### 5c. Computed, and it answers a different question than deficit

    38_5   dim M_{3/2}(ρ_L^*) = 1594     (deficit.m measures 1)
    146_1  dim M_{3/2}(ρ_L^*) =  888     (deficit.m measures 0)

Both in the thousands, both wildly larger than the tracked deficit. **This is not a bug in the
formula — it is the correct dimension of the wrong space.** `deficit.m` does not measure
`dim S_{3/2}(ρ_L^*)`; Serre duality (§1) says the obstruction to prescribing an *arbitrary* principal
part is dual to the *whole* space `S_{3/2}(ρ_L^*)`, but the pipeline only ever asks to hit a small,
fixed target set `T` of specific CM-divisor classes (the handful of cosets the model construction
actually needs). The right statement is

    deficit  =  rank of the pairing  S_{3/2}(ρ_L^*) → T*   ≤  min(dim S_{3/2}(ρ_L^*), dim T)

so the deficit is bounded by `dim T` (small, by construction — a dozen or so classes), not by the
ambient dimension, and the Riemann–Roch number above is essentially irrelevant to it: making the
haystack bigger does not change whether this specific needle is in it.

### 5d. Why this resists a closed form — and it is not merely unfinished work

Whether the pairing degenerates on `T` depends on whether *specific Fourier coefficients* of a basis
of `S_{3/2}(ρ_L^*)` vanish at the *specific discriminants* in `T`. That is a "hard" arithmetic
question (Waldspurger's theorem ties half-integral-weight coefficients to central values of the
Shimura/Shintani lift), not a "soft" dimension-counting one — no Riemann–Roch or trace-formula
argument can answer it, by construction, since those only ever see the size of a space, not which
specific coefficients within it vanish. This is consistent with the data: most bases screen clear
(the generic, full-rank outcome), and the 132 known-obstructed bases are presumably each sitting on
a specific algebraic coincidence rather than a smooth trend in `(D,N)`.

**What survives**: a genuine, validated, `O(n)`-time computation of `dim M_{3/2}(ρ_L^*)` (equivalently
`S_{3/2}`, since weight 3/2 is cuspidal here by §3) — a real invariant, and a correct upper bound on
the deficit, banked as `vvdata/weyl-campaign/gksz-dim-formula.m` on the campaign branch. It is not,
and per §5d cannot be turned into, a predictor for the small number `deficit.m` actually measures.

## Sources

* R. E. Borcherds, *The Gross–Kohnen–Zagier theorem in higher dimensions*, Duke Math. J. 97 (1999)
  219–233 (+ erratum, Duke 105 (2000) 183–184). `https://math.berkeley.edu/~reb/papers/gkz/gkz.pdf`
* J. H. Bruinier, *On the converse theorem for Borcherds products*, arXiv:1210.4821.
* J. H. Bruinier, *Borcherds products on O(2,l) and Chern classes of Heegner divisors*, LNM 1780.
* R. E. Borcherds, *Automorphic forms with singularities on Grassmannians*, Invent. Math. 132 (1998)
  — Theorem 13.3 is the lift itself, cited as `[Bo2]` by Bruinier and as `[B]` above.
