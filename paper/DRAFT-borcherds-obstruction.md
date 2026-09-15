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

## 5. ⚠ The open piece, and where it will bite

The predictor needs `dim S_{3/2}(ρ_L^*)` in closed form. Two warnings:

* The relevant dimension formula is the Riemann–Roch / trace-formula statement (Borcherds' §4;
  Bruinier, *Borcherds products on O(2,l)*, LNM 1780, §2). It is **not** in §4 of "Reflection groups
  of Lorentzian lattices" — checked, that section is about Lorentzian reflection groups.
* **Weight 3/2 is the singular case.** Riemann–Roch computes `dim M_k − dim S_{2−k}`; at `k = 3/2`
  the partner is `M_{1/2}`, spanned by unary theta series (Serre–Stark). The formula does not
  separate the two, so the `M_{1/2}` contribution must be computed explicitly. This is where an
  unchecked formula returns a plausible wrong number.

⚠ CONJECTURAL, NOT ESTABLISHED: by the Eichler/Shimizu and Shimura correspondences one expects
`S_{3/2}(ρ_L^*)` to be related to weight-2 forms on the Shimura curve (in the GKZ case the weight-3/2
forms are the Kohnen plus-space partners of `S_2(Γ_0(N))`). If that holds there may be a route to the
dimension that avoids the singular-weight formula entirely. **Nothing here has been checked.**

**Calibrate before trusting.** Reproduce the KNOWN deficits first — `38_5` (1), `146_1` and `194_1`
(0, both build), `58_13` (2), `26_31` (3) — before believing any new number. The repo's own record of
what happens otherwise is `HANDOFF.md`, and this session alone retracted three readings that were
correct arithmetic about the wrong object.

## Sources

* R. E. Borcherds, *The Gross–Kohnen–Zagier theorem in higher dimensions*, Duke Math. J. 97 (1999)
  219–233 (+ erratum, Duke 105 (2000) 183–184). `https://math.berkeley.edu/~reb/papers/gkz/gkz.pdf`
* J. H. Bruinier, *On the converse theorem for Borcherds products*, arXiv:1210.4821.
* J. H. Bruinier, *Borcherds products on O(2,l) and Chern classes of Heegner divisors*, LNM 1780.
* R. E. Borcherds, *Automorphic forms with singularities on Grassmannians*, Invent. Math. 132 (1998)
  — Theorem 13.3 is the lift itself, cited as `[Bo2]` by Bruinier and as `[B]` above.
