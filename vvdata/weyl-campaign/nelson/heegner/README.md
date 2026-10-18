# Heights of the branch points on an elliptic quotient

Scripts run from a checkout of `main` carrying `ExactDoubleCover.m` (PR #123), 2026-10-18.

* `h2.m`: X_0^35(1) and X_0^51(1). The canonical heights of the branch points of X -> E = X/w on
  E are exactly rational multiples of L(E,1) L'(E^d,1) sqrt|d| / Area(E) (Area = covolume of the
  period lattice): 2 and 1/2 at D = 35 (d = -7; the -28 point is -2 times the -7 point, a_2 = 0),
  3/2 and 0 at D = 51 (d = -3; one point is torsion). Thirty digits.
* `scan.m`, `scan.log`: the same over every stored genus-3 hyperelliptic quotient whose completed
  octic is even and whose elliptic quotient has a small rational point: 70 cases. The ratio is
  rational in every case where it is defined (L'(E^d,1) ≠ 0), ranging from 1/4320 to 18, with the
  square of the unit index visible at d = -3 (9/32, 9/14, 9/22, 9/26, 18 = 2·9). Branch points over
  REAL quadratic fields and rational branch points occur (then the comparison quantity must be
  built from L'(E,1) and the twist differently; not done). The rational factor is not yet derived:
  it should be the Gross-Zagier/Zhang constant for the parametrisation X_0^D(N) -> X/W -> E, with
  the unit index, the degree 2|W|, the conductor of the CM order, and the quaternionic Petersson
  norm. Deriving it and checking it on these 70 cases is what turns the heights into an a priori
  bound, which is what the exact route to the level-5460 cases (210_13) needs.

## The constant, pinned (later on 2026-10-18)

Cai–Shu–Tian, "Explicit Gross–Zagier and Waldspurger formulae", Algebra & Number Theory 8 (2014),
arXiv:1408.1733, Theorem 1.5 (special case 2, F = Q, c = 1): for f : X_0^D(N) → E of degree deg f and
a CM point P of fundamental discriminant d,
  L'(1, E, χ) = 2^{-#Σ_D} · 8π² (φ,φ)_{Γ_0(N_E)} / (u² √|d|) · ⟨P⁰_χ(f), P⁰_{χ̄}(f)⟩_K / deg f,
P⁰_χ(f) = Σ_{t ∈ Pic(O_d)} χ(t) f(P)^{σ_t}, Σ_D the primes dividing both N_E and d, u = |O_d^×|/2, the
height ⟨,⟩_K over K (twice the absolute height), and 4π²(φ,φ) = deg(π_opt)·Area(E_opt)/c² with
π_opt the optimal parametrisation by X_0(N_E) and c the Manin constant. The test vector condition
(f must be χ_v(ϖ_v)-eigen under the Atkin–Lehner w_v at each v ∈ Σ_D) is where our f = X → X/W → E
enters: f∘w_v = ε_v f up to a rational point, ε_v = a_v(E) for v | D and −a_v(E) for v | N
(`conv.m`: this convention is the only one for which E's signs are +1 on W', in all 48 cases), and
Shimura reciprocity (σ_{𝔭_v} acts as w_v on the CM points of disc d) makes P⁰_χ(f) torsion unless
χ(𝔭_v) = ε_v for every v ∈ Σ_D. Hence, with deg f = 2|W| and genus characters χ = (d1, d2), for
which L(s, E, χ) = L(E^{d1}, s) L(E^{d2}, s):

  ĥ_abs(f(P)) = (1/2)(1/h²) Σ_{χ admissible} u² √|d| · 2|W| · 2^{#Σ_D} · L'(1,E,χ) / (2 deg π_opt Area(E_opt)).

`cst3.m` (log `cst3.log`): with the TRUE CM discriminants of the branch points (the fixed-point data,
not the field of definition — they differ as soon as h > 1), **82 of 92 branch points on the stored
models match to six digits**, class numbers 1, 2, 4, 8. The 10 mismatches are all on E of rank 1,
where the formula concerns P − ξ (Hodge class) and the trivial-character component is shifted by a
rational point; the character components are unaffected. Earlier attempts (`cst.m`, `cst2.m`,
`genus.m`) are kept: `cst2.m` shows the admissibility rule, `genus.m` what goes wrong without it.

`pipeline.m` is the route end to end (E by signs, branch discriminants and fields, predicted
heights, Mordell–Weil of E and its twist over F, lattice points of the predicted height,
`ExactDoubleCover` against the trace formula, comparison with the stored model). Control
X_0^10(11)/{1,5} (`pipe_10_11.log`): the heights exclude 110a1 and select 110a2 within the isogeny
class, each branch point is exactly twice the generator of E(F), one configuration survives, and
its model is Q-isomorphic to the stored one. Paper copy: `CaiShuTian-explicit-GZ-arXiv1408.1733.pdf`
in the scratchpad (not committed).

## Intrinsics (PR #124, branch heegner-heights on #123)

`HeegnerHeights.m` on `main`'s side carries the route as intrinsics; `tools/heegner_model.m` reruns a
search; `tests/HeegnerHeights.m` holds the calibration, the two controls end to end, and the two models
as regression values. `repro_1326.log` and `repro_210.log` are the reproductions of both models through
the intrinsics (one curve up to Q-isomorphism each; every other curve of the isogeny class excluded by
the lattice test). Only the class mod 2 of the origin-dependent component matters (translation by a
rational point of E gives an isomorphic cover), so the enumeration is complete.
