# Special-fiber non-hyperellipticity — handoff

Status snapshot for the special-fiber (reduction-mod-p) non-hyperellipticity test and its
remaining frontier. Written for another Claude instance picking up the work.

## 1. What this is

Goal of the project: classify which Atkin–Lehner quotients `X_0(D,N)/W` of Shimura curves
(`D` = quaternion discriminant, `N` = level) are **geometrically hyperelliptic**. The repo
runs a pipeline of filters; `FilterBySpecialFiber` is one of them, marking curves **not**
(geometrically) hyperelliptic. The whole project is committed to *geometric* hyperellipticity
(see `ShimuraQuotients.m` line 4); the special-fiber test proves geometric non-hyperellipticity
because the hyperelliptic involution of a genus ≥ 2 curve is unique, hence defined over Q, hence
reduces to an F_p-rational Frobenius-compatible Möbius involution — and the test rules such out.

## 2. The method (Furumoto–Hasegawa Section 5, generalized)

For `X_0(D,N)/W` and a prime `p | N` (with `p ∤ D`): the special fiber at `p` of `X_0(D,N)`
degenerates into copies of `X_0(D, N/p)` crossing at the **supersingular points** of `X_0(D,N/p)`
mod p. After the quotient by `W`, the genus-0 special-fiber **component** is `X_0(D, N/p)/W''`,
where `W''` is the p-coprime image of `W` (the "component group"). The test needs this component
to be **genus 0** (a rational curve, so the supersingular points are P^1 coordinates).

- **case 1** (no `w_m ∈ W` with `p | m`): two copies of the component crossing at the ss points;
  ALL supersingular points are constrained.
- **case 2** (some `w_m ∈ W` with `p | m`): a single self-intersecting copy; only the
  non-F_p-rational ss points are constrained. Computed as `case2 := exists{w in W | w mod p eq 0}`.

`has_compat(qss, p, F, case2)` (in `special_fiber_modular.m` and `special_fiber_cm.m`) then asks:
is there a Möbius involution on the genus-0 component (over F_{p^2}) that acts as Frobenius on the
constrained ss points and is F_p-rational? Returns `"no"` → NOT hyperelliptic; `"yes"`/`"underdet"`
→ inconclusive. This logic is shared across all discriminants.

**The hard part for each case is COORDINATES**: getting the supersingular points of `X_0(D,N/p)/W''`
as explicit P^1 coordinates. Everything below is about that.

## 3. What's implemented (all wired into `FilterBySpecialFiber`)

| D  | intrinsic | component coordinates | coverage | status |
|----|-----------|----------------------|----------|--------|
| 1  | `SpecialFiberNotHyperelliptic` | classical `ss_points` via `SupersingularPolynomial` on genus-0 `X_0(M)`; `push_group` to `X_0(M)/W''` | composite N; `\|W''\| ≤ 4` (Klein-four) | done |
| 6  | `SpecialFiberNotHyperellipticD6` | (2,4,6) **hypergeometric** `d6_hypg` (star is a triangle group) + cone vertices | N=p prime, intermediate quotients | done |
| 10 | `SpecialFiberNotHyperellipticD10` | (2,2,2,3) **Heun** closed form `d10cm_heun_generic` + cone points; base/intermediate/star lifts | N=p prime; full prime coverage (Heun) | done; **proved new: X_0(10,71)/⟨w_10,w_71⟩** |
| 22 | `SpecialFiberNotHyperellipticD22` | **general CM-value** (rational CM disc → s table); intermediate/star | N=p prime; only p≤17 (no Heun, see §6) | feasibility only, 0 new; base case TODO |

Files: `special_fiber_modular.m` (D=1 classical + D=6 hypergeometric + the Brandt layer
`SupersingularALData` + `FilterBySpecialFiber` dispatch), `special_fiber_cm.m` (D=10 Heun + D=22
general CM-value + shared `cm_inter_conic_P1`). Tests: `tests/SpecialFiberD{1,10,22}.m`. Cost proxy
for the parallel filter: `Utils.m` (`CurveCostProxy`, the `FilterBySpecialFiber` stage).

### The supersingular-point sources (three ways, increasing specificity)
1. **Brandt module** (`SupersingularALData(D,p)` in `special_fiber_modular.m`): ss points of
   `X_0(D,1)` mod p = left-ideal classes of the definite quaternion algebra of disc `D·p`. Gives
   the COUNT, Frobenius (= AL `w_p`), and AL action — but NOT P^1 coordinates. The count equals
   `genus(X_0(D,p)/full-AL) + 1`; the production code now uses that **genus formula** instead of a
   Brandt module call (much faster — see the optimization that cut the D=10 test 310s→1.7s).
2. **CM-value coordinates** (general): the Hauptmodul value at a CM point, from the repo's
   Borcherds/Schofer pipeline (`ValuesAtCMPoints` in `SchoferFormula.m`), reduced mod p. Works for
   any D the pipeline handles, but COVERAGE is limited to the rational CM discriminants (quadratic
   discs crash — see §6).
3. **Heun/hypergeometric closed form** (D=6, D=10 only): the supersingular polynomial as the
   mod-p reduction of the Picard–Fuchs ODE. Full coverage, no CM-coverage limit. Only exists when
   the star quotient is a triangle group (D=6) or single-accessory Heun (D=10).

### The D=10 Heun derivation (the most reusable closed-form idea)
The (2,2,2,3) star has singular points `s = 0, 1, 2/27, ∞`. The **generic** supersingular points
(away from the cone points) are the roots of the polynomial solution of the Heun ODE
`s(s-1)(s-2/27) y'' + Q y' + (αβ s − q) y = 0`. Derived facts (all in `special_fiber_cm.m`):
- exponent at each order-2 cone point = **1/2 if ordinary, 3/2 if supersingular** (Kronecker rule);
- `αβ` fixed by the degree N (= `star_n − #supersingular cone points`) via the Fuchs relation;
- the accessory `q` is an **eigenvalue** of the Heun operator `A: f ↦ P f'' + Q f' + αβ s f`
  (square on degree ≤ N once αβ is fixed), selected by one known generic root from the CM data.
The **cone points** (small-CM elliptic points) are added separately, exactly as `d6_ss_star` adds
the elliptic vertices to the hypergeometric roots.

## 4. Running things

- Tests: `magma "target:=SpecialFiberD10" "exitsignal:=true" run_tests.m` (substring match; e.g.
  `target:=SpecialFiber` runs D1/D6/D10/D22). Each asserts 0 contradictions + effectiveness.
- The pipeline: `run_pipeline.sh` (parallel) / `workingcode.m` `GetHyperellipticCandidates`
  (sequential). `FilterBySpecialFiber` runs first in the all-quotients block (it is the cheapest
  filter, so its determinations prune the rest). Inputs/outputs chain via `run_parallel_filter.sh`.
- polymake: the Borcherds-form basis (D > 1) shells out to polymake (`BorcherdsForms.m`). It was
  broken and is now fixed (perl + tap-trust + `brew reinstall polymake` 4.15). **Magma runs that
  hit an UNCACHED polymake LP must run outside the sandbox** (`dangerouslyDisableSandbox`), else
  `pwd: Operation not permitted`. Cached LPs are fine sandboxed.

## 5. The frontier: composite-N (Phase B)

Scoping over the undetermined g≥3 curves (`data/curves_after_UpdateCurves8.dat`):
- **827 undetermined**; **397 reachable** by composite-N special fiber (a genus-0 component at some
  `p|N` with `N/p > 1`) — ~48%.
- Split: **23** have the component *base* `X_0(D,M)` itself genus 0 (all D=1 — Phase A territory);
  **374** have only the *quotient* `X_0(D,M)/W''` genus 0 (16 D=1, **358 D>1**).
- By D: D=6→138, D=10→57, D=14→37, D=15→22, D=22→16, D=1→39, plus ~20 more discriminants.

**Phase A (done):** the 23 + the D=1 group-quotient cases. `SpecialFiberNotHyperelliptic` now
loops over `p|N` (composite N) and handles `|W''| ≤ 4` via `push_group`. Proved **17 new D=1
determinations**, 0 contradictions.

**Phase B (open, the bulk — 374 curves):** genus-0 *quotients* of positive-genus bases. For each
`X_0(D,N)/W` with a genus-0 component `X_0(D, M)/W''` (M = N/p > 1), need:
1. supersingular points of `X_0(D,M)` mod p = ideal classes of the **Eichler order of level M** in
   the definite algebra of disc `D·p`: `BrandtModule(D*p, M)` (note the level argument). Count =
   `genus(X_0(D, pM)/W) + 1`-style node count; use a genus formula for the coverage check.
2. **coordinates** on the genus-0 quotient `X_0(D,M)/W''`. This is the real work. The level-1
   coordinate methods generalize: CM points of `X_0(D,M)` have Hauptmodul values (repo's
   Borcherds/Schofer pipeline already computes level-N>1 models — the dataset has them), reduce mod
   p, lift/push to the quotient. No Heun shortcut in general (only D=6/D=10 level-1 had one).
3. `has_compat` (shared), with `case2` as in §2.

Suggested Phase B order: start with the discriminants with the most reachable curves (D=6: 138,
D=10: 57) and the smallest M, reusing the level-1 lift/push machinery one level up. The genus-0
quotient models are the new ingredient — investigate `EquationsCovers.m` /
`AllEquationsAboveCovers` at level M > 1 (it produced the level-1 cover conics `f_q` we used).

## 6. Known issues / TODOs (carry these forward)

- **Quadratic-CM-points bug** (the coverage limiter). `ValuesAtCMPoints` with quadratic CM discs
  (class number ≥ 2, value = a degree-2 min poly) hits a hard **Magma internal error** (uncatchable)
  and an 80 GB blowup on large ones. Rational discs only → coverage caps (D=10 ~p≤61 patchily,
  D=22 only p≤17). The quadratic points SHOULD work (some repo models were built from them), so this
  is a **bug to fix**; fixing it unblocks full coverage for every D without a Heun ODE. To reduce a
  quadratic value mod p: factor its min poly over F_{p^2} (roots = the Frobenius-conjugate ss pair).
- **Genus-0-quotient criterion** (already exploited in Phase A's design, central to Phase B): the
  test applies whenever `GenusShimuraCurveQuotient(D, M, W'') = 0`, NOT only when the base
  `X_0(D,M)` is genus 0. Enumerate on the quotient genus.
- **D=22 base case** (`W''={1}`): TODO. The Klein-four base lift works (`f_2 f_11/f_22 = (11(s-1))^2`,
  base conic `256 z^2 + 1331 w^2 + 1936 t^2 = 0` with `z=y_2/s, w=y_11/y_2`) but the special-point
  casework wasn't finished. The D=10 base (`d10cm_base_P1`, conic `z^2+25w^2+2t^2=0`) is the model.
- **D=1 larger group quotients** (`|W''| > 4`): `push_group` handles up to Klein-four; bigger AL
  groups (order 8) would need more push-through-generator iterations.
- **Other genus-0 Shimura discriminants**: only D=6, 10, 22 have a genus-0 base. All other D appear
  only through the genus-0-quotient criterion (Phase B).

## 7. Pointers

- Memory notes (agent memory, `.claude/`, may not travel): `special-fiber-genus0-quotient-criterion`,
  `quadratic-cm-points-bug` — both summarized above.
- Validation oracle for D=1: `verify_d1_hyperelliptic.m` (Petri/canonical-model, decides geometric
  hyperellipticity independently — use to cross-check D=1 results).
- The new D=10 result and the full method validation live in `tests/SpecialFiberD10.m`.
