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
| 22 | `SpecialFiberNotHyperellipticD22` | **general CM-value** (rational+quadratic CM disc → s table); intermediate/star | N=p prime; **p≤43** (12 rational + 12 quadratic discs, added 2026-06 via the now-fixed `ValuesAtCMPoints`) | proves 33 (0 contradictions) — but **0 NEW**: all 23 undetermined D=22 curves are composite-N (Phase B), unreachable by this N=p method. base case still TODO |

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
   p, lift/push to the quotient. No Heun shortcut in general (only D=6/D=10 level-1 had one). NB: the
   quadratic-CM-disc enumeration crash that would otherwise block this (and limit coverage to rational
   discs) is now **fixed** — see §6; re-running `ValuesAtCMPoints` now yields the quadratic-disc values.
3. `has_compat` (shared), with `case2` as in §2.

Suggested Phase B order: start with the discriminants with the most reachable curves (D=6: 138,
D=10: 57) and the smallest M, reusing the level-1 lift/push machinery one level up. The genus-0
quotient models are the new ingredient — investigate `EquationsCovers.m` /
`AllEquationsAboveCovers` at level M > 1 (it produced the level-1 cover conics `f_q` we used).

**Phase B status (2026-06): the machinery is built, validated, and demonstrated.**
- Generic driver `SpecialFiberNotHyperellipticCM(D,N,W)` + per-`(D,M)` `CM_BASE_DATA` registry
  (`special_fiber_cm.m`, committed 60bac2b), dispatched in `FilterBySpecialFiber` after the M=1 tests.
  Reuses every level-1 helper; validated to reproduce `SpecialFiberNotHyperellipticD22` exactly
  (400 curves, 0 mismatches).
- **First base `<6,5>` registered** (staged): X_0(6,5)* — 5 rational discs + 6 genus-0 conic covers.
  D=6 M=5 prototype: **reached=391, proven=3, contradictions=0, 3 NEW** undetermined curves resolved.
- Per-base data generator: replicate `EquationsOfCovers`'s body but filter `Xstar`CoveredBy` to
  genus-0 covers (drops the data-starved higher-genus cover that forces `num_vals` too high) and set
  `Include := {}` (the Hauptmodul-divisor Include triggers a giant `bd=include_bd*2` CM search — ~80s+).
  Script lives at `/tmp/gen_base.m` in the working session; X_0(6,5) generated in 89s.
- TWO Phase-B cost limits to keep in mind: (i) **coverage** — bases have few CM discs (X_0(6,5)* only
  6), so the prime ceiling is low (add discs via higher MaxNum / the quadratic disc); (ii) **the Embed
  hang** (now fixed by the hybrid `FindLambdas`, §6) — without it, generation stalled 13–38 min.

## 6. Known issues / TODOs (carry these forward)

- **Quadratic-CM-points bug — FIXED (2026-06; staged on `SchoferFormula.m`, uncommitted).** This was
  the coverage limiter: `ValuesAtCMPoints` blew up (80 GB / "Magma internal error") on large quadratic
  CM discs, capping the general CM-value method to rational discs (D=10 ~p≤61 patchily, D=22 only p≤17).
  Root cause was NOT the min-poly recognition (`find_hauptmodul_signs_quadratic` was fine) but
  `ElementsOfNorm`/`FindLambdas`: they found the CM element by brute-forcing `CartesianPower([-bound..
  bound], n)` with `bound` starting at `d/2` and **doubling** — `O(d^3/d^4)` on the **indefinite**
  ternary norm-form (a solution has coordinates only `~sqrt(d)`). Fixed by replacing the box with
  Magma's `Embed` (an exact optimal embedding of the disc-(-d) CM order into the quaternion order `O`)
  in all four functions `FindLambdas`, `FindLambda`, `ElementsOfNorm`, `ElementOfNorm`. Subtleties that
  bit: the indefinite form rules out `ShortVectors`; the CM order is **disc -d** (`Z[sqrt(-d/4)]` for
  `d = 0 mod 4`, the maximal order for `d = 3 mod 4`), NOT disc `-4d` (else `elt = lam/2` leaves `O` and
  the Smith-form optimality test rejects it); the coordinate vector is `v = `coords of `lam` directly
  (`v.Q.v = 2d`); a non-embeddable disc (e.g. 1000: conductor 5, 5 ramified in B) now **errors cleanly**
  instead of looping forever. Sub-second even for disc 10^6. Validated: `Kappa0` + `InternalBorcherds`
  regression tests pass (identical Schofer values to the old box); `ElementsOfNorm([40,235,10^6])` correct
  in <1 s. **UPDATE — pure-Embed HANGS; now HYBRID (2026-06, staged on `SchoferFormula.m`).** During
  Phase-B prototyping on X_0(6,5)* the Embed-only `FindLambdas` was found to be non-deterministic and to
  HANG: `ShimuraCurveLattice` returns a RANDOM order representative, and for some, Magma's `Embed` enters
  a non-terminating CPU-bound search (100% CPU, flat ~64 MB RSS — a hang, not a memory blowup; e.g. disc
  -24 on X_0(6,5)*). That, not the Schofer/Borcherds work, made Phase-B generation stall 13–38 min.
  `FindLambdas` is now HYBRID: PRIMARY is a BOUNDED box (start 16, double to a hard cap 128, so no 80 GB
  blowup — the old box exploded only because it started at max_d/2 and doubled unbounded), with Embed as
  FALLBACK for discs whose optimal vector exceeds the cap. This routes the small hang-prone discs through
  the deterministic box. Re-validated: Kappa0/Schofer(D6,10,26,142)/CMPoints/X0_6_11/X0_26_1 all pass;
  D=22 `ValuesAtCMPoints` regenerated → 24 discs, 0 mismatches vs committed `d22cm_s`; X_0(6,5) now 89 s.
  (A cleaner long-term fix is our OWN bounded indefinite-ternary `Embed`.)
  **PAYOFF — DEMONSTRATED end-to-end for D=22 (2026-06, staged on `special_fiber_cm.m`):**
  `ValuesAtCMPoints(D=22 star, curves : MaxNum:=24)` now runs (209 s, no crash; polymake LPs cached so
  sandboxed was fine) and yields **12 quadratic discs** alongside the 12 rational ones. `d22cm_discs`/
  `d22cm_s`/`d22cm_star_coords` extended to use them: a quadratic s-value is stored as its degree-2 min
  poly and reduced mod p by `cm_red_minpoly` (clear denominators to a binary form, factor over `F_{p^2}`
  = the Frobenius-conjugate ss pair, with a root at infinity when p|denom). This lifts D=22 from p≤17 to
  **p≤43** (SpecialFiberD22 test: reached=90, proven=33, **contradictions=0**). CAVEAT: it gives **0 new**
  determinations — all 23 undetermined D=22 curves are composite-N (Phase B), which the N=p method can't
  touch. So the real value is the **validated end-to-end CM-value pipeline** (the exact machinery Phase B
  needs). To push D=22 further: larger `MaxNum` adds more discs (gap grows ~8 by p=97); D=10 could be
  extended past p≤61 the same way (less pressing — it has the Heun closed form). The generic recipe to
  reduce a quadratic value mod p in any special-fiber test: factor its min poly over `F_{p^2}`.
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
