# Design: VERIFICATION.m for genus-2 quotient models + lowgenus ports

Date: 2026-06-22
Branch: main

## Goal

Bring two features from `lowgenus` to `main`, grab the polymake solution files the
relevant models need, and add a standalone `VERIFICATION.m` that reconstructs the
genus-2 models in two paper tables and checks each is isomorphic to the published
`y^2 + h(x)y = f(x)`.

## Scope (confirmed)

- Port: **only** `Targets` (target-restricted `EquationsOfCovers`) and the **normalize**
  (libnormaliz) polymake change. **Do not** bring the `LP_SIZE_CUTOFF` /
  `BORCHERDS_POINT_CUTOFF` safety caps that are bundled in the same lowgenus diff.
- Verify: **both** tables — the 15-row genus-2 non-bielliptic table (`modelsgen2`) and
  the 2-row genus-2 bielliptic table (`modelsbielliptic`, both rows on `X(6,17)`).

## Part 1 — Port `Targets` + normalize (surgical)

### `BorcherdsForms.m`
1. **normalize** — in `write_polymake_scriptfile`, replace
   `print $p->LATTICE_POINTS;` with:
   ```
   prefer_now "libnormaliz";
   print $p->LATTICE_POINTS_GENERATORS->[0];
   ```
   Identical output, ~18x faster; only affects newly generated scripts (cached
   solutions are unchanged).
2. **Targets** — add `Targets := {}` to the `BorcherdsForms` intrinsic signature and,
   when non-empty, drop every `rams` key whose `curves[k]`W` is not in `Targets`
   (the two hauptmoduls, keys -1/-2, are added later and always kept).
   `require #Keys(rams) gt 0`.
3. **Do NOT** add the `LP_SIZE_CUTOFF` / `BORCHERDS_POINT_CUTOFF` default blocks or
   the cutoff `error` checks inside `get_integer_prog_solutions`.

### `EquationsCovers.m`
- Replace the `EquationsOfCovers(Xstar, curves : Prec := 100)` body with the lowgenus
  version that takes `Targets := {}` and: passes `Targets` to `BorcherdsForms`;
  restricts `num_vals` (CM-point demand) to `target_keys`; and, when `Targets` is
  non-empty, filters `schofer_tab`K_idxs` down to the targets before the final solve.
  (The lowgenus version also drops the `[1/6]…[6/6]` vprintf timing lines; acceptable.)

No `.sig` regeneration is required by us — Magma regenerates signatures on attach.

## Part 2 — polymake solution files

Copy from `lowgenus` into `main`'s `polymake/` the `polymake_solution_*` (and matching
`polymake_script_*`) files for the 13 levels the table models need:

| level M | (D,N) |
|---|---|
| 348 | (6,29) | 364 | (14,13) | 340 | (10,17) | 812 | (14,29) | 444 | (6,37) |
| 492 | (6,41) | 516 | (6,43) | 1060 | (10,53) | 708 | (6,59) | 732 | (6,61) |
| 212 | (106,1) | 236 | (118,1) | 204 | (6,17) bielliptic |

Only copy files **not already present** on main (348/444/204 partially exist; the
lattice-point enumeration is deterministic so existing files are content-identical —
do not clobber them). This lets `VERIFICATION.m` run without polymake/libnormaliz.

## Part 3 — `VERIFICATION.m`

Standalone script (`magma VERIFICATION.m`), modeled on lowgenus
`eqns_genus2_bielliptic.m` but **verification-only** (no logging/resume/cost machinery).

**Data:** inline list of rows, each `<D, N, gens, h, f>` where `gens` is the AL-subscript
set (e.g. `{3,29}`), and `h,f` are the table polynomials. Target curve to match:
`HyperellipticCurve(f, h)` (`y^2 + h*y = f`). 15 + 2 = 17 rows.

**Per row:**
1. Build `W := AllALsFromGens(gens, D*N)`.
2. Find the star curve `Xstar` for `(D,N)` in `curves := GetHyperellipticCandidates()`.
3. Decide immediate vs deep: `W in { curves[i]`W : i in Xstar`CoveredBy }`.
   - **Immediate:** `crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := {W});`
     then pick the cover with `curves[k]`W eq W`.
   - **Deep** (the two `(6,17)` rows): `all_eqns, _ := AllEquationsAboveCovers(D, N, curves);`
     find key `k` with `curves[k]`{D,N,W}` matching, take
     `all_eqns[k][Representative(Keys(all_eqns[k]))]`.
4. `Ctab := HyperellipticCurve(f, h)`. Assert `IsIsomorphic(C, Ctab)` (over Q). The
   table uses different hauptmoduls per row, so only the isomorphism class is expected
   to match — not the literal polynomials.
5. Print `MATCH` / `MISMATCH` (and on mismatch, whether `G2Invariants` agree).

**Output:** one line per row plus a final summary `N/17 matched`. Wrap each row in
try/catch so one failure doesn't abort the sweep; exit nonzero if any row is not MATCH.

`(D,N)` → level grouping is incidental; rows are processed independently. Caching means
each star curve's Borcherds work runs once per `(D,N)` even across rows.

## Out of scope / non-goals

- Not porting `EnoughCMPointsForTargets` (lowgenus-only pre-check; verification just
  tries and catches).
- No changes to the pipeline, analysis, or ratpts drivers.
- No git commit (user commits manually).

## Testing / success criteria

`magma VERIFICATION.m` reports `MATCH` for all 17 rows using only the copied cached
polymake solutions (no polymake invocation needed).
