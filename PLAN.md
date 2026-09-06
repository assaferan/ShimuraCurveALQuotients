# Working plan — what to do next

Assembled 2026-09-02 from the branch audit, `HANDOFF.md` on `main`, and the GATE3B measurement on
`10_61`. `HANDOFF.md` remains the record of *what happened*; this file is the record of *what to
do*. When they disagree about state, `HANDOFF.md` wins.

Five tracks. One is the main line; the rest run in parallel and **none of them blocks it**.

## Picking this up cold

> **⇒ Spend your effort on WHICH OBJECT a claim is about, not on whether the computation is right.**
> The 2026-09-04 session's most expensive errors all had correct arithmetic about the wrong object.
> Reproduce a KNOWN value before trusting a new one; draft an edit rather than applying it.
> Full account: `HANDOFF.md`, "READ THIS FIRST".

As of **2026-09-04** everything is committed and pushed, both branches and lovelace are in sync,
and the housekeeping list is empty. **The ordering below CHANGED on 2026-09-04** — the per-coset
`tau` fix landed, and it moved the frontier. In decreasing order of value:

1. **Run the full model pipeline on `34_11`** (NEW, and the best-value thing here). It was
   blocked at class-constancy; with the `tau` fix it now clears **every** gate in
   `M0MultiplierExact` and returns multipliers confirmed 9/9 against the closed form. Nothing is
   known to stand between it and models. **The project has produced models from exactly ONE base
   ever** (`58_5`, 4 models) — a second would be the largest single change to the state of play in
   months. Bounded, needs no theory. Cost: `AllEquationsAboveCovers` on a base this size ran 18–37
   min for its siblings, and its `M0MultiplierExact` alone took 2.4 h, so budget hours not minutes.
2. **Re-test `10_61` and `14_43`** (cheap, and directly implied by the fix). Both died at the
   gate-3 slash-constant check — *exactly* the failure the per-coset `tau` addresses, and `10_61`'s
   two poisoned cosets (`wi = 2`, `wi = 1221`) sit at precisely the two `Im(z)` extremes the
   fallback now catches. They may simply run. If they do, the "two runnable candidates of 122"
   count and the 49-base obstructed class both need revisiting.
3. **The `A_m` theorem** (MAIN LINE). Still the biggest prize — 49 bases — but ⚠ **its "concrete
   next step" as written was FALSIFIED on 2026-09-04**: the missing cusp classes are NOT a dump
   bug, so there is no quick implementation task waiting. See that section; it is a search problem
   again.

⚠ **Also newly in doubt: "gate 4 is a genuine violation."** `34_11` now passes class-constancy
with the guard untouched, which suggests that 43% deviation was the same evaluation-point defect
as gate 3, one gate later. Not yet confirmed — see REPAIR.

Worktrees: `main` (here, `.`) and `worktrees/campaign` (`m0-theta-campaign`, research data);
`tier1-models` and `whbasis-speedup` are retired, 10 `archive/*` tags on origin. Probes live at
`vvdata/weyl-campaign/tau-precision/` on the campaign branch.
⚠ Before trusting either branch's code, run the shared-path check in `CLAUDE.md` — two capability
gaps hid there for days.

## State of play

     4   verified models, from one base (58_5, partial set by construction)
     1   more base now CLEARING every M0MultiplierExact gate (34_11) -- pipeline not yet run
    49   Borcherds-obstructed bases (unchanged; A_m still gates these)
    81   of 122 still unclassified
   ---   gate-3 failures at 10_61: FIXED. Confirmed 2026-09-04/05 -- 10_61 and 14_43 both
         completed their a0 tables (27 and 22 fallback points of 64), so the two-point check
         never fired. Do NOT quote the old "924 failures, zero near-misses" as current.
    34   Guo-Yang CM-value tables now verified as offline tests, 254 checks, 0 failures
    30   pages of paper; the theory arc is closed

Two arcs running at different speeds. The theory arc is essentially closed — the multiplicity
formula end to end, `kappa_mu(0)` derived, N>1 subsuming N=1, the support rule correctly reframed
as a gauge. The model arc is blocked, and not by engineering.

Across the record, most proposed levers were refuted once measured: the 336x hoist worth 2%, the
odd-D eta-quotient explosion that wasn't, three successive symbol laws, and the scalar `a_E` at
1 of 13. That is discipline, not failure.

⚠ **But the conclusion this file used to draw from it — "the engineering levers are spent and the
binding constraint is mathematical" — took a real hit on 2026-09-04.** The per-coset `tau` fix is
an *engineering* fix, and it unblocked `34_11`, a base this file had recorded as failing on
"real mathematics" (a GENUINE 43% class-constancy violation, explicitly "NOT another tolerance").
That diagnosis now looks wrong. The honest version: **the binding constraint is mathematical for
the 49 Borcherds-obstructed bases, where the deficit is invariant and `A_m` is genuinely missing.
It was NOT binding for the gate-3/gate-4 bases — those were numerical all along.** Do not
generalise "it must be mathematics" from one class to the other; that framing cost `34_11` a
month.

---

## SPEEDING UP THE EXPENSIVE BASES — measured 2026-09-05

The five odd-`D` level-1 bases (`93_1 95_1 159_1 111_1 119_1`) are the largest block of
unreproduced Guo-Yang bases and need 8+ hours each. Three constant-factor wins shipped, and two
algorithmic hypotheses were tested and REFUTED. Recording both so neither is re-attempted.

### Shipped: 41% off `BorcherdsForms` (137.3 s -> 81.1 s on a full `genmodels` run of `51_1`)

    c04f938  &+ single-pass          the fold copied the whole accumulator per addition and
                                     re-reduced -- O(k^2).            137.3 -> 107.6
    1d8b2fc  reduce stops sorting    it went through Exponents (which sorts) to find zeros;
                                     removal order is irrelevant.     107.6 -> 103.4
    33379fa  power cache             qExpansionAtoo recomputed nor_eta_ds[i]^r[i] per quotient
                                     -- 3451296 series '^' calls.     103.4 ->  81.1

⚠ The linear algebra was NEVER the bottleneck: `EchelonForm` is 6 s of 222 s. It was all
data-structure overhead in the eta-quotient layer. Every step verified by byte-identical
regeneration of `51_1` and `14_3` plus `ModelChecks` 8767/0.

### ⚠ REFUTED (1): the basis pool is NOT redundant

`BFPROGRESS=1` now reports `BFPOOL pole_order=... pool=... rank=...`. Measured at `51_1`,
`pool ~ rank` at every call (`159/154`, `37/37`, `149/149`, `173/173`, `253/253`, `312/307`).
There is no repeat of the 66x `WeaklyHolomorphicBasis` win (a rank-258 space echelonised as
12784 rows) hiding here. Spanning a rank-253 space genuinely needs 253 q-expansions.

### ⚠ REFUTED (2): "find a SHORT eta-quotient combination instead of a full basis" — base-dependent

The idea, from the observation that Guo-Yang found their forms almost by hand: skip the complete
basis up to `pole_order` and search directly for a short combination with the required principal
part. Measured how many eta-quotient terms our produced forms ACTUALLY use:

    X_0^6(1)   : 4 to 7 terms      <- hand-scale, matches the literature's small bases
    X_0^51(1)  : 433 terms         <- not hand-scale

⇒ **It would work on the bases we already solve in seconds and not on the expensive ones.** The
"almost by hand" framing holds at small `D` (`6_1`, `10_1` -- Errthum's bases) and breaks by
`D = 51`.

⚠ ONE THING STILL OPEN, and it is the only way this idea survives: 433 is the length of OUR
representation. If the eta quotients used are linearly DEPENDENT, a shorter representation may
exist and our construction is simply not finding it. `FindMinimalEtaQuotient` exists for exactly
that question. Settle that before discarding the approach entirely — but do not assume a short
form exists just because one does at `6_1`.

### Where the remaining time goes

After the three fixes, re-profiled: `Constructor (sub)` 113 s / 1.23M calls (Magma-internal
allocation), `qExpansionAtoo` 53 s / 12531 calls, `&*` 13.8 s / 262860. These look like
allocation rather than anything with an obvious algorithmic fix, so the next constant-factor
increment is likely smaller and harder-won than these three.

### How the cost SCALES — measured 2026-09-05, and it corrects a claim above

`basis_of_weakly_holomorphic_forms` at `51_1`, timing its three parts against `pole_order`
(`scratchpad/poolcost.m`; pool grows linearly in `pole_order` because `full_basis` is the
t-EXPANDED pool):

| `pole_order` | pool | qexps | `EchelonForm` | `ech_etas` | total |
|---|---|---|---|---|---|
| 100 | 86 | 0.40 s | 0.01 s | 0.42 s | 0.83 s |
| 200 | 186 | 2.75 s | 0.26 s | 2.89 s | 5.90 s |
| 400 | 386 | 22.09 s | 5.49 s | 16.26 s | 43.84 s |
| 800 | 786 | 213.84 s | 112.94 s | 90.38 s | 417.16 s |

⚠ **"`EchelonForm` was only 6 s of 222 s" does NOT extrapolate.** That was measured at a small
pole order. `EchelonForm` has the STEEPEST growth of the three (~x20-26 per doubling, ~`PO^4.4`,
vs ~`PO^3` for qexps): 1% of the cost at `PO=100`, 27% at `PO=800`, and on this trend it overtakes
qexps around `PO ~ 1200-1600`. Any future profile must say at which `pole_order` it was taken.

### The t-ladder leaves the pool ALREADY TRIANGULAR — structure confirmed, exploitation REFUTED

Measured at `51_1`, `87_1`, `55_1`, `15_2` (`scratchpad/poolstruct.m`), with `k` = 32, 56, 37, 8
and `r` = 11, 6, 10, 49 — so this is not one base's accident:

* every pool element has a **distinct valuation** (e.g. 386 distinct over 386 forms);
* sorting rows by valuation puts **100% of rows in triangular position**, every consecutive pair
  strictly increasing.

This is structural, not luck: `WeaklyHolomorphicBasis` already returns an echelon `E`, and
multiplying a block by `t^j` shifts every valuation by exactly `-jk`, so the size-`k` blocks tile
the pole range without collision. So `EchelonForm` never SEARCHES for pivots — all its time is
elimination above pivots that are already on the diagonal.

⚠ **The obvious inference — "so skip the reduction" — is REFUTED by measurement.** At `51_1`,
`PO=800` (`scratchpad/echcost.m`):

    EchelonForm (full RREF)          : 112.910s
    sort rows into valuation order   :   0.030s
    pivot columns identical          : true (786 vs 786 pivots)
    max |num|,|den|: RREF 1 digits, sorted pool 33 digits

The RREF turns **33-digit** rational entries into **1-digit** ones. It is not wasted elimination;
it is buying a vastly better-conditioned basis. Feeding the raw sorted pool downstream would push
33-digit rationals through `ech_basis * mat` and the linear solve, RELOCATING the cost and most
likely increasing it. Do not re-propose this without first measuring the downstream cost on
33-digit input.

What survives as a lever: the pivot order is known in advance, so a solver specialised to a
triangular system with a fixed pivot sequence (or a fraction-free/modular one) could beat the
generic call — but the target is the coefficient growth, not the pivoting.

### The untested lever: q-expansions of the t-expanded pool

`qexps` is still the largest single term (213.84 s of 417 s at `PO=800`). We already bootstrap the
BASIS the way Guo-Yang do — `BorcherdsForms.m:600`, `full_basis := [t^r*f : f in init_basis]`, with
the basis computed only at `n0+k` — but line 607/609 then pays a from-scratch eta-quotient
expansion for EVERY pool element. `t^j*f` has the same number of eta-quotient terms as `f` (the
exponent vectors merely shift), so each costs full price, and there are `r` times as many.
The alternative is `qexp(t^j*f) = qexp(t)^j * qexp(f)`: expand only `init_basis`, `basis_n0` and
`t`, then multiply series. ⚠ NOT free — the products need higher ABSOLUTE precision on the factors
(`1 + jk` and `1 + v` respectively), so measure before believing it. NOT YET TESTED.

---

## PRODUCING MORE MODELS RELIABLY — the plan, 2026-09-05

*(adopted after the session that recovered `39_2` and `14_3`. Priority set by assaferan: make the
REAL REPAIR first, not the workaround.)*

### Where we are

Coverage **34 of 43**. Both recoveries came from the same place, and it was not mathematics: the
**coprime-to-level CM filter** was starving the CM pool (`39_2`: 3 points against demand 19; with
the filter off, 24). Three blocked bases examined, three whose recorded diagnosis had gone STALE.

### The problem to solve

`CMNONCOPRIME=1` is the lever that works and it has **no theoretical guarantee**: the
`p | gcd(d,N)` local factor has NO live implementation (`kappaminuszero`, `SchoferFormula.m:569`,
is dead code and does not even cover that case — it loops over `PrimeDivisors(N div GCD(d,N))`,
the coprime part), and at `26_3` two non-coprime discriminants give provably wrong values. So its
output has only ever been trusted where Guo-Yang publishes an equation to check against.
**That does not scale**: 47 of our models are beyond Guo-Yang's list, and any new base has no
published answer.

⇒ **THE REAL REPAIR IS THE `p | gcd(d,N)` SCHOFER VALUE.** Fix it and the filter can be removed,
which turns an oracle-dependent hack into a sound method and unblocks bases nobody has published.

### Why this is tractable now: there is a per-value oracle

Guo-Yang's **CM-value tables** give the correct `s` at every listed discriminant — a per-value
oracle, far sharper than a per-model one. At `26_3` we know exactly which values are wrong and
what they should be:

    disc    ours     correct   relation
    -267    8/25     17/25     ours = 1 - correct   (s and s~ exchanged)
    -708    11/49    38/49     ours = 1 - correct
    other 12 discs: correct

### What is already ruled out — do not re-derive these

* NOT the sign tie-break. `find_signs_hauptmodul` never reports ambiguity at those discs (checked
  with an unconditional marker, so the silence is real and not a run that never reached it). The
  absolute values arrive already exchanged from `AbsoluteValuesAtCMPoints`.
* NOT seed or quaternion-order representative dependence. Stable across 3 seeds — so the `15_2`
  root cause ([[embedding-selection-root-cause]]) does NOT transfer.
* NOT a function of `d`'s standard invariants. `-267` matches the GOOD `-123`, and `-708` matches
  the GOOD `-132`, on class number, fundamentality, every valuation and every Kronecker symbol.
  ⇒ so a `Keep`-list keyed on `d` cannot work.
* NOT duplicate rows: each disc appears exactly once in our table.
* Coprimality alone does not predict it: 11 of the 14 discs are non-coprime to `N=3` and only 2
  are wrong.

### The live hypothesis, and the experiment that tests it

Both bad discs have `h(d) > 1`, so SEVERAL CM points share the discriminant. "We hold both values"
does not rule this out — evaluating `s` and `s~` at DIFFERENT points of the same discriminant
produces exactly the observed exchange. The good discs also have `h > 1`, so the hypothesis needs
a second ingredient: perhaps at the good ones the several points give the SAME value (conjugates
with a rational common value) while at the bad ones they differ.
⇒ EXPERIMENT RUN 2026-09-05. **The hypothesis is REFUTED and the problem is split.** Every
discriminant contributes **exactly ONE** rational CM point — `-267` (h=2) and `-123` (h=2) alike,
`-708` (h=4) and `-132` (h=4) alike. `h(d) > 1` does NOT put several points in the pool, so this
is **not point selection**.

### Where the defect is NOT — ruled out 2026-09-05, do not re-run these

    one CM point per disc          both bad and good; not selection
    ScaleForSchofer                -1/4 UNIFORMLY for -123/-267/-132/-708
    the optimal embedding lambda   all have the correct norm -2d, all with v3 = 1;
                                   no structural difference between bad and good
    find_signs_hauptmodul          never reports ambiguity at the bad discs
    seed / order representative    stable across 3 seeds

⚠ Incidental but useful: `ElementOfNorm` returns a DIFFERENT `lambda` for the same `d` on
different runs (`-123` gave `[-93,-10,-122]` then `[-67,-10,-90]`) while the computed values stay
fixed — so PR #18's isometry-invariance fix is holding, and lambda-dependence is NOT a live
suspect.

⇒ **THE DEFECT IS IN THE KAPPA / LOCAL-FACTOR SUM.** It is the only stage left between inputs that
are identical across the bad/good pairs and outputs that differ.

⚠⚠ **AND THE MISSING THEORY IS NOT MISSING — IT IS IN OUR OWN §10.** *(assaferan, 2026-09-05;
this supersedes the "the fix was never found" framing carried from
[[embedding-selection-root-cause]].)* Schofer Thm 4.1 assumes the local lattice is unimodular at
unramified primes, which FAILS at `p | N` for Eichler orders — and `sec:determined` (§10) supplies
the modification. From the proof of `thm:noresplit` (paper line ~1140): at each `p | DN` the local
discriminant group is `(Z/p)^2` and the two relevant forms are

    Eichler level p (UNRAMIFIED, non-unimodular):  hyperbolic plane  Qbar(b,c) = -bc/p
                                                   from Q = -x^2 - pyz
                                                   Gauss sum  p,  signature 0
    ramified (p | D):                              anisotropic plane (uc^2-b^2)/p
                                                   Gauss sum -p,  signature 4

    c_p^Eich = 1 if p|c, else +1/p          c_p^ram = 1 if p|c, else -1/p

so `c_p^ram = eps_p c_p^Eich`, and `rho_0 = prod_{p|D} c_p^ram * prod_{p|N} c_p^Eich`
(`thm:Ngeneral`).

**This is exactly the term the `15_2` investigation concluded was missing** — it ended at "the
+4log2 must come from a NON-kappa level-`N` bad-place factor that Schofer Thm 4.1 LACKS", and
noted that the DEAD `kappaminuszero` (`SchoferFormula.m:569`, `log_coeffs[p] := 1` for
`p | N/gcd(d,N)`) is "the vestige of this missing `p|N` factor". §10 is what that vestige was
waiting for.

⚠ **NOT a drop-in.** §10 gives these as MULTIPLICATIVE factors on the rho / constant-term side,
while the Schofer side is ADDITIVE in logs. So it supplies the missing local STRUCTURE and its
values; translating that into the `kappa` log contribution is the work. But this is now an
implementation-and-translation task against a written formula, NOT open theory.

### The next probe

Instrument `kappaminus` for `d = -267` (BAD) against `d = -123` (good) — an invariant-identical
pair — and compare **which primes vanish, and whether they vanish ALONE**. That is the quantity
that decides whether a prime contributes a log term, and at `15_2` it was exactly where the p=2
analogue went wrong (prime 2 co-vanishing with an odd prime -> double zero -> no contribution).
⚠ This needs the Borcherds forms for `26_3` (to supply the `m` values), so it is not a two-minute
probe; budget a run.

### Ordering, and why

1. **The experiment above** — cheap, and it splits the problem in two.
2. **Fix whichever half it indicts.**
3. Only then, sweep `CMNONCOPRIME` across the remaining blocked bases.
⚠ An earlier draft put "quantify how discriminating `ModelChecks` is" first. It does NOT accelerate
this: `ModelChecks` adjudicates whole models, while this defect needs per-value truth, and the
published CM tables already supply that. Keep it as the tool for validating models on bases with
no published equation — a separate job.

---

## MAIN LINE — the `A_m` theorem

One object blocks disproportionately much. Everything below this section is secondary to it.

- [ ] **State it precisely before attempting anything.** The vector-valued weight-3/2 Eisenstein
      coefficient `b^{eta*}_0(r)` at a *nonzero isotropic coset*, general `m`, supported on
      `N | r` — with `A_r = -b(r)/4`. A different object from the scalar, not a special case.
      ✅ **CHECKED 2026-09-04 against `sec:determined`, and the gap is CONFIRMED REAL.** The
      scorecard's phrase "Section `sec:determined` determines the functional outright" invites the
      hope that `A_m` is already in the paper and merely needs extracting. It is not. Two distinct
      things live there, and neither is `A_m`:
      * `thm:Ngeneral` determines `rho_mu` at isotropic cosets — but that is the **CONSTANT TERM**
        (`ct_L = glob * prod_p loc_p`), i.e. `m = 0`. Its own proof says the coset "contributes an
        **indicator, not a phase**": `rho_mu = kappa * rho_0` with `kappa = prod_{p in S}[p nmid c]`,
        because `Q(y) = 0` for isotropic `mu`. There is no `m`-dependence anywhere in it.
      * the "canonical representative — the Fourier coefficients of an explicit weight-3/2 form"
        is `-a_E`, the **SCALAR** of `prop:closedcoef` — which reproduces 39 forms' multipliers but
        only **1 of 13** measured `A_m`.
      ⇒ `m = 0` at those cosets: solved. All `m` for the scalar: solved. **All `m` AT a nonzero
      isotropic coset: genuinely absent.** `A_m` needs new mathematics, not extraction. This
      is now a checked claim rather than an assumption — do not re-open it hoping for a shortcut.
- [ ] **Load the regression set first.** 13 exact values already measured at `15_2`, `6_5`,
      `10_3`, `21_2`. Any candidate formula gets checked against all 13 before it is believed.
- [x] **Confirm the scalar route is closed.** *(done 2026-09-02)* Re-run with the repo's own
      `Hurwitz`: still 1 of 13. And structurally dead — at `15_2`, `m=10` and `m=30` give the
      *same* `a_E` but different `A_m`, so no sign convention or rescaling can rescue it.
- [x] **Three more routes closed, 2026-09-03 — all now dead ends, with reasons on record.**
    - *Insert KY Prop 5.4/5.5 into a level-`N` analogue of Theorem 8.1.* REFUTED, and now a
      permanent CI check (`tests/KudlaYangLocal.m`, `test_prop55_nonzero_isotropic_coset`, 1440
      checks): the level-prime local Whittaker *value* at a nonzero isotropic, `N`-only-supported
      coset is the constant polynomial `1` — no pole, no `m`-dependence — for every `m` and every
      base tested, including `N = 2`. Whatever produces `A_m`, it is not a local-factor product at
      `N`.
    - *Reduce `E_{A,eta*}` to zero-coset data via Schwagenscheidt's oldform relation
      (`eq:oldform`, `sec:ident`).* Already tried **in the paper itself** — re-read `sec:ident`
      before re-attempting this. It fails at weight 3/2 (needs analytic continuation, hits a
      non-holomorphic Zagier-type term; witnessed numerically as `-5` where the constant-term
      identity forces `0`).
    - *Push the `s`-law / genus-theta closed form (`sec:slaw`, `Theorem thm:rhoentry` +
      Siegel–Weil) to a closed form for `A_m`.* This is the derivation **behind** `prop:closedcoef`
      (`-a_E(m)`) — already refuted above under a different name. Don't re-derive it a second time.
- [ ] **The gauge ambiguity (`rem:gauge`) is real but resolvable — resolving it needs one more
      piece of data.** `-a_E(m)` and "Table A" (the `N|r`-restricted oracle fit) disagree at 6
      indices on `15_2` (`m = 1,2,3,10,15,30`) yet both fit the 9-form panel exactly, because that
      panel's principal-part matrix has rank 4 there (2-dim kernel; `gauge152.py`). **New
      2026-09-03: the 158-monomial data already computed in `cusp7_15_2.out` (campaign branch)
      gives that same 6-index matrix full rank 6** — the ambiguity is *not* one of the "50
      unremovable" directions `sec:exact` describes; it is a real gap the 9-form panel just didn't
      probe.
      ⚠ Using raw monomials as test inputs does NOT work directly: `c_{eta*}(0)` sums over *every*
      cusp class (`sum_w (rho(w^-1)e_0)_{eta*} a_0(f|w)`), and unlike a genuine Borcherds
      principal part (protected by [GY, Lemma 24] / `prop:nohalf`), an individual eta-monomial can
      have a nonzero constant term at *intermediate* cusps that a plain `(A_m, B_j)`-only model
      never sees. Restricting to monomial combinations with zero constant term at the intermediate
      classes actually on record in `cusp7_15_2.out` (`g = 3,4,5,12,15,20`, via its `PP` dump)
      still gives rank 6 at the gauge indices (152-dim reachable subspace) — but the numeric solve
      against real `c_{eta*}(0)` values is *still* inconsistent, because that dump is missing
      constant-term data at four more intermediate classes: `g = 2, 6, 10, 30`.
      ⚠⚠ **THAT "CONCRETE NEXT STEP" WAS FALSIFIED, 2026-09-04. Do not spend a session on it.**
      The step used to read: "re-run a `cusp7.m`-style dump so every intermediate cusp class gets
      its constant term recorded — not just the first one encountered — a scoped implementation
      task, not another open-ended search." Both halves are wrong:
      * The first-encountered-per-class dedup in `cusp7.m` IS a real latent bug and is now fixed
        (campaign `3d0720f`, keyed on `<class, monomial>` instead of class alone) — **but fixing
        it changed nothing.** The re-run produces a BYTE-IDENTICAL `PP` section, still missing
        `g = 2,6,10,30`.
      * The actual reason those classes are absent: **every coset word in them has `lead > 0` for
        every one of the 158 monomials** — measured, `minlead` 12–36, `#le0 = 0` at every such
        word. The monomials are simply holomorphic there. All 12 divisor-classes of `M = 60` DO
        occur among the 144 words, so the classes exist; the POOL has no support at them.
      ⚠ **The "provably resolvable" premise is NOT established** *(checked 2026-09-04)*. The
      rank-6 figure is real but is computed over MONOMIALS, which are not valid probes of
      `c_{eta*}(0)` — exactly the caveat the paragraph above states. Whether the six-index gap is
      resolvable over the FORM space is uncomputed. Treat this section's premise as open, not
      verified, and see SHIP.
      ⇒ It is a **search problem, not a dump fix**: it needs eta-monomials that actually carry a
      pole or constant term at cusp classes 2, 6, 10 and 30. `extraseqs`' 69 hand-picked vectors
      came from an earlier session with no surviving generator, so that pool has to be
      reconstructed or extended first. Write-up:
      `vvdata/weyl-campaign/note-cusp7-missing-classes.md` on the campaign branch.

**Why this one.** It *is* the `Kappa0` log-`N` defect at `SchoferFormula.m:589` — the same defect
`coprime_to_level` (`ShimuraQuotients.m:1420`) papers over with what the code itself calls "a blunt
instrument". Closing it unblocks the 49 obstructed bases in a single move.

**Honest risk.** No product of local densities reproduces `b` under any sign convention, with KY
Props 5.3/5.4 validated at 470+150 checks (now including Prop 5.4/5.5 at nonzero isotropic cosets,
2026-09-03). Three separate derivation routes are now closed (above). This may need new
mathematics — or, more likely given the gauge finding, one more piece of cusp-class data.

---

## SHIP — the paper

Runs in parallel. The fast arc is currently hostage to the slow one — that is the thing to fix.

- [ ] **Decide explicitly: does the paper need model results?** The blocking decision. If it does,
      scope exactly which claim — what exists is four models from one base, a partial set.
      ⚠ **This decision is time-sensitive as of 2026-09-04**: `34_11` is being run through the
      pipeline right now and may produce a second base. Wait for it, or decide conditionally —
      "four models from one base" reads very differently from "two bases".
- [x] **Re-read `conj:sq` under the gauge finding — DONE 2026-09-04, and the paper is already
      consistent.** It attaches the caveat in BOTH places it needs to: after the scorecard
      (~line 673) and in `rem:gauge` (~line 1573), each saying `w_square` "is a coordinate of a
      gauge-fixed representative rather than an invariant" and that `conj:sq` "is to be read as a
      statement about that representative". So the answer to "is it still well-posed?" is: it is
      well-posed *as a statement about the gauge-fixed representative*, the paper says so, and its
      tested predictions stand unfalsified (`14_5` and `34_5` both hit).
- [x] **The `rem:gauge` / MAIN LINE tension: CLOSED — the paper is right. My 2026-09-04 "the paper is wrong"
      claim was RETRACTED the same day — do not act on it.** What was measured is solid:
      `vvdata/weyl-campaign/rankcheck_gauge.py` (campaign) shows the `oo`-side matrix over the
      **158 MONOMIALS** at `m = 1,2,3,10,15,30` has **rank 6**, against **rank 4** for the
      nine-form panel, with two hard-assert validations (form `-1`'s principal part; the panel's
      rank 4).
      ⚠ **But it does not bear on the paper's claim — and the REASON is not the one MAIN LINE
      gives.** *(reasoning corrected 2026-09-04, second pass)* MAIN LINE objects that monomials
      can carry constant terms at intermediate cusps. `sec:exact` says the opposite and is
      explicit about it: each monomial "is itself a weakly holomorphic weight-`1/2` form of the
      same multiplier character, **holomorphic at every intermediate cusp**", and the joint
      `(oo,0)` factorisation is measured to hold to `10^{-11}` on that space. So intermediate
      cusps are not the obstruction.
      **The actual obstruction is that I computed in the wrong MODEL:**
      * for genuine **FORMS** the `oo`-ONLY model is valid — `(1/2) sum_m c(-m) W(m)` reproduces
        the measured multiplier on all **39 forms** across nine bases;
      * for **MONOMIALS** the `oo`-only model provably FAILS — `sec:exact` measures residual
        **2.08**; the `0`-side genuinely enters.
      `rem:gauge`'s rank-4 claim lives in the `oo`-only model **on forms**, where it is
      legitimate. My rank-6 lives in the `oo`-only model **on monomials**, where the paper has
      already measured that model as inadequate. Ranking an inapplicable model says nothing.
      ⇒ **`rem:gauge` is NOT known to be wrong, and the paper should NOT be edited on this basis.**
      The question it answers — can more FORMS separate `-a_E` from Table A? — is one nobody has
      computed.
      **RESOLVED 2026-09-04 (third pass), IN THE PAPER'S FAVOUR — `rem:gauge` looks CORRECT.**
      The tempting next step was "redo the rank over forms, e.g. a weakly-holomorphic basis".
      That would have been the SAME MISTAKE a third time. `sec:exact` states that the 158
      monomials ARE "weakly holomorphic weight-`1/2` form[s] of the same multiplier character" —
      and the `oo`-only model demonstrably FAILS on them (residual `2.08`). So being such a form
      is **not sufficient** for the `oo`-only functional to be meaningful; a weakly-holomorphic
      basis is not the right domain either.
      The `oo`-only formula `(1/2) sum_m c(-m) W(m)` is validated only on **genuine Borcherds
      forms from the divisor construction** (39 forms, nine bases). At `15_2` those are the NINE
      panel forms, whose `oo`-principal parts have **rank 4** at the six indices. Within the
      domain where the model means anything the 2-dimensional ambiguity is therefore REAL, and
      enlarging the probe set does not help because every enlargement leaves that domain —
      which is exactly what the paper's hedge "further panels **of the same kind**" already says.
      ⇒ **Do not edit `rem:gauge`, and do not re-run this as a rank computation over any larger
      space.** The one honest caveat left is narrow: this assumes the nine panel forms are ALL the
      Borcherds forms at `15_2`; if that base admits more (more covers/keys) their principal parts
      would add rows and the rank could rise. That is the only version still worth asking.
      **Method note, three passes in:** monomials, then weakly-holomorphic forms, then the right
      answer. Every time the arithmetic was fine and the DOMAIN was wrong — and the validations
      that made the numbers trustworthy (they caught Magma line-wrapping and dropped rational
      coefficients) said nothing about which object was being measured. On this problem spend the
      effort on "which objects is this functional even defined for", not on the linear algebra.
- [ ] **Rebuild the PDF, confirm refs resolve, submit.** *(partly done 2026-09-04)* The PDF is
      **current** — `.tex` and `.pdf` were last committed in the same commit and the `.tex` has not
      moved since — and **refs already resolve clean**: 0 undefined references, 0 undefined
      citations, 0 LaTeX warnings in the committed build log, and 0 occurrences of `??` in the
      PDF text. So no rebuild is needed unless the `.tex` changes; what remains is the decision to
      submit.

---

## REPAIR — the `k0` coset bug

New as of 2026-09-02, sharply localized, and it may free the last two runnable bases.

- [x] **Kill the `10_61` run.** *(done 2026-09-02)* Stopped at 17 h 26 m, RSS 5.6 GB and still
      climbing. Nothing flushed on exit, so the 924 recorded failures are all we have — expected,
      per the Magma-buffering trap.
- [x] **Find what distinguishes `wi = 2` and `wi = 1221`.** *(largely answered — see
      "What this turned out to be", below)* Both sit at the **extremes** of the `Im(z0)`
      distribution: `wi = 1221` is rank 1 of 2232 (smallest, 8.81e-7), `wi = 2` is rank 2231
      (0.7229, nearly the largest).
- [x] **Why does `wi = 2` fail? CLOSED, quantitatively.** *(done 2026-09-03)* `eta` is evaluated
      at `d*z0` for `d` in `ds = Divisors(M)`, so `d` runs to 1220. A LARGE `Im(z0)` is therefore
      just as fatal as a small one: `Im(d*z0)` reaches `1220 * 0.7229 = 882`, and since
      `log10|eta(z)| ~ -0.11374*Im(z)`, that predicts `-100.31`. **Measured `-100.272`** — a
      four-figure match. The `eta` values span `[10^-100.27, 10^-0.08]`, i.e. a **100.19-digit
      dynamic range at `Prec := 80`**, so the ratio is destroyed in floating point even though it
      is mathematically fine.
- [~] **Why does `wi = 1221` fail? — LARGELY MOOT as of 2026-09-04.** The per-coset `tau`
      fallback now *routes around* both extremes rather than needing them explained: `wi = 1221`
      sits at `Im(z) ~ 1e-6`, below the `1e-5` floor, so it gets a corrected evaluation point
      automatically. The question below is still of scientific interest but is no longer on the
      critical path, and the measurement it asks for is no longer a prerequisite for anything.
- [ ] *(original item, kept for the caveat it carries)* **why does `wi = 1221` fail?** The roles have swapped. Its measured
      spread is **1.39 digits — identical to typical cosets** (`wi` 300, 611, 900 all give
      1.39436), so the dynamic-range mechanism does *not* explain it.
      ⚠ **But that number is not trustworthy**: at `Im(z0) = 8.8e-7`, *both* `z` and `-1/z` have
      small imaginary part, so the one-`S`-step approximation used in `tauwindow.m`'s `logeta` is
      invalid there. **Re-measure with a real `DedekindEta` (full `SL2` reduction) before
      concluding anything about this coset.**
- [x] ~~Cross-check against the `39_2` malformed-form pathology.~~ **Superseded.** This is not a
      malformed form. It is an evaluation-point failure of the harness — see below.
- [ ] **Re-run `14_43` AND `10_61` — but now with the fix, not with GATE3B.** *(reframed
      2026-09-04)* The cause IS understood: extreme `Im(z)` at a handful of cosets, and the
      per-coset `tau` fallback addresses exactly it. So the useful run is no longer a GATE3B
      report-and-continue measurement — it is the ordinary pipeline, to see whether they now pass.
      `10_61`'s two poisoned cosets are at the two extremes the fallback catches. ⚠ `10_61` is
      `M = 1220` with 2232 cosets and previously ran 17 h before being killed; budget accordingly,
      and use `M0PROGRESS=1` so a stall is diagnosable instead of a silent buffered log.
- [ ] **Reclassify both bases in the sweep record — but wait for the re-run above.** `10_61` was
      called "not runnable, a real upstream defect". With the `tau` fix that verdict may simply be
      wrong; do not rewrite the sweep record until the re-run says which.

Raw evidence: `~/shimura/models/10_61.gate3b.log` on lovelace (353 KB, 926 lines).

### What this turned out to be: precision is `M^2`

`M0MultiplierExact` evaluates the slash constant at two HARDCODED points
(`VectorValuedForm.m:419-420`) and pushes every coset word through them. Measured 2026-09-02:

    base     M      #words   min Im(z0)   Im(tau0)/M^2   ratio    known precision
    15_2      60      144     3.7222e-4     3.6389e-4    1.0229   exact
    58_5     580     1080     3.9034e-6     3.8942e-6    1.0024   33 digits
    34_11    748     1296     2.3457e-6     2.3414e-6    1.0018   18 digits
    10_61   1220     2232     8.8114e-7     8.8014e-7    1.0011   CATASTROPHIC

**`min Im(z0) = Im(tau0)/M^2 * (1 + O(1/M))`.** Elementary cause:
`Im((a*tau+b)/(c*tau+d)) = Im(tau)/|c*tau+d|^2`, and the coset reps carry `c, d` up to `O(M)`.
Not word length — `maxlen` saturates at 13 across all three large bases while `Im` keeps falling.

The precision column is this repo's own recorded numbers, measured independently and earlier. They
track the `Im` collapse in lockstep, so **"achievable precision is base-dependent" now has a cause:
`M = 2DN` (`4DN` odd), through the coset denominators.** Three consequences:

* **Predictive with no pipeline run.** `VVCosetReps`/`VVSTWord`/`slashdata` depend only on `M`, not
  on the Borcherds forms — 7 min at `M = 1220`. A new cheap triage axis.
* **The fix is not a tolerance.** The slash constant is *independent of the evaluation point* —
  that independence is exactly what the two-point check tests — so `tau` is a FREE CHOICE. Pick it
  per coset so `Im(w*tau)` stays in a safe band, instead of forcing two global constants through
  2232 words.
* **`10_61` is not defective arithmetic**, so the "two runnable candidates of 122" count is wrong
  on other grounds too.

### The safe window, measured

Worst `eta` dynamic range per `Im(z0)` band, over all 2232 cosets of `M = 1220`:

    Im(z0) band          #cosets   worst spread (digits)
    [0,      1e-5)         1137          92.17
    [1e-5,   1e-4)          669          37.45
    [1e-4,   1e-3)          307           9.90
    [1e-3,   1e-2)           91           3.55   <- SAFE
    [1e-2,   1e-1)           22           9.25
    [1e-1,   1e0 )            5         100.26
    [1e0,    1e2 )            1         181.56

**Both ends are fatal and the middle is safe** — the usable band is roughly
`Im(z0)` in `[1e-3, 1e-2)`, where the worst case over 91 cosets is 3.55 digits against `Prec := 80`.
Only **six** cosets have `Im(z0) > 0.1`, and they carry the worst spreads in the whole set.

⇒ **The fix, concretely.** `tau` is a free choice (the slash constant is `tau`-independent — that
independence is exactly what the two-point check tests). So instead of forcing two global constants
through 2232 words, choose `tau` per coset so that `Im(w*tau)` lands in the safe band. This is a
local change to `M0MultiplierExact` and needs no theory.

- [x] **Implement per-coset `tau` selection — DONE 2026-09-04 (`475e72b`).** Validated three
      ways: `15_2` stays exact (9/9), `58_5` keeps its models (full suite clean, `ModelChecks`
      8241 checks 0 failures), and BOTH `34_11` and `58_5` match the Prop 9.15 closed form 9/9 —
      an independent route with no complex arithmetic, no evaluation points and no `tau` at all
      (`58_5`'s was a PRE-REGISTERED prediction). `34_11` went from failing to passing. Details
      and the two hypotheses it refuted along the way are below.
      **2026-09-03: first attempt tried and reverted (not committed).** Picked two fixed target
      points `z0targ, z1targ` in the safe band and set each word's `tau0, tau1` to be the
      PREIMAGE of those targets under the word's own SL2(Z) matrix (`tau := g^-1.ztarg`), reusing
      `slashdata`/`sfun` unchanged. Passed `15_2` 9/9 exactly. But on `58_5` the run (just
      `M0MultiplierExact`, not the full pipeline) did not finish in 35+ minutes where the OLD
      fixed-tau version completed as part of an 18-minute full pipeline run — a real slowdown, not
      just noise. Suspected cause, UNVERIFIED: `sfun`'s own argument is `(tri[i][1]*tau+tri[i][2])/
      tri[i][3]`, evaluated at the SAME `tau` (not at `z`) — fixing `z` in the safe band by taking a
      preimage under the word's full matrix can push `Im(tau)` itself to extreme values (small or
      huge) depending on that word's own `(c,d)`, and `DedekindEta` at extreme `Im` may need much
      more internal work to hit `Prec := 80`. **Next attempt should control `Im(tau)` directly**
      instead of deriving it as a free variable: fix `Re(tau) = x0` and solve the exact quadratic
      `t = y/((c x0+d)^2 + c^2 y^2)` for `y = Im(tau)` on the LARGE-`y` branch (`y = [1 +
      sqrt(1-4c^2t^2(cx0+d)^2)] / (2c^2t)`, `t` the target `Im(z)`), so both `Im(z)` (via `t`) and
      `Im(tau)` stay in a controlled, moderate range simultaneously — never solved this way and
      never measured. Reverted file kept at hand in case useful:
      `/private/tmp/claude-501/.../scratchpad/VectorValuedForm.m.new` (session-local, may not
      survive) — re-derive from this note rather than relying on that path.

      **2026-09-03, SECOND attempt — WORKS, and the diagnosis above was WRONG.** Two corrections
      to the paragraph above, both measured:
      1. **`DedekindEta` at small `Im` is NOT slow.** Timed directly at `34_11`'s own points:
         `DedekindEta` at the old `tau` (Im 1.31) and at the new one (Im 0.0057) both take
         0.000s, individually and across all 12 divisors. The whole "extreme `Im(tau)` makes
         `DedekindEta` expensive" theory is refuted; do not re-derive it.
      2. **The real cost is ONE word, and it is `tau`-INDEPENDENT.** With `M0PROGRESS=1` (new,
         see below): at `34_11`, word `wi=1` — the IDENTITY coset — has `W=1, depth=13993`,
         while *every other selected word* has `depth=1`. With `W=1` the triangularisation does
         no reduction, so `leads = sum r_i d_i` runs to -13992 (divisors to 748, exponents to
         29) and the code raises ~14000-term power series to powers up to 29, once per monomial,
         over 799 monomials. `W`, `leads`, `depth`, `units` are all computed from `tri` BEFORE
         any `tau` logic runs, so this cost is identical in the unmodified code. `58_5` is the
         same shape (`wi=1`, `W=1`, `depth=10585`, 603 monomials).
         ⇒ The first attempt was never "slow because of the preimage"; it was slow because
         **fixing the evaluation points lets the computation actually RUN TO COMPLETION**, where
         the old code aborted early on the two-point check. There is an obvious optimisation
         here (`units[i]^(r[i])` is recomputed per monomial; cache by `<i, exponent>`), NOT done.

      **What the working fix is.** Keep the global `tau0/tau1` for most words; fall back to the
      per-word preimage ONLY where the default lands `z = w.tau` outside the safe range. The
      threshold must SCALE WITH `M` on the upper side: `Abs(Im(z)) ge 1e-5 and M*Abs(Im(z)) < 100`.
      Using the raw `M = 1220`-calibrated absolute band `[1e-5, 1e-1)` is a real bug — at `15_2`
      (`M = 60`) it flags the DEFAULT `tau0` itself (`Im = 1.31`) as unsafe, on a base that has
      always been exact with that `tau0`. Mechanism: `eta` is evaluated at `d*z` for `d` up to
      `M`, so large-side risk is about `M*Im(z)`, and the `M = 1220` measurement's boundary
      `Im(z) ~ 0.1` is `M*Im(z) ~ 122`. The small side is about the word's own S-structure, not
      `d`-scaling, so it stays an absolute floor. Fallback rates measured: `15_2` 1/64,
      `58_5` 8/64, `34_11` 13/64 (triggered on BOTH mechanisms — `wi=1` at Im 1.31/1.73 on the
      high side, `wi=376,749,936,1123,1296` at Im ~ 2-9e-6 on the low side).

- [x] **`34_11` PASSES, and its multipliers are confirmed by an INDEPENDENT route.**
      *(2026-09-03)* With the fix, `34_11` completes `M0MultiplierExact` (8604 s) and returns
      `-1/2, 1/2, -5, -2, -3, -8, 0, -2, -6` for keys `-2,-1,8..14`, clearing all five internal
      guards. Cross-checked against the Prop 9.15 closed form (`closedcoef.py`'s `W(D,N,m)`,
      `mult = (1/2) sum_m c(-m) W(m)` on the oo-side principal parts): **9/9 exact, fractions
      included.** The closed form is exact rational arithmetic with no complex numbers, no
      evaluation points and no `tau` at all, so this is a genuinely independent confirmation that
      the fix produces CORRECT values, not merely self-consistent ones. Method validated first on
      `15_2`, where the same comparison also gives 9/9 against the known panel.
      ⚠ Scope: the oo-only closed form is itself "verified on nine bases / 39 forms but NOT known
      to be valid in general" (`note-39_2.md`), so this corroborates `34_11`, it does not prove
      the fix universally correct.

- [ ] **⚠ REVISIT: "gate 4 is a GENUINE violation" is probably WRONG.** This file and `HANDOFF.md`
      both record `34_11` as clearing gates 1-3 and then failing class-constancy at
      `dev 0.0186, scale 0.0430` — "a 43% deviation... NOT another tolerance... Loosening it would
      manufacture a multiplier wrong by 43%." But with per-coset `tau` and the SAME `1e-15`
      tolerance (no guard was touched), `34_11` passes class-constancy and lands on the closed
      form 9/9. Class-constancy compares up to 3 sampled cosets per class, and at `34_11` several
      sampled cosets sit at `Im(z) ~ 1e-6`, where the measured loss is ~92 digits — i.e. the old
      fixed `tau` fed the check numerically meaningless values at exactly those representatives,
      which is precisely how a spurious O(scale) spread appears. Working hypothesis: gate 4 was
      never real; it was the gate-3 evaluation-point defect resurfacing one gate later. Confirm
      (or refute) before relying on the "gate 4 is real" framing anywhere else.

- [ ] **New diagnostic: `M0PROGRESS=1`.** Env-gated, silent by default, `WriteStderr` (NOT
      `printf` — buffered output is lost when a run is killed, which cost this session a 78-byte
      log after 90 minutes). Prints words selected, per-word `W`/`depth`, every fallback with its
      `Im(z)`, per-monomial progress, and a final fallback tally. Use it before assuming a long
      `M0MultiplierExact` run is stuck.

Probes: `vvdata/weyl-campaign/tau-precision/` on the campaign branch — `tauprobe.m` (per-coset
`Im`), `tauscale.m` (the `M^2` law), `tauwindow.m` (the dynamic range and the band table), plus
`amtest.m` (the scalar-`a_E` comparison). ⚠ `tauwindow.m`'s `logeta` is a **one-`S`-step
approximation**: valid when `Im(z)` is large or one `S` makes it large, INVALID when both `z` and
`-1/z` sit near the real axis. That is exactly the `wi = 1221` regime, so its small-`Im` numbers
are provisional.

---

## MEASURE — route B's `k = 3/2` phase

The only cheap-predictor path left, and it converts 81 unknowns into data.

- [ ] **Fix the half-integral convention *from theory*, not by fitting.** Do not tune the constant
      — a fitted phase would silently corrupt all 81 classifications.
      **2026-09-04, MEASURED — the target is much sharper than "the error is nearly `-d/6`".**
      Ran `weildim2.m` (cross-check at `6_1` still PASSES, all 10 traces). Decomposing
      `dimM = d + dk/12 - a1 - a2 - aT` by its leading rationals, at `k = 3/2`:

          base      d       a1/d      a2/d      aT/d     dimM   truth
          34_3    20808   0.37495   0.41671   0.49986   -3465     0
          38_5    72200   0.37499   0.41668   0.50237  -12204     1
          38_7   141512   0.37499   0.41667   0.50009  -23598     0

      `a1/d -> 3/8` and `aT/d -> 1/2` on the nose; `a2/d -> 5/12`. And
      `1 + k/12 - 3/8 - 5/12 - 1/2 = -1/6` **exactly** — so the `-d/6` is not a vague "phase is
      off", it is **entirely attributable to `a2`, the ST / order-12 term**. `a1` (the order-8
      S-term) and `aT` are already right.
      Sharper still: had `a2/d` been `3/12` instead of `5/12`, the leading term would be
      `1 + k/12 - 3/8 - 1/4 - 1/2 = 0` **exactly** — i.e. `a2`'s index looks offset by **2 units of
      the 12-cycle**. Since `a2 = sum_j m_j (j/12)` with `sum_j m_j = d`, relabelling `j -> j+2`
      moves it by `2d/12` MINUS the wrap at `j = 10, 11`, which is base-dependent — and that is
      exactly why the residual is "nearly" `-d/6` and differs per base (`-0.16652`, `-0.16904`,
      `-0.16676`) rather than being constant.
      ⚠ **What this does NOT establish.** It shows the LEADING term vanishes under that shift. The
      actual answer (`0, 1, 0`) lives in the `O(1)` deviations from those rationals — `a1/d`
      differs from `3/8` only in the 5th decimal — so the systematic `O(d)` offset currently
      swamps the entire signal. Getting the leading term to cancel is necessary, not sufficient,
      and **applying a `+2/12` shift by hand would be exactly the fitting this item forbids**: it
      would make the big number vanish while leaving the `O(1)` content unvalidated. Use this to
      aim the derivation — the metaplectic base point for the order-12 `ST` lift — and then check
      the `O(1)` part against `38_5 -> 1`, `38_7 -> 0`, `34_3 -> 0`.
- [ ] **Validate against the known deficits:** `38_5 -> 1`, `38_7 -> 0`, `34_3 -> 0`.
- [ ] **Classify the 81**, including whether the obstructed class is larger than 49.

The machinery is already sound: the trace formulas reproduce 11 classical `dim M_k(SL2(Z))` values
exactly. Only the half-integral phase is wrong.

---

## COVERAGE — reproducing Guo-Yang's published equations

*(stock-take RE-MEASURED 2026-09-05, superseding the 2026-09-04 counts. Distinct from the
CM-value tables in `tests/_offline/` — this is about the paper's headline output, the EQUATIONS.)*

    43   (D,N) bases with published equations in Guo-Yang
    36   we have a model for            <- was 35; 93_1 RECOVERED 2026-09-05 (the vx fix)
    35   ...and a test comparing it to Guo-Yang   <- GuoYangEquations.m + GuoYangCurve_14_3
     1   model exists but NO GY comparison test   <- 22_5 only
     7   no model: the real blockers   <- 15_4 26_3 69_1 95_1 111_1 119_1 159_1
         ⚠ 15_4 is NOT one of them -- see below; it is outside the method by the authors' own
         statement, so the honest target is 42, not 43.
    47   models we have BEYOND Guo-Yang's list entirely (83 model files total)

⚠ **TWO DIFFERENT COUNTS, do not conflate them** (this bit on 2026-09-06 — `93_1` was first
reported as "34 -> 35", which was wrong). **34** is the number of GY **CM-VALUE tables** we have as
offline tests (`tests/_offline/GuoYang_*.m`); that set does not even contain `93_1`. **43** is the
number of bases with published **EQUATIONS**, which is what this section counts. The 43 are
reproducible from the source: extract every `$X^D_0(N)$` label in lines 3150-3800 of
`vvdata/weyl-campaign/guoyang/ShimuraCurves-arxiv.tex` (⚠ starting at 3300 silently drops five
bases — `35_1 38_1 39_1 51_1 55_1` — whose table row begins at 3264) and intersect with
`data/models/`.

**⇒ TIER 1' IS EXHAUSTED AS A TRANSCRIPTION TASK.** Eight of the original ten are done
(`14_5 15_2 21_2 22_3 51_1 55_1 57_1 87_1`, in `tests/GuoYangEquations.m`). The remaining **two are
not transcribable at all** — we do not have the object to compare:
* `14_3` — GY publishes a PAIR (`z^2=-9x^2-2`, `y^2=-7x^4+22x^2+1`). Our `models[[1]]` entry is
  literally `[* *]`, **empty**. A fresh plain `genmodels` run on 2026-09-05 reproduced the
  committed file byte-for-byte in 37 s, so the empty `{1}` is a REPRODUCIBLE outcome of
  `AllEquationsAboveCovers`, not a stale or truncated artifact.
* `22_5` — GY publishes `y^2 = -11x^12-80x^10-240x^8-362x^6-240x^4-80x^2-11` (degree 12, genus 5).
  Our file has **no `[1]` key at all**, only three quotients; its commit message (`08ce5fa`) says
  "3 covers", i.e. it was committed knowing the full curve was absent.
⇒ Both are **model-generation** items, not test-writing items. Do not re-file them under TIER 1'.

#### `15_4` IS OUT OF SCOPE -- BY GUO-YANG'S OWN STATEMENT (settled 2026-09-06)

**Do not attempt to make the pipeline produce `15_4`.** The published version of the paper
(Compositio Math. 153 (2017) 1-40) carries a remark that **arXiv v1 does not have**:

> *Remark 39.* Note that there is a curve, namely, `X = X_0^15(4)`, whose equation is not obtained
> using our method. This is because **the normalizer of the Eichler order in this case is larger
> than the Atkin-Lehner group**. For this special curve, we use the result of Tu [Tu14].

They then take a Hauptmodul `t_4` on `X/<w_3,w_5>` from [Tu14, Lemma 13] (values `+-1/sqrt(-3)`,
`+-sqrt(-15)/5`, `(+-1 +- sqrt(-15))/8` at discriminants `-12, -15, -60`), read off
`y^2 = a(4t_4^2-t_4+1)(4t_4^2+t_4+1)(5t_4^2+3)` and `z^2 = b(3t_4^2+1)` from the ramification, and
fix `a = b = -1` using `t_2 = (5t_4^2+2t_4+1)/(7t_4^2-2t_4+3)` on `X_0^15(2)/<w_3,w_5>` plus
Schofer. `X_0^15(4)` is also one of the hyperelliptic Shimura curves that are **not hyperelliptic
over R** [Ogg83].

⇒ Since `N^+_B(O)` strictly contains `W_{15,4}`, the star quotient our pipeline forms is the WRONG
OBJECT. This is the hypothesis flagged in memory as never checked (`O^+_{L,F} = O^+_L`), and here
it demonstrably fails.

**What we measured before finding the remark, kept because it is still true of the CODE:** two
independent layers encode squarefree-`N`. (1) `assert IsSquarefree(N)` in `get_D0_M_g`. (2) Behind
`NONSQFREE=1` the basis builds fine (4.4 s), and the run then stops at `ShimuraQuotients.m:1431`,
whose star check compares `W` against `Set(Divisors(D*N))` -- but the AL group is indexed by HALL
divisors, and `DN = 60` has 12 divisors to 8 Hall divisors, so it CANNOT pass (same at `10_9`,
`21_4`). Fixing that indexing would still not produce `15_4`, because of Remark 39.

#### CAN THE `15_4` ROUTE BE GENERALISED? — scoped 2026-09-06

**The target set is 9 bases, not hundreds.** Of 246 `D>1` non-squarefree-`N` entries in
`GetHyperellipticCandidates`, only **9 carry a `W={1}` full curve** and are therefore real model
targets:

    6_25 (g5)  6_49 (g9)  10_9 (g5)  14_9 (g7)  15_4 (g5, DONE)  15_8 (g9)  21_4 (g7)
    22_9 (g11) 33_4 (g11)

**What generalises: the VALIDATION half, and it is worth having.** `ComputePointsViaTrace` and the
Shimura genus formula run fine on non-squarefree bases — that is exactly what pinned `15_4`'s two
constants, and it is an oracle for a whole class the pipeline cannot otherwise touch. Any candidate
model for these 9 can be checked the same way.

**What does NOT generalise: the GENERATION half.** `15_4` worked because Tu (Pacific J. Math. 269
(2014), Lem. 13) published a Hauptmodul with explicit CM values, which Guo-Yang quote. Nothing
supplies that for the other 8. Two possible routes:
* **the Hall-divisor refactor** (index `W` by Hall divisors; lift `assert IsSquarefree(N)`), which
  only helps where Guo-Yang's obstruction is ABSENT;
* **external Hauptmoduls** — Tu's paper is about genus-zero Shimura curves generally, so it may
  cover more of these. We do not have it; it is not in the user's Dropbox.

**⚠ A PREDICTION, NOT A MEASUREMENT — which of the 8 are only CODE-blocked.** Guo-Yang's
obstruction is that `N^+_B(O)` strictly contains the Atkin-Lehner group. By analogy with
Atkin-Lehner-Newman for `Gamma_0(N)` — where the normalizer exceeds the AL group exactly when
`h > 1` for `h` the largest divisor of **24** with `h^2 | N` — one expects:

| base | `N` | `h` | expected |
|---|---|---|---|
| `6_25` | 25 | 1 (5 does not divide 24) | **no obstruction — possibly only code-blocked** |
| `6_49` | 49 | 1 (7 does not divide 24) | **no obstruction — possibly only code-blocked** |
| `15_4`, `21_4`, `33_4` | 4 | 2 | obstruction — matches Guo-Yang's Remark 39 ✓ |
| `15_8` | 8 | 2 | obstruction |
| `10_9`, `14_9`, `22_9` | 9 | 3 | obstruction |

The single confirming data point is `15_4` itself, where the criterion agrees with Guo-Yang.
⚠ **I could not COMPUTE this.** `TwoSidedIdealClassGroup(O)` returns 1 for all nine INCLUDING
`15_4` — it fails the positive control, so it is not measuring `N^+_B(O)/Q^*O^*` and its output must
not be quoted. Until someone computes the normalizer properly (or reads Michon/Ogg on normalizers
of Eichler orders), the table above is an **analogy with one confirmation**, not a result.

⇒ **Recommended if this is pursued:** test the prediction at `6_25` or `6_49` — they are the only
two where the payoff (a base reachable by fixing code) justifies the Hall-divisor refactor. Confirm
the normalizer claim FIRST; the refactor is wasted if the obstruction is present anyway.

#### ⚠ WE HAVE BEEN READING THE SUPERSEDED VERSION

arXiv has exactly **one** version (v1, 2015-10-21). The paper of record is **Compositio Math. 153
(2017) 1-40**, substantially revised and **not on arXiv**; the copy at
`vvdata/weyl-campaign/guoyang/ShimuraCurves-arxiv.tex` is v1. Differences found so far:

| | arXiv v1 | journal |
|---|---|---|
| `93_1` equation | `3s^3 - 7s^2 - 3t - 1` (typo) | `3s^3 - 7s^2 - 3s - 1` |
| `39_2` involutions | `w_4, w_3, w_5` (copy-paste from `15_4`; `4` does not divide `78`) | `w_2, w_3, w_39` |
| Remark 39 (`15_4` out of scope) | absent | present |
| Remark 38 (`10_19` NOT hyperelliptic over Q; `14_5` is) | absent | present |

The journal PDF is at `~/MIT Dropbox/Eran Assaf/Research/HyperellipticQuotients/Guo, Yang -
equations-of-hyperelliptic-shimura-curves.pdf`. **All ten equations in
`tests/GuoYangEquations.m` were re-verified against it on 2026-09-06 and agree.** Remark 38 also
bears on the recorded `10_19` anomaly (its v1 table contradicting the paper's own worked example)
-- re-check that against the journal before trusting either.

#### What blocks `14_3` and `22_5`, measured 2026-09-05

Cheap and now on record, so nobody re-runs it. Both bases fail the same way and it is **not** the
back-fill stage:

* **The covers are UNDER-DETERMINED, not short of CM points.** `EquationsCovers.m:157-158` accepts
  a cover only when `solve_quadratic_constraints` returns a *unique* solution; `#coeffs ne 1`
  raises, is caught at `:166`, and the cover is deferred. A deferred cover is then recoverable
  ONLY as a quotient of an already-determined cover above it (`:701`). At `14_3`, 6 of the 7
  index-2 covers defer, nothing above them is determined, and all 6 emit
  `could not recover deferred cover` — so `W={1}` has no substrate to be built on. `W={1}` is
  never itself mentioned in the log: not determined, not deferred, not recovered.
* ⚠ **`cmsupply` says `OK`, and that verdict DOES NOT APPLY to the full curve.** Measured:
  `BASE 14 3 demand 7 genera [0,0,0,0,1,1,1] CMVERD OK margin 0` and
  `BASE 22 5 demand 9 genera [0,0,0,1,1,1,2] CMVERD OK margin 0`. Those genus lists top out at
  **1 and 2**, but Guo-Yang's published `14_3` is genus **3** and `22_5` is genus **5** — the 7
  curves it measures are the *index-2* covers of the star curve (their genera match the seven
  `W={1,6,7,42}`-type rows in the genmodels log exactly, which is what pins the identification),
  and `demand 7 = 2*1+5` comes from a max genus of 1. This is structural, not a coincidence of one
  base: the supply check iterates `Xstar`CoveredBy` (`ShimuraQuotients.m:1526`), the *immediate*
  covers, exactly as `AllEquationsAboveCovers` does (`EquationsCovers.m:740`). So `OK` means "supply is adequate for the
  easy targets, with **zero** margin"; it is silent on `W={1}`. Do not quote `CMVERD OK` as
  evidence that these bases have enough CM points. (Same wrong-object family as the rank-over-
  monomials error of 2026-09-04.)
* **`INTSOL=1` is REFUTED as a lever at `14_3`** — plain and `INTSOL=1` runs produce **byte-identical**
  model files (37 s / 38 s), and both are byte-identical to the committed `data/models/models_14_3.m`.
  That byte-identity is also the reproduction check: the committed file is current, not stale.
* ⚠ **`22_5`: do not overwrite `data/models/models_22_5.m` — but NOT for the reason first recorded.**
  *(diagnosis corrected 2026-09-05 after measuring; the first version of this bullet was wrong and
  is kept only as the retraction.)* A fresh plain run (680 s) yields **8 cover-keys, 2 populated**,
  against the committed **3 keys, 3 populated**, dropping `[1,2,5,10]`
  (`P![-1, 4755/1024, -8267/1024, 797/128, -115/64]`).
  **The retracted explanation:** "the two runs used different target cover sets, so `08ce5fa` came
  from a path this checkout cannot reproduce." **That is false.** `{1,2,5,10}` (label 7584, g=1) is
  present in `GetHyperellipticCandidates()` AND in `Xstar`CoveredBy` today; the cover set is
  unchanged. Current code enumerates the cover and then withholds it **on purpose**.
  **What is actually happening:** commit `1768517` (2026-08-24 19:20) — *"Force-defer covers with an
  unpinned y2-scale instead of trusting their twist"*, i.e. the PR #38 / issue #36 guard — makes the
  run emit `Cover W={1,2,5,10} (g=1) has an unpinned y2-scale; deferring to recover as a quotient
  (twist untrusted)`, and back-fill then fails for want of a determined parent. The committed file
  is `08ce5fa`, 2026-08-24 **07:30** — **twelve hours older than the guard** — so its entry was
  produced by exactly the code path the guard was later written to distrust.
  ⇒ **But the committed entry is CORRECT, and that is measured, not assumed.** `VerifyModelSet`
  check [4] — the Eichler-Selberg trace-formula point count, independent of the Borcherds/Schofer
  machinery that generated the model — passes it (24 checks, 0 failures) and **discriminates the
  twist**: negative-controlled against six quadratic twists `d = -1, 2, -2, 5, -5, 11`, every one
  fails (3-5 failures each) while the committed curve passes. So this is a case where the guard is
  **CONSERVATIVE**: it suppresses a cover whose twist is in fact right.
  ⇒ **The lever this exposes.** The `y2`-scale is unpinnable from sparse CM data, but the twist is
  **independently decidable** by the point count already implemented in `ModelVerification.m`. So
  instead of dropping an unscaled cover, the pipeline could emit each candidate twist and SELECT
  the one matching `ComputePointsViaTrace`. Shipped behind `Y2TWIST=1` (`36ac71e`).

#### THE STALE-MODEL FINDING, and what `Y2TWIST` actually buys — measured 2026-09-05

**THREE committed model files do not regenerate from current code**, and all three for the same
reason: the y2 guard (`1768517`, 2026-08-24 19:20) POSTDATES them, so regeneration now withholds
covers they contain. Found by an 8-base regeneration sweep (6 identical, 2 differ) plus a
fresh-vs-fresh baseline run that proved the difference was staleness, not the vx fix.

    models_22_5.m   08ce5fa, 2026-08-24 07:30   12 h older than the guard
    models_22_3.m   the base issue #36 was actually filed about
    models_15_2.m

⚠ **They are NOT wrong.** All three pass `ModelChecks` and their Guo-Yang comparisons. They are
*unreproducible*. Nothing in CI regenerates a model and compares, which is why this went unseen —
now covered by **`tests/_offline/ModelRegen.m`** (`afa0412`), auto-discovering, no per-base
authoring, with these three recorded in `MR_KNOWN_DRIFT`.

**What `Y2TWIST=1` restores, measured against the committed files:**

    22_5    FULLY   3/3 populated covers, coefficient-for-coefficient
    15_2    FULLY   12/12 keys; one cover ([1,2,3,6]) differs in presentation but IsIsomorphic
    22_3    13/14   missing [1,3,22,66], and [1,66] loses 1 of its 3 entries

⇒ **The residual gap is GENUS 0** — both `22_3` losses are conics, which `select_y2_twist` skips
(`X`g lt 1`, since `HyperellipticCurve` needs degree >= 3).
⚠⚠ **"Extend twist selection to conics" was recorded here as the next step. It is IMPOSSIBLE, and
the guard is CORRECT — do not attempt it.** *(measured 2026-09-05, before implementing)* Every
quadratic twist of a conic has **exactly `p+1` points over `F_p`**, for every `p`: probed
`d = 1,-1,2,-2,3,-3` at `p = 3..23` on `22_3`'s own `[1,66]` conic `y^2 = 4x^2+1`, and every row is
constant. The reason is structural — a smooth conic over a finite field always has a rational point
(Chevalley-Warning), so it is isomorphic to `P^1` and its count is `p+1` regardless of twist. The
trace formula returns `p+1` too, so **the discriminator is vacuous in genus 0**. Extending the
selector would make every candidate "match", give `#good ne 1`, and defer exactly as now — same
behaviour, more work.
⇒ **What would actually pin a conic's twist is a different mechanism.** Over `Q` conic twists are
separated by rational solubility / Hilbert symbols, not by reductions, and we have no independent
source of truth for that here (the trace formula is blind). The route the pipeline already has is
`backfill_deferred`: a genus-0 cover under a DETERMINED higher cover is pinned as its quotient. So
this is a **back-fill** problem, not a twist-selection one.
⚠ Also measured: at `22_3` the selector reported `0 twist(s) matched` for `W={1,2,3,6}` — a clean
negative, not an ambiguity. Worth understanding before widening the candidate set.

**Still 0 new Guo-Yang equations.** `W={1}` at `22_5` remains empty, and `14_3`'s six lost covers
are plain under-determination, which this lever does not touch. The value here is
**reproducibility of the model corpus**, not coverage.
* In both bases `models[[1]]` now exists as a key and is **empty**, so the full curve is genuinely
  not produced — consistent across every variant tried.

⚠ **How the "43" was confirmed, because the obvious grep gets 41.** Two rows write the label
without braces round `D` — `$X^6_0(17)$` and `$X^6_0(29)$` — so a pattern anchored on `X^{D}_0(N)`
silently drops exactly those two and returns a plausible 41. Cross-check with the equation cell
instead: `multirow{1}{*}{\text}` occurs 43 times in the table range, once per base. (Same
read-a-fragment-and-generalise family as the four extraction traps in `GuoYangEquations.m`'s
header.) Also note `6_17` and `6_29` appear ONLY in CM-value captions elsewhere — having a
`tests/X0_6_17.m` does not imply a published equation, and `15_1` has a test but is **not** a GY
equation base at all.

**The 9 remaining blockers, with the classification CORRECTED 2026-09-04:**
* `15_4` — **structural**, and the only one of its kind: `N = 4` is not squarefree, and the method
  boundary is literally `assert IsSquarefree(N)`. The sole non-squarefree `N` in Guo-Yang's list.
* `93_1`, `95_1`, `159_1` — **RE-CLASSIFIED**. Previously filed under the `IsSquarefree(N)`
  boundary, which **cannot** be right: all three have `N = 1`, which IS squarefree. Measured
  traceback at `93_1` is Magma's OWN assert:
  `GalFldFun.m:305  assert vx ge 0`, reached via `AbsEltseq` on a `q^-60` pole — the **vx class**.
  `genmodels.m` already fails fast on it via `vx_skip = {<95,1>,<115,1>,<123,1>,<129,1>}`, but
  **that list is incomplete — `93_1` belongs in it** (`159_1` under test). See
  [[assert-failed-is-squarefree-n]], now corrected. ⚠ `genmodels.m` cites "memory
  vx-laurent-n0-circular", which **does not exist** — dangling reference.
* `39_2` — ✅ **RECOVERED 2026-09-05. No longer a blocker.** It was never the NONINTEGRAL failure
  its record claimed; it was starved by the coprime CM filter. With `CMNONCOPRIME=1` it sees 24 CM
  points against demand 19 and builds cleanly (15 keys, **0 empty**), and its `W={1}` genus-7 curve
  is **`IsIsomorphic` to Guo-Yang's published equation** (verified 0.06 s; now pinned in
  `tests/GuoYangEquations.m`, 9 bases). `ModelChecks` also passes it (82 files, 8573 checks, 0
  failures), which is an independent check since it uses trace-formula point counts rather than the
  Borcherds/Schofer path that produced the model.
  ⚠ `data/models/models_39_2.m` does NOT regenerate by default — see its header, and its entry in
  `ModelRegen`'s `MR_KNOWN_DRIFT`. The flag is NOT known to be safe in general; what makes this file
  trustworthy is the independent oracle, not the flag.
⚠ **On the two re-classifications below: those records were CORRECT WHEN WRITTEN, not careless.**
The intervening fixes (Schofer cusp-0 isometry PR #18, per-coset `tau` `475e72b`, the y2 guard
`1768517`, the vx shift `d9b52d0`) moved where these bases fail. So the remedy is **re-measure a
base whose record predates a fix that could move its failure** — not "triage more carefully". This
is the same phenomenon as the stale committed models above, one level up: a model goes stale when a
GUARD postdates it, a triage record when a FIX postdates it, and both are invisible to any check
that reads the artifact instead of regenerating it.

* `69_1` — ⚠ **RE-CLASSIFIED 2026-09-05: an exponent OVERFLOW, not a non-rational value.** The old
  entry called it "non-rational value (`RationalNumber` failure), the embedding-selection class".
  It is indeed inside `RationalNumber`, but the actual error is
  `LogSum.m:142  ret := &*[Rationals() | p^(Integers()!s`log_coeffs[p]) : ...]` ->
  `Runtime error in '^': Argument 2 is too large`, i.e. a `log_coeffs` exponent so large that
  Magma refuses the power. That is a different defect from the `15_2` embedding-selection story,
  so **do not assume the `15_2` root cause applies**.
  ✅ **PINNED 2026-09-05.** `RationalNumber` now names the offending prime instead of letting
  Magma emit its opaque `'^': Argument 2 is too large` (`LogSum.m`). The actual failure is:

      6*Log3 + 826241926712017437948244622352640031335552334419770916034895634120110682322*Log23

  i.e. the coefficient at **23 has diverged** (~`8.26e74`) while `Log3`'s is an ordinary `6`.
  `D = 69 = 3*23`, so BOTH primes are **ramified** — this is a runaway at a RAMIFIED prime, not a
  level-prime effect, and NOT related to the `p | gcd(d,N)` work.
  ⚠ Note the sum is supported ONLY on `3` and `23`, both dividing `D`. Whatever diverges is
  concentrated on the ramified places.
  ✅ **ROOT-CAUSED 2026-09-05 — the overflow is a SYMPTOM, not the disease.** `Kappa0` is fine:
  its per-`m` coefficients at 23 are all small fractions (`-1/2, -2/3, -3/2, -1, -5/6, -1/3`)
  across every CM discriminant. The runaway enters through the OTHER factor of
  `log_coeffs += Coefficient(f,-m) * Kappa0(m)` — **the Borcherds forms themselves have
  astronomical principal parts**. Measured, all five keys:

      key -1:  valuation -62, 29 nonzero principal coefficients,
               MAX |coeff| = 1.17e73 / 3072
      key 10669: MAX |coeff| = 1.96e76 / 1536      (others 10^73..10^76 likewise)

  ⚠ And they are NON-INTEGRAL (denominators 512, 384, 1536, 3072), which puts `69_1` in the known
  **non-integral-forms class** ([[nonintegral-forms-root-cause]]) — NOT the `p | gcd(d,N)`
  level-prime story, and not the embedding-selection story its original record named. Three
  different diagnoses, three different classes; only this one is measured.
  ⇒ NEXT: `IntegralSolution` is the documented rescue for that class (it "rescues 7 of 18"), so
  `INTSOL=1` is the obvious experiment.
  ⚠ **Attempted 2026-09-05 LOCALLY and the result is INCONCLUSIVE, not negative.** It died with
  `Bus error` after 1013 s at `Memory usage: 11335.62MB` — which is exactly the documented
  "Magma dies ~11 GB" ceiling on this Mac ([[tshift-ladder-gap-filler]]), and the machine had
  ~2.4 GB free at the time. **That is a machine limit, not a verdict on `INTSOL`.** Re-launched on
  lovelace (2 TB, 1.8 TB free); do not record `INTSOL` as refuted for `69_1` until that returns.
  ⚠ General lesson: `69_1` is in the class where the forms themselves are ~`1e75`, so ANY variant
  of it is memory-hungry. Run this base on lovelace, not locally.
  Also observed in the same run: the y2-scale `IsSquare` check fails for a cover at `69_1`, so
  this base additionally exercises the unpinned-y2 path.
* `111_1`, `119_1` — the odd-`D` basis ceiling; killed at **17.5 h CPU each**, `119_1` peaking at
  40 GB, inside `BorcherdsForms`.
  ⚠ **"Ceiling" may be the wrong word — measured 2026-09-05 with the new `BFPROGRESS=1`.** On
  `95_1` (same class: odd `D`, level 1, previously undiagnosable) the `m`-search costs
  **4653 s = 77.6 min for the FIRST `m` alone**, of **6**, and `m` grows in magnitude
  (`-7 -> -35 -> ...`) so later values should cost more. That projects to **8+ hours for one
  base**, and it is PROGRESSING, not wedged: `/proc` shows 100% CPU, no `normaliz` child, no I/O,
  flat 0.6 GB.
  ⇒ So 17.5 h may simply not have been long enough, rather than evidencing a hard wall. Before
  writing these two off again, **re-run with `BFPROGRESS=1` and read `m_idx=k of N`** — that
  distinguishes "slow" from "stuck", which nothing previously could. (`119_1`'s 40 GB is a
  separate and more serious signal; the memory-side claim is not affected by this.)
  ⚠ Note the earlier diagnosis was made from `/proc` because the logs were empty: `BorcherdsForms`
  reports each `m` at `vprintf` LEVEL 2 while `genmodels.m` defaults to `verb := 1`, and even at
  `VERB:=2` those go to buffered stdout and are lost on a kill. `BFPROGRESS` uses `WriteStderr`.
* `26_3` — **the anomaly is EXPLAINED as of 2026-09-05: an `s` <-> `s~` SWAP.** At discs `-267`
  and `-708` Guo-Yang's `s` sits in **our `s~` row** (`sIndex = 7` is `s`, row 9 is `s~`); the
  other 12 of 14 discs are correct. `s + s~ = 1` holds at EVERY disc (`8/25 + 17/25`,
  `11/49 + 38/49`, ...), and the exact `z -> z/(z-1)` is simply how an `s -> 1-s` swap looks after
  the checker's cross-ratio normalisation — the involution was the shadow, not the cause. The
  Mobius map from GY's coordinate to ours is `phi(w) = (1-w)/2`, verified on all nine
  non-reference rows.
  ⚠ **It is NOT a CM-point selection ambiguity** (the old description). Both values are the same
  point, different hauptmodul; and our table has each disc exactly ONCE (multiplicity checked).
  **Root cause:** the pair is pinned by `s + s~ = 1` (`find_signs_hauptmodul`), and that constraint
  is **symmetric under exchanging `s` and `s~`** — `8/25 + 17/25 = 1` reads the same either way —
  so the relation meant to resolve the ordering structurally cannot. The signs themselves are
  forced (only `+ +` sums to 1); the freedom is purely in the LABELLING.
  ⇒ **Before fixing: identify what pins the ordering at the other 12 discs.** Something is doing
  real work there, and adding a tie-break without knowing the intended invariant is guessing.
  ⚠⚠ **THE SWAP IS *NOT* WHY `26_3` HAS NO MODEL — measured 2026-09-05, and this closes the
  question flagged as open here.** Re-run against current code, `26_3` dies much earlier:
  `Computing absolute values at CM points...Runtime error: Could not find enough points, sorry!`
  It is **CM-STARVED**, confirmed by `cmsupply`:

      BASE 26 3 demand 9 genera [ 0, 0, 0, 1, 1, 1, 2 ]
      POOL rat 3 quad 0 include 3
      CMVERD 26 3 SHORT margin -5

  ⚠⚠ **BUT "26_3 IS CM-STARVED" IS WRONG — RETRACTED SAME DAY. It is the COPRIME FILTER, and the
  swap is ON the critical path after all.** The `SHORT` verdict measures the FILTERED pool, not the
  available points. Measured:

      bd | coprime_to_level | #rat | #quad | total      (demand 9)
       2 | true             |   3  |   0   |   3
       2 | false            |   7  |   0   |   7
       4 | true             |   3  |   0   |   3
       4 | false            |  14  |   7   |  21   <-- at the DEFAULT bd
       8 | false            |  15  |  22   |  37

  At the default `bd := 4`, dropping `coprime_to_level` gives **21 points against demand 9** — and
  those 14 rational discs are exactly Guo-Yang's 14-row table. `26_3` has ample CM points.
  **The causal chain runs: non-coprime discs misbehave -> the filter excludes them
  (`SchoferFormula.m:1063` passes `coprime_to_level := true`) -> pool 21 -> 3 -> "Could not find
  enough points".** So the earlier claim here that "fixing the swap would NOT unblock the model"
  was backwards: the swap-class misbehaviour is precisely WHY the filter exists.
  **And the swap is confined to that class.** Of Guo-Yang's 14 discs only `-8, -11, -20` are coprime
  to `N = 3` — exactly the 3 the filter admits — and BOTH swapped discs (`-267`, `-708`) are
  non-coprime, while 9 other non-coprime discs give correct values. The filter is blunt: it drops
  11 points to avoid 2 bad ones (`ShimuraQuotients.m:1420` says so in as many words).
  ⇒ **Route to unblocking `26_3`: fix the non-coprime misbehaviour, then relax the filter** (or
  extend the `Keep` exemption). ⚠ NOT established: that the swap is the ONLY misbehaviour, or that
  relaxing the filter yields a CORRECT model — the bad points would poison the solve, which is what
  the filter is protecting against. Verify against Guo-Yang's table before trusting any model so
  produced.
  ⚠ Note `BorcherdsForms.m:709` ALREADY falls back to `coprime_to_level := false` for CM-starved
  bases, but `AbsoluteValuesAtCMPoints` does not — an asymmetry worth understanding.
  Write-up: memory `26-3-hauptmodul-swap`.

**⚠ The lesson from `51_1`/`57_1`: "no recorded failure" was being read as "blocked".** Neither had
ANY triage record — no `INTSOL`, no `BASEVERD`, nothing. They had simply never been run. One
`genmodels` run each produced 4 cover-keys apiece, and `51_1` reproduces Guo-Yang's
`y^2 = -(x^2+3)(243x^6+235x^4-31x^2+1)` **exactly** under `x -> 3x, y -> (27/16)y` (over-determined
and consistent across all five coefficients). **Before filing a base as blocked, check whether it
was ever attempted.**

### The two-tier test plan (adopted 2026-09-04)

- [~] **TIER 1 — ATTEMPTED 2026-09-04 AND ABANDONED. The "cheap" framing was WRONG.**
      The Magma comparison is indeed cheap; **acquiring the data reliably is not**, and that was
      the whole cost. Automated extraction of Guo-Yang's equation tables from the arXiv LaTeX hit
      FOUR distinct silent-corruption bugs, every one of which yields a plausible wrong polynomial
      rather than an error:
      1. an equation may **wrap across `\\` into several `$...$` groups** (6_11 loses four terms
         if you split on `\\`) — the same content-split-with-no-marker class as the Magma
         line-wrapping trap, third instance this session;
      2. a leading `-` may belong to the first TERM, not be an overall factor (26_1 came out
         sign-flipped);
      3. the base label `$X^{D}_0(N)$` is itself inside `$...$`, so segmenting AT it flips the
         parity of every subsequent `$` pairing;
      4. newlines inside the math break term parsing.
      **And the target was wrong twice over.** Guo-Yang's tables are heterogeneous: some rows are a
      single `y^2=f(x)`, others a PAIR (`82_1`: `y^2=f(s)` AND `x^2=g(s)`), `15_4` is a conic in
      `z`, and `93_1` mixes two variables (`3s^3-7s^2-3t-1`) — almost certainly a typo in the
      paper.       ⚠ **One claim made here on 2026-09-04 was WRONG and is retracted:** that the published
      equation does not uniformly correspond to a key of our model files, "evidenced" by
      `X0_26_1.m` comparing a degree-3 `s`-equation. That came from grepping only the FIRST
      `HyperellipticCurve` in the file. `X0_26_1.m` actually carries THREE cover entries and its
      `{1}` entry IS Guo-Yang's degree-6 equation, with the identity matrix. **The mapping to
      `W={1}` is uniform as far as checked**, and the apparent `82_1` genus contradiction was the
      parser taking only the first of a PAIR, not a mapping failure.
      ⇒ **Do not retry the automated extraction.** The existing 24 tests were hand-transcribed, and
      that was the right call. Scripts kept for reference at
      `vvdata/weyl-campaign/guoyang/extract_equations.py` (campaign) with the traps documented in
      its header; the generated `gy_equations.m` is **NOT trustworthy** and was removed from `main`.
- [x] **TIER 1' (replacement) — hand-transcribe the ~10 bases that matter. DONE 2026-09-05 for all
      that are transcribable: 8 of 10**, in `tests/GuoYangEquations.m` (`ed38ed5`, `41da3fb`,
      `64d9316`). `14_3` and `22_5` are NOT transcribable — see the stock-take above; they need
      models generated. The hand route worked exactly as predicted: slower per base than a parser,
      and the only approach with a working track record here.
      **Method that made each entry trustworthy, worth reusing:**
      * **Read the raw `.tex` row, not a parse of it**, and read the neighbouring rows too — `51_1`
        and `55_1` sit either side of `57_1` and are already-passing hand-transcribed entries, so
        matching them byte-for-byte calibrates the reading before the new row is trusted. This is
        the "reproduce a KNOWN value before trusting a new one" habit applied to transcription.
      * **Negative-control every entry before committing.** For `57_1`, five single-coefficient
        perturbations of the published pair, each chosen to KEEP genus 3 so that `IsIsomorphic`
        does the discriminating rather than the genus assert — all five correctly rejected. A
        perturbation that changes the genus proves nothing about the comparison.
      * `57_1` cost ~13 s; `21_2` still dominates the file at ~100 s (112 s total, 8 bases).
- [ ] ~~TIER 1 (superseded)~~ **one cheap test, all 34 bases at once** `ModelChecks` validates the
      79 model files structurally (genus, Weil-polynomial divisibility, point counts) but **never
      compares them to Guo-Yang**. So nothing in CI would notice a committed model silently
      disagreeing with the paper. A test that checks the **STORED** model against the **published**
      equation, up to the hyperelliptic isomorphism, needs **no pipeline run** — milliseconds — and
      covers every base where both exist. Best value per unit of CI time available here.
- [ ] **TIER 2 — pipeline tests for the FAST bases only.** The existing `X0_D_N.m` pattern calls
      `AllEquationsAboveCovers` and re-derives the curve, which is the stronger test (it exercises
      the computation, not the stored data). Measured: `51_1` takes **240 s** — comfortably CI-able.
      But `87_1` is likely hours (its CM-value run was ~7000 s) and GitHub caps a job at 6 h, so
      slow bases go to `tests/_offline/`, reusing the split already set up for the CM-value tables.
      ⚠ **Per-base cost is NOT just CI time.** Each `X0_D_N.m` needs an **isomorphism matrix**
      between our model and the published one — that is **not in the paper** and must be derived
      per base (`51_1`'s is `x -> 3x, y -> (27/16)y`) — plus the Atkin-Lehner involutions as
      matrices in our coordinates. Budget 30-60 min of careful work each, and note that a
      plausible-but-wrong matrix would make the test **vacuously pass**. Derive it the way `51_1`
      was done: solve the transformation from two coefficients and CHECK the rest agree.

---

## HOUSEKEEPING — about an hour, once

- [x] **Delete the 12 merged local branches.** *(done 2026-09-02)* Local 23 -> 11. Used `-d`, so
      git would have refused anything not genuinely merged. Was: `even-control`, `even-correction`,
      `fix-pointless-conics-empty`, `intsol-optin`, `m0-exact-oracle`,
      `m0exact-relative-tolerance`, `odd-d-invariant-hoist`, `odd-d-zeroskip`,
      `preprint-assembly-theorem`, `preprint-two-channel`, `tier1-model-data`,
      `y2-unscaled-deferral`.
- [x] **Delete the merged remote branches.** *(done 2026-09-02)* Remote 31 -> 13. This split was
      wrong in the original plan: there were **17**, not 7 — the other 10 were the remotes of the
      local branches above, which only became remote-only once those were deleted. Was: `add-external-cm-value-tests`,
      `add-kudla-yang-local-tests`, `add-m0-local-density-tests`, `fix-cm-supply-divisor-discs`,
      `fix-wpoly2-p2-local-density`, `m0-vv-constant-term`, `preprint-n5-validation`.
- [x] **Retire or tag the three stale ones.** *(done; confirmed 2026-09-04)* `non_optimal`,
      `odd_DN` and `pointlessconics` are gone as branches and preserved as `archive/<name>` tags,
      exactly as recommended — **all 9 archive tags exist on `origin`, not just locally**.
- [x] **The six stale local branches: DONE, and nothing was lost.** *(confirmed 2026-09-04)*
      `m0-hauptmodul-ground-truth`, `m0-prop53-anisotropic`, `m0-level-rule-e2e`,
      `fix-m0-multiterm`, `valuesatcmpoints-diagnostic` and `fix-15-2-find-signs` are all deleted
      AND preserved as `archive/<name>` tags on `origin`. That matters for `fix-15-2-find-signs`
      in particular, which had 3 commits that never reached its own remote — the tag holds them.
      `whbasis-speedup` is deliberately kept as the CI-green record for a cherry-picked commit;
      it is still there, correctly.
- [x] **The `polymake/nmzsolve.err` stray is gone.** *(confirmed 2026-09-04: campaign has 0 untracked files)*
- [x] **Retire the merged worktrees.** *(done 2026-09-02)* `-fix`, `-hoist` and `-control` removed;
      worktrees 6 -> 3 (`tier1-models`, `-campaign`, `-mainport`). `-control` needed `--force` for
      its 23 untracked files, all previously verified as duplicates; its one single-copy item was
      carried over as `f3fcf1e` on the campaign branch.
- [x] **Resolve the `deficit.m` contradiction.** *(done 2026-09-02)* `HANDOFF.md` was right and
      memory was stale. Campaign `0bc7f28` already made the fix; the file has no `Probe*` reference
      and calls only public intrinsics, so the three validated results stand unpatched. Memory
      corrected. Still **even-D only**.
- [x] **Pull the lovelace clone up to `main`.** *(done 2026-09-04)* It was **52 behind**, not
      current — a plain `git rev-list HEAD..origin/main` there reports 0 because its `origin` ref
      is stale, so **fetch first**. Now at `05471c8`, 0 behind, carrying the tau fix and the
      t-shift `nmzsolve.py`. What was discarded, and why it was safe:
      * `EquationsCovers.m` — was byte-identical to `origin/main`; the pointless-conics guard
        really was redundant, as this item originally claimed.
      * `VectorValuedForm.m` — was NOT that guard but the **GATE3B report-and-continue
        diagnostic**. Preserved on campaign as `vvdata/weyl-campaign/gate4/gate3b-gate4-probes.patch`
        (diffed to confirm), so nothing was lost — but note the pull DID remove GATE3B-measuring
        capability from that machine. `git apply` the patch there to get it back.
      * its 109 untracked solves — checksummed against the now-committed copies before deleting,
        and restored by the pull (907 polymake files there now).
      ⚠ Its five root `.m` triage scripts stay untracked by design (that is how they are run
      there). Two were STALE and have been refreshed from campaign: `deficit.m` still called
      `ProbeWHBasis`/`ProbeDivisorMatrix`, the instrumentation-only functions campaign's `0bc7f28`
      removed so the file runs on an UNPATCHED tree — so it would have failed there; and
      `genmodels.m` was missing the INTSOL toggle. `cmsupply.m`, `ppint.m`, `gencapped.m` already
      matched. **Re-check these against campaign whenever you use that machine** — they are copies,
      and nothing keeps them in sync.
- [x] **Harvest the 109 stranded Normaliz solves from lovelace.** *(done 2026-09-04, `e6d8ff2`)*
      They existed on lovelace only, committed nowhere, `M = 708..3372` — past the old frontier.
      Cache 797 -> 906 files. Verified before trusting: not t-shift output (lovelace runs the
      99-line `nmzsolve.py`, no `tshift` files, so the fallback could not have fired), all
      well-formed, and the 16 files lovelace shared with `main` agree as **vector sets** (byte
      comparison is invalid — committed files use escaped `\[` line starts).
- [x] **`tau-precision` probes are committed.** *(confirmed 2026-09-04)* Five files at
      `vvdata/weyl-campaign/tau-precision/` on the campaign branch since `c672464`.
