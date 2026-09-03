# Working plan — what to do next

Assembled 2026-09-02 from the branch audit, `HANDOFF.md` on `main`, and the GATE3B measurement on
`10_61`. `HANDOFF.md` remains the record of *what happened*; this file is the record of *what to
do*. When they disagree about state, `HANDOFF.md` wins.

Five tracks. One is the main line; the rest run in parallel and **none of them blocks it**.

## Picking this up cold

As of 2026-09-03 (evening) everything is committed and pushed, no worktree is dirty, and nothing
is running. Three items are ready to start immediately, in decreasing order of value:

1. **The `A_m` theorem** (MAIN LINE). Three more routes were closed today (KY Prop 5.4/5.5
   insertion, Schwagenscheidt's oldform relation, the `s`-law/genus-theta closed form) — see that
   section for what NOT to re-attempt. What's actually open now is narrower than "state and prove
   a theorem": the `rem:gauge` ambiguity between `-a_E` and Table A is *provably resolvable* (full
   rank 6 on the 158-monomial data), and the concrete blocker is a **missing data dump** —
   constant-term data at 4 intermediate cusp classes (`g = 2,6,10,30`) that a `cusp7.m`-style pass
   never recorded. Get that, and the reachable-subspace solve for `A_m` should go through. This is
   the item that unblocks 49 bases — the other two are tidy-up by comparison.
2. **Per-coset `tau` in `M0MultiplierExact`** (REPAIR). Fully specified and needs no theory: target
   `Im(w*tau)` in `[1e-3, 1e-2)`, gate on `15_2` 9/9 and `58_5` `ModelChecks` 48/0. Self-contained.
3. **Re-measure `wi = 1221`** (REPAIR). Small, and the last gap in an otherwise closed result.
   Needs a real `DedekindEta` with full `SL2` reduction — do **not** reuse `tauwindow.m`'s
   one-`S`-step `logeta`, which is invalid in exactly that regime (see its README).

Worktrees: `tier1-models` (here, `.`), `worktrees/campaign` (`m0-theta-campaign`, research data),
`worktrees/mainport` (`main`) — nested under this checkout, not siblings; see CLAUDE.md for the
convention on new worktrees.
Probes live at `vvdata/weyl-campaign/tau-precision/` on the campaign branch.

## State of play

     4   verified models, from one base (58_5, partial set by construction)
    49   Borcherds-obstructed bases
    81   of 122 still unclassified
   924   gate-3 failures at 10_61, zero near-misses
    30   pages of paper; the theory arc is closed

Two arcs running at different speeds. The theory arc is essentially closed — the multiplicity
formula end to end, `kappa_mu(0)` derived, N>1 subsuming N=1, the support rule correctly reframed
as a gauge. The model arc is blocked, and not by engineering.

Across the record, most proposed levers were refuted once measured: the 336x hoist worth 2%, the
odd-D eta-quotient explosion that wasn't, three successive symbol laws, and now the scalar `a_E`
at 1 of 13. That is discipline, not failure — but the cumulative message is that the engineering
levers are spent and the binding constraint is mathematical. Plan around that.

---

## MAIN LINE — the `A_m` theorem

One object blocks disproportionately much. Everything below this section is secondary to it.

- [ ] **State it precisely before attempting anything.** The vector-valued weight-3/2 Eisenstein
      coefficient `b^{eta*}_0(r)` at a *nonzero isotropic coset*, general `m`, supported on
      `N | r` — with `A_r = -b(r)/4`. A different object from the scalar, not a special case.
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
      **Concrete next step:** re-run (or extend) a `cusp7.m`-style dump so every intermediate cusp
      class gets its constant term recorded — not just the first one encountered while iterating
      cosets — then redo the reachable-subspace solve. That is a scoped implementation task, not
      another open-ended search.

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
- [ ] **Re-read `conj:sq` under the gauge finding.** It is stated in terms of `w_square`, and
      `w_square` is a gauge — a choice of representative, not arithmetic. The paper was updated
      for that finding (`ee94175`); confirm the conjecture is still well-posed rather than a
      statement about a normalisation.
- [ ] **Merge `paper/` into `main`, or write down that the split is deliberate.** The ambiguity is
      what produced the 27/27 drift and cost a working branch the pointless-conics guard.
      ⚠ Merge trap: anything cut from `tier1-models` drags the whole paper rewrite with it — run
      `git log origin/main..origin/<branch>` first.
- [ ] **Rebuild the PDF, confirm refs resolve, submit.**

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
- [ ] **NOW THE OPEN HALF: why does `wi = 1221` fail?** The roles have swapped. Its measured
      spread is **1.39 digits — identical to typical cosets** (`wi` 300, 611, 900 all give
      1.39436), so the dynamic-range mechanism does *not* explain it.
      ⚠ **But that number is not trustworthy**: at `Im(z0) = 8.8e-7`, *both* `z` and `-1/z` have
      small imaginary part, so the one-`S`-step approximation used in `tauwindow.m`'s `logeta` is
      invalid there. **Re-measure with a real `DedekindEta` (full `SL2` reduction) before
      concluding anything about this coset.**
- [x] ~~Cross-check against the `39_2` malformed-form pathology.~~ **Superseded.** This is not a
      malformed form. It is an evaluation-point failure of the harness — see below.
- [ ] **Re-run GATE3B on `14_43` once the cause is understood.** Its log ends identically to
      `10_61`'s, so expect the same defect; one root cause probably covers both.
- [ ] **Reclassify both bases in the sweep record.** `10_61` was one of the "two runnable
      candidates of 122". It is not runnable — it has a real upstream defect, not a threshold
      problem.

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

- [ ] **Implement per-coset `tau` selection** in `M0MultiplierExact`, targeting
      `Im(w*tau)` in `[1e-3, 1e-2)`. Validate on `15_2` (must stay exact, `tests/M0MultiplierExact.m`
      9/9) and `58_5` (must keep its four verified models, `ModelChecks` 48/0) before trying
      `34_11` / `10_61`.
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

Probes: `vvdata/weyl-campaign/tau-precision/` on the campaign branch — `tauprobe.m` (per-coset
`Im`), `tauscale.m` (the `M^2` law), `tauwindow.m` (the dynamic range and the band table), plus
`amtest.m` (the scalar-`a_E` comparison). ⚠ `tauwindow.m`'s `logeta` is a **one-`S`-step
approximation**: valid when `Im(z)` is large or one `S` makes it large, INVALID when both `z` and
`-1/z` sit near the real axis. That is exactly the `wi = 1221` regime, so its small-`Im` numbers
are provisional.

---

## MEASURE — route B's `k = 3/2` phase

The only cheap-predictor path left, and it converts 81 unknowns into data.

- [ ] **Fix the half-integral convention *from theory*, not by fitting.** The error is very nearly
      `-d/6`. Do not tune that constant — a fitted phase would silently corrupt all 81
      classifications.
- [ ] **Validate against the known deficits:** `38_5 -> 1`, `38_7 -> 0`, `34_3 -> 0`.
- [ ] **Classify the 81**, including whether the obstructed class is larger than 49.

The machinery is already sound: the trace formulas reproduce 11 classical `dim M_k(SL2(Z))` values
exactly. Only the half-integral phase is wrong.

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
- [ ] **Retire or tag the three stale ones:** `non_optimal` (Oct 2025, 429 behind main), `odd_DN`
      (279 behind, 12 ahead), `pointlessconics` (277 behind, 2 ahead). Still open as of 2026-09-03;
      `odd_DN`'s non-mergeability was already investigated (memory: `odd-dn-branch-not-mergeable`).
      Recommend: `git tag archive/<name> origin/<name>` for each (preserves the commits under a
      permanent pointer), then delete the branch, local and remote.
- [ ] **NEW 2026-09-03: six more local branches are stale and unmerged into either `tier1-models`
      or `main`**, all predating the 09-02 cleanup and missed by it (`git branch -d` correctly
      refused them — they're genuinely unmerged): `m0-hauptmodul-ground-truth`,
      `m0-prop53-anisotropic`, `m0-level-rule-e2e`, `fix-m0-multiterm`,
      `valuesatcmpoints-diagnostic` (all 2026-08-14 to 08-21, commit messages mostly prefixed
      `SCRATCH:`), and `fix-15-2-find-signs` (2026-08-19, and **has 3 commits not even on its own
      remote** — `origin/fix-15-2-find-signs` is 3 commits behind the local branch). All look like
      superseded intermediate stages of the m=0 investigation that shipped a different way (via
      `prop:kappa0` in the paper) — worth a skim before deleting, since `m0-level-rule-e2e` in
      particular touches the same KY Prop 5.3/5.4 territory as this session's CI test and might
      have relevant scratch notes. `whbasis-speedup` is *not* in this list — it's deliberately kept
      as the CI-green record for a cherry-picked commit (see HANDOFF.md); leave it.
- [ ] **`worktrees/campaign` has one untracked stray file:** `polymake/nmzsolve.err` (136 bytes,
      2026-08-29, Normaliz solver log debris). Harmless; safe to `rm` whenever convenient.
- [x] **Retire the merged worktrees.** *(done 2026-09-02)* `-fix`, `-hoist` and `-control` removed;
      worktrees 6 -> 3 (`tier1-models`, `-campaign`, `-mainport`). `-control` needed `--force` for
      its 23 untracked files, all previously verified as duplicates; its one single-copy item was
      carried over as `f3fcf1e` on the campaign branch.
- [x] **Resolve the `deficit.m` contradiction.** *(done 2026-09-02)* `HANDOFF.md` was right and
      memory was stale. Campaign `0bc7f28` already made the fix; the file has no `Probe*` reference
      and calls only public intrinsics, so the three validated results stand unpatched. Memory
      corrected. Still **even-D only**.
- [ ] **Pull the lovelace clone up to `main`.** Its hand-applied pointless-conics guard is now
      upstream as `127b044`, byte-identical, so the local modification is redundant.
- [ ] **Commit the `tau-precision` probes to the campaign branch.** Staged at
      `vvdata/weyl-campaign/tau-precision/`; they are currently scratchpad-only, which is the exact
      failure mode that lost the original `genmodels.m` to a nightly `/tmp` purge.

---

## DO NOT — each of these already cost a session

* **No third slash-constant tolerance change.** Settled twice over: 924 failures, every one at
  `reldiff > 0.999`, zero near-misses. No threshold in (0,1) passes any of them.
* **No deeper pole orders.** The deficit is exactly 1 and *invariant* under enlargement; refuted
  at `38_5` (bump 0→8: rows 164→172, cols 36→38, rank 35→37).
* **No re-attempting the even correction as an implementation task.** It needs the theorem.
* **No tuning route B's phase constant to fit.**
* **No more odd-D constant-factor work.** Exhausted; the ceiling is
  `basis_of_weakly_holomorphic_forms(... : Zero)`.
* **Do not treat `ppint`/`cmsupply` as cheap triage on large bases** — 81 of 122 gave no verdict
  in a full hour.

## TRAPS — the ones that keep biting

* **Tilde plus a remote path.** It bit the same run twice in opposite directions:
  `OUTDIR:=~/shimura/models` did *not* expand (bash only expands after `:` in real assignments), so
  `Write` would have died on a literal `~` path; while `LOG=~/...` in a watcher *did* expand — to
  the local Mac home, which was then sent to a Linux box. Single-quote remote paths or use `\$HOME`.
* **`magma | tail` hangs** — an error drops it to the interactive prompt where it blocks on stdin.
  Redirect `< /dev/null`. One such hang burned 4 h 23 m at 0.02 s CPU.
* **Magma buffers stdout to a file.** An unchanged log is *not* evidence of no output — the
  `10_61` log sat at 154 bytes, then flushed 353 KB at once.
* **`pkill -x magma` matches nothing** (the binary is `magma.exe`). Verify a kill with a
  *different* pattern than the one used to kill.
* **`git log --all -- file.m` misses `vvdata/weyl-campaign/file.m`.** Search by basename:
  `'*file.m'`. Triage tooling lives on the campaign branch, never at the root.
* **No `timeout` on the Mac**; lovelace has GNU `timeout` and it truncates cleanly.
