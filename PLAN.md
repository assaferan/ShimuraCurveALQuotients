# Working plan — what to do next

Assembled 2026-09-02 from the branch audit, `HANDOFF.md` on `main`, and the GATE3B measurement on
`10_61`. `HANDOFF.md` remains the record of *what happened*; this file is the record of *what to
do*. When they disagree about state, `HANDOFF.md` wins.

Five tracks. One is the main line; the rest run in parallel and **none of them blocks it**.

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
- [ ] **Work the intersection of the two propositions you already have.** `prop:kappa0` gives
      `m = 0` at a nonzero isotropic coset; `prop:closedcoef` gives all `m` for the scalar.
      `A_m` is the missing corner of that square.
- [ ] **Hold the success criterion: all 13, or it isn't the theorem.**

**Why this one.** It *is* the `Kappa0` log-`N` defect at `SchoferFormula.m:589` — the same defect
`coprime_to_level` (`ShimuraQuotients.m:1420`) papers over with what the code itself calls "a blunt
instrument". Closing it unblocks the 49 obstructed bases in a single move.

**Honest risk.** No product of local densities reproduces `b` under any sign convention, with KY
Props 5.3/5.4 validated at 470+150 checks. This may need new mathematics rather than a derivation
from what is in hand.

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
- [ ] **Retire or tag the three stale ones:** `non_optimal` (Oct 2025, 429 behind), `odd_DN`,
      `pointlessconics`.
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
