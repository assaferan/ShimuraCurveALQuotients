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

- [ ] **Kill the `10_61` run** (PID 2803742 on lovelace). The result is in; it is 10 h deep, RSS
      climbing 940 MB → 3.6 GB, and now grinding with poisoned `a0tab` entries, so anything it
      writes is untrustworthy.
- [ ] **Find what distinguishes `wi = 2` and `wi = 1221`.** Exactly two cosets fail, each at *all*
      462 discriminants (924 = 2 x 462), and they fail in *opposite* directions: `k0` overflows to
      `+inf` at one and collapses to `~1e-10` at the other, while `k1` stays sane
      (1.1e4 – 2.4e12) at both. Chase that reciprocal symmetry — an exponent sign flip, or a coset
      and its negative.
- [ ] **Cross-check against the `39_2` malformed-form pathology.** Same family — a non-integral
      principal part poisoning `M0MultiplierExact` by catastrophic cancellation — but far more
      extreme: there the magnitude was 0.026, here it overflows the real field.
- [ ] **Re-run GATE3B on `14_43` once the cause is understood.** Its log ends identically to
      `10_61`'s, so expect the same defect; one root cause probably covers both.
- [ ] **Reclassify both bases in the sweep record.** `10_61` was one of the "two runnable
      candidates of 122". It is not runnable — it has a real upstream defect, not a threshold
      problem.

Raw evidence: `~/shimura/models/10_61.gate3b.log` on lovelace (353 KB, 926 lines).

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

- [ ] **Delete the 12 merged local branches:** `even-control`, `even-correction`,
      `fix-pointless-conics-empty`, `intsol-optin`, `m0-exact-oracle`,
      `m0exact-relative-tolerance`, `odd-d-invariant-hoist`, `odd-d-zeroskip`,
      `preprint-assembly-theorem`, `preprint-two-channel`, `tier1-model-data`,
      `y2-unscaled-deferral`.
- [ ] **Delete the 7 merged remote-only branches:** `add-external-cm-value-tests`,
      `add-kudla-yang-local-tests`, `add-m0-local-density-tests`, `fix-cm-supply-divisor-discs`,
      `fix-wpoly2-p2-local-density`, `m0-vv-constant-term`, `preprint-n5-validation`.
- [ ] **Retire or tag the three stale ones:** `non_optimal` (Oct 2025, 429 behind), `odd_DN`,
      `pointlessconics`.
- [ ] **Clean the `-control` worktree:** `git clean -fd && git checkout VectorValuedForm.m`.
      Everything in it is provably duplicated elsewhere; the one single-copy item was carried over
      as `f3fcf1e` on the campaign branch.
- [ ] **Resolve the `deficit.m` contradiction.** Memory says broken, `HANDOFF.md` says repaired and
      validated. The handoff is newer and authoritative — but memory is what a fresh session reads
      first, so fix whichever is wrong.
- [ ] **Pull the lovelace clone up to `main`.** Its hand-applied pointless-conics guard is now
      upstream as `127b044`, byte-identical, so the local modification is redundant.

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
