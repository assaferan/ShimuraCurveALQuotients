# Handoff — 2026-09-23

**The newest section is the FIRST one below; everything after it is older and kept for provenance.**
Earlier material still says things like "34 of 43", "23 of 34 tests check involutions", or that
`95_1` is the only killed base with an oracle — those counts are STALE.

✅ **2026-09-13: everything is COMMITTED AND PUSHED on both branches**, and the branch-divergence
invariant prints nothing against `origin`. ⚠ lava's clone is still stale at `8dac84c` — `git fetch
&& git reset --hard origin/main` there before any new run, but NOT while a job is alive.

**➡ For what to do next, see `PLAN.md`.** This file records *what happened*; when the two disagree
about state, this file wins.

## Handoff — 2026-09-23 — THE ORACLE AUDIT: a wrong cohort, a base nobody checked, and a guard that could not fail

**Theme, because all five items share it: every defect found today was a WRONG OBJECT, not a wrong
computation.** A cohort labelled by the wrong criterion, a coverage question asked of the wrong
screen, a ratio computed against the wrong denominator (twice), a guard counting the wrong thing,
and an instrument answering about the wrong ambient. Nothing arithmetical was broken anywhere.

### ✅ `95_1`, `119_1`, `159_1` ALL have published Guo-Yang equations — `PLAN` said only `95_1` did

`PLAN.md` called `95_1 115_1 123_1 119_1 159_1` "the five Guo-Yang bases" and said `95_1` was the
only one carrying an oracle. Both wrong. The five are the jobs `earlyoom` reaped — a KILL cohort.
The Guo-Yang label came from `genmodels.m`'s `vx_skip = {95_1,115_1,123_1,129_1}`, which groups by
the **vx defect**, and nobody checked it against the paper. Against the 43 equation cells:

    95_1  119_1  159_1   published equation -> RETURNS WITH AN ORACLE
    115_1 123_1           NOT in Guo-Yang at all

⇒ the ~400 CPU-hour loss cost **three** oracle-bearing bases, not one. All three are now
transcribed in `tests/GuoYangEquations.m` (`31bd605`), read three independent ways (journal page,
journal PDF text layer, arXiv v1 TeX, all agreeing), reporting `PENDING` until a model lands and
then comparing automatically. Each is checked WITHOUT a model via the involutions the same table
publishes, plus the genus column; 6 symmetry controls keep that non-vacuous.
⚠ **Both checks are needed.** Truncating at the `\\` wrap — the documented trap — leaves all three
invariant under their own involution and is caught ONLY by genus; a mistyped middle coefficient
preserves degree and is caught ONLY by the involution.

### ✅ `69_1` closed — it had a published equation and NOTHING external checking it (`11532a3`)

Neither an entry in `GuoYangEquations.m` nor any `X0_` test; only `ModelChecks`, which is
structural and passes on a wrong curve of the right genus. Cause was structural, not an oversight:
`models_69_1.m` landed 2026-09-14 (`4bfb859`), AFTER both `X0_*` batches (2025-11-19, 2026-09-06/07)
and after the table was last extended. **A model that arrives after a sweep is never swept.**
Now has both halves: `IsIsomorphic` against their degree-8 curve (2 perturbation controls fail
correctly), plus `tests/X0_69_1.m` re-deriving all four covers with `w_3`/`w_69`. Negative-controlled
(swapping the matrices goes red, 2 involutions compared, 5 torsor maps tried). **117 s, so it is in
`tests/` and CI-visible** — the only one of `69_1/87_1/39_2/111_1/93_1` that CI sees.

### ✅ `tests/OracleCoverage.m` — so this cannot recur (`f0c3c67`)

    ok (40 of 43 published bases have a model; 12 via the equation table, 39 via an X0_ test;
        1 exempt; 3 not yet built: 119_1 159_1 95_1)          0.03 s, no pipeline run

Fails when a Guo-Yang base has a model and no oracle. Exemptions must carry a reason (only `15_4`,
Remark 39), and a STALE exemption is itself a failure, so the list cannot become a hiding place.
⚠⚠ **Its third negative control FAILED against the first version, and that is the real lesson.**
With the search pattern deliberately broken the test still printed `ok`: `93_1`'s hardcoded special
case held the count at 1, so the `eq 0` non-vacuity guard could not fire. **A guard incapable of
failing, inside the file whose whole purpose is catching that.** Reading the code it looks correct;
only running the control exposed it. Fixed by counting pattern hits separately from the special
case, plus a `51_1` canary.

### ✅ PLAN item 2 CLOSED — the coprime fix cannot move any recorded verdict (`c7108fb`)

Not "false clears", not "false obstructions" — **the question was asked of the wrong screen.**
`screened-2026-09-14.txt` was produced by `deficit.m`, and `deficit.m` says of itself that it
computes the deficit *"WITHOUT the CM points ... therefore skips `RationalandQuadraticCMPoints`"* —
exactly the call `0ca6e37` changed. Skipping it is what makes it a fast predictor. ⇒ **none of the
99 recorded verdicts (77 at `N>1`) can move.** All 99 are even `D`, consistent with `deficit.m`
being even-`D` only.
Measured anyway on `deficit_odd.m`, which DOES consume the pool, via the `PTSCOPRIME=1` control:

    38_5  tgt 4->9  wdef 1 every rung  obstructed->obstructed    34_11  tgt 3->5  wdef 0  clear->clear
    15_2  ladder numerically IDENTICAL both ways                 134_3  identical CM-supply error both ways

At `38_5` **only `tgt` moves**; `rows cols nds rank deficit wdef` are identical across all five
rungs — the extra targets land INSIDE the image, so `dim W` and `dim(W meet Im)` grow together.
⚠ Side finding for the screen, not this item: **`15_2` is a SECOND odd-`D` counterexample** — it has
a committed model and BUILDS, yet exhausts its `all_ms` ladder at `wdef >= 1` and prints
`obstructed`, identically in both modes. `deficit_odd.m`'s header names only `21_2`.

### 🔄 The involution gaps — all four now have matrices, NONE yet verified end to end

`87_1` and `111_1` by transcription, `39_2` and `93_1` by work:

    87_1   w_3,w_87    diagonal; our stored poly is EVEN in x, so the change from GY is diagonal
    111_1  w_37,w_111  cover_data IS their curve, so their formulas apply verbatim
    39_2   w_2,w_3,w_39  TRANSPORTED: phi := IsIsomorphic(ours,theirs) (0.06 s), then read
                         DefiningPolynomials off directly. w_2 = (x+z,-16y,x-z) was NOT in the
                         obvious candidate set. Canonical because #Aut = 8 and ABELIAN, so the
                         Isom-torsor choice of phi does not change the answer -- checked, not assumed
    93_1   w_3,w_31    + THE FULL CURVE, built from their published PAIR in P(1,3,1,1), genus 5.
                         No manual_isomorphism: the helper's construct-the-CRV-isomorphism branch
                         is hundredths of a second vs the 10 h+ IsIsomorphic

⚠⚠ **`IsIsomorphism` reported FALSE for `93_1`'s two involutions, and they are CORRECT.** Toric /
weighted-projective breakage (Magma #123 territory). Taken at face value it would have rejected two
correct transcriptions from the paper. Substituting into the defining polynomials leaves both
literally unchanged, with a control `x -> x+z` that breaks. ⇒ **On `P(1,w,1,1)`, `IsIsomorphism` is
not the instrument to check with — and the fallback must be a level BELOW the instrument, not a
sibling at the same level, since every map-level predicate routes through the same machinery.**
⚠ **Status: matrices verified as automorphisms of the right curves; LABELLING is unverified.** Only
a pipeline run tests that, and that is the claim that matters. `87_1`/`111_1` running, `39_2` queued,
`93_1` not started.

⚠ **CORRECTION, same day.** This block first said "these four files are UNCOMMITTED on purpose".
They are not: `git add -A` in THIS commit (`638223e`) swept all four in and they were pushed. The
intent was real and was stated twice, but the command did not implement it. **They are offline
tests, so CI does not run them and nothing went green on unverified matrices** -- but the repo
carried four labelling claims whose verification was still in flight.
⇒ Two things that follow. **`git add -A` does not respect an intention held only in your head**;
stage explicitly when part of a tree is deliberately held back. And **a claim about repo state is
as checkable as any other claim** -- one `git show --stat` would have caught it at the time, and it
was caught only when a later `git status` looked surprisingly clean.

**Verification status of the four, which is the part that matters:**

    87_1   ✅ GREEN 4660 s, AND negative-controlled: swapping the matrices goes RED in 3008 s
           naming "all 2 labelled involution(s), 5 candidate map(s) tried" -- so the {1} key WAS
           matched and the pass is not a lucky phi
    111_1  🔄 running, 3 h+ at 98% CPU
    39_2   ⏳ queued
    93_1   ⏳ not started

⚠ **A bare `Success!` from these files is NOT evidence on its own.** `run_tests.m` does not enable
`ShimuraQuotients` verbosity and the helper prints its comparison counts only under
`vprintf ... 1`, so a green run is equally consistent with "the involutions matched" and "the `{1}`
key was never matched, so `ws_data` was skipped in silence" (`if not ws_def then continue`). The
swapped-matrix control is what separates those, and it must be run per file.

### ⚠ `X0_*` re-derivation is 43%, and BOTH of my "corrections" to the recorded 41% were wrong

    152 cover_data keys / 355 populated model keys across the 40 tested bases = 43%

The recorded 41% (2026-09-09) was sound. I first got 17% by using the denominator across ALL 113
model files rather than the 40 tested bases — the note says "across the 34 tests". Then 44%, still
short, because **three model files write keys as `models[[ 1 ]]` rather than
`models[[Integers()|1]]`** (`21_2`, `58_5`, `34_3`) and my regex silently missed them. Eight bases
still check 1 of 15. ⇒ two wrong-denominator errors in one afternoon; the number was only ever
flagged as uncertain, never asserted.

### ⚠ A PEER SESSION `pkill`ed EVERY Magma PROCESS — and the misdiagnosis cost more than the kill

`ModFrm-CrvMod` ran `pkill -f "magma.exe"`, which kills every Magma the user owns. It took out
`X0_87_1`, `X0_111_1` and a `39_2` probe. **Signature: exit code 144 (NOT 143/SIGTERM) and a log
holding only the banner, ~20 bytes.** I diagnosed it as "three concurrent launches" and told the
user — a wrong cause that would have gone into the notes permanently, and was corrected only
because that session volunteered the disclosure. ⇒ **before blaming your own launch for a
simultaneous multi-job death, ask: `ListAgents` lists the peers.**
⚠ `pgrep -f magma.exe` then killing PIDs is NOT the remedy — that is a pattern kill with extra
steps, and there is no safe pattern here because `build/debug/native/magma.exe` is shared by every
worktree. Use `TaskStop` by task id, or an explicit PID with its cwd confirmed in the same command
(`lsof -a -p <pid> -d cwd -Fn`).

### lava: all three oracle-bearing bases running, and lava is NOT the safe harbour PLAN assumed

`95_1` and `119_1` from 11:40, `159_1` added 16:47, `$HOME/lavarun`, 0 `EXIT` lines.
⚠⚠ **lava runs the SAME `earlyoom --prefer (...|magma)`.** The move buys headroom, not immunity.
But the 13–82 GB per-job figure that justified holding the third back is from **lovelace and did
not transfer**: measured here, `95_1` 2.2 GB and `119_1` 8.7 GB against 125 GB, and `119_1`'s spike
to 17.4 GB came back down — a build phase, not a growth rate. ⇒ do not re-derive a resource ceiling
from another machine's numbers.
`run_base.sh` records `128+N` for a signalled child AND verifies the model file, because a bad
Magma command line **exits 0** (`OUTDIR=` instead of `OUTDIR:=` is read as a filename). Check with
`ssh -J lovelace lava 'grep "^EXIT" ~/lavarun/out/DRIVER.log'`.

## Handoff — 2026-09-22 — SIX MODELS, `358_1` RESOLVED, AND ~400 CPU-HOURS LOST TO `earlyoom`

### ✅ Six models collected and committed

    14_37    7 of 14 keys    9bb63ec      6_89    11 of 15    7796a99
    6_73     8 of 15         88634f9      6_113    8 of 15
    6_107    8 of 15                      6_137    8 of 8  (a COMPLETE set)

All verified: `ModelChecks` 12414/0 over **113 model files**, and each base individually
negative-controlled (corrupt one leading coefficient -> red). ⚠ None is a Guo-Yang base, so all six
are at `10_61` evidence level -- no external oracle.

### ✅ `358_1` RESOLVED -- the last base carrying "no verdict"

    ladder 5 3 2 2 2 2 2  at P = 144 190 266 350 450 550 700

Flat at **deficit 2 across five rungs** (P = 266 -> 700), the same invariance signature as the
`158_1` positive control. ⇒ OBSTRUCTED with a **2-dimensional obstruction space**, the EIGHTH such
base (`166_3 22_19 74_7 10_67 58_13 302_1 334_1` + `358_1`). **Known obstructed 133 -> 134**, still
a lower bound, and **no base is left in the unresolved state**.
⚠ `WeaklyHolomorphicBasis` alone took 29152 s (8 h), total wall 41240 s -- which is why this base
sat header-only for days and looked stalled. Buffered stdout + a slow first stage reads exactly
like a dead job; it was not.

### ⚠⚠ ~400 CPU-HOURS LOST: `earlyoom` KILLS LONG MAGMA JOBS ON lovelace

`bk3/DRIVER.log` records `EXIT 143` (SIGTERM) for **all five re-runs within 41 seconds**:

    123_1  2026-09-18T23:24:10     119_1  23:24:17     115_1  23:24:49
     95_1  2026-09-18T23:24:50     159_1  23:24:51
    314_1  2026-09-18T17:41:55     (62_7 died 17:48, not even recorded in DRIVER.log)

Each had run **66-87 hours** and was still in `Computing Borcherds forms`; `62_7` died at
`Computing equations of covers`, i.e. the LAST stage. No reboot (uptime 20 days), no stop script
touched since 09-15, and lovelace runs **`earlyoom -r 3600 --prefer ^(...|magma)$`** -- killing
preferred processes one after another under memory pressure is exactly this signature, on a shared
box where other users can squeeze 2 TB.

⇒ **Multi-day Magma jobs on lovelace are not survivable as currently launched**, and the loss is
silent: `/usr/bin/time` still reports `Exit status: 0`, so only the `Command terminated by signal 15`
line and `EXIT 143` in the driver log reveal it.
⇒ Use **lava** (`ssh -J lovelace lava`, 32 cores, near-idle) for runs of that length, or expect to
lose them. See [[remote-machines-lovelace-lava]].

⚠ **CORRECTED 2026-09-22, same day: this block first called the five "the Guo-Yang re-runs" and
said `95_1` "was the one with oracle value". Both are wrong.** The five are the jobs `earlyoom`
reaped -- a kill cohort. The Guo-Yang label came from `genmodels.m`'s
`vx_skip = {95_1,115_1,123_1,129_1}`, which groups by the vx defect and not by the paper. Against
the equation tables (43 `\multirow{1}{*}{\text}` cells, label two lines above each; 43 cells to 43
distinct labels, no repeats): **`95_1`, `119_1` AND `159_1` all have published equations**, while
**`115_1` and `123_1` are not in Guo-Yang at all** -- no `X^{115}_0` or `X^{123}_0` label anywhere
in the source. So the expensive loss is three oracle-bearing bases, not one.
⚠ `119_1`/`159_1` are not transcribed in `tests/GuoYangEquations.m` (it has ten + `93_1`); both are
degree-20 and wrap across `\\`, so transcribe by hand from the JOURNAL version --
[[guoyang-journal-version-differs]]. A textbook instance of the repo's own rule: the arithmetic
about *which machine killed what* was right, and the object -- *which bases carry an oracle* -- was
never checked.

### The even-correction / quadratic-CM line: where it stopped

Full account: `vvdata/weyl-campaign/even-correction/HATCH-EXISTS.md` §5-§7 (campaign), probes in
`probe-witness-guard.patch`. In brief:

* **Condition 4 is the binding filter.** Conditions 1-3 are cheap lattice tests with cheap
  candidates (floors 2-10 across controls), but the cheapest 1-3 candidate fails condition 4 at
  both bases tried (`34_3` cost 8, `14_3` cost 2). ⇒ the 1-3 cost floor does NOT predict legality.
* **`QUADCONSTRAINTS.md` §9's "minimum cost 24 at `34_3`" is VINDICATED** and now explained:
  condition 4 wants amounts `= 0 mod 6`, condition 3 kills the cheap ones that are. My same-day
  claim that cost 8 refuted it was wrong.
* ⇒ **The quadratic-CM route cannot be validated at `34_3`**: legal cost 24 and a rational supply
  hard-capped at 10 give `dim P(B) ~ 19`, where the Groebner step ran 6 h for nothing.
* ⚠ **A cheap condition-4 screen does NOT work by calling `SchoferFormula` directly** -- its
  baseline control (the UNPERTURBED target, at a base that builds) also came back FAIL, which is
  what caught it. `RationalNumber` is applied to the table entries AFTER `row_scales`, so the raw
  value is the wrong object. Any future screen must test the unperturbed target first and refuse
  to report if that fails.

### ✅ Also landed: the coprime filter fix (`0ca6e37`, main)

`BorcherdsForms.m:863` still passed the `coprime_to_level` default of `true`, though the same
filter was flipped OFF by default at `SchoferFormula.m:1153` on 2026-09-07. These points are
divisor support only and never Schofer-evaluated -- the `#pts<3` fallback says so itself -- so the
filter's purpose does not apply. Pool grows at **all 16** `N>1` bases (15_2 2->10, 21_2 2->9,
38_5 4->9); **14 of 16 re-derivation tests byte-identical both ways**, including `26_3` (pinned
base_label) and `22_5` (drift flags). `PTSCOPRIME=1` restores the old behaviour at THIS site only
-- deliberately not `CMCOPRIME`, which `SchoferFormula` also reads, so a control using it is
confounded (it re-enables the Schofer-side filter and alone breaks `14_3`/`39_2` -- which it did,
and nearly read as evidence FOR the change).
⚠ NOT a rescue for CM-starved bases: `134_3` has 0 coprime and 2 total points, still under the bar.

## Handoff — 2026-09-16 (evening) — THE EVEN-CORRECTION HATCH EXISTS EVERYWHERE SURVEYED, AND IT HAS AN UNGUARDED HAZARD

Full account, with every measurement and control:
`vvdata/weyl-campaign/even-correction/HATCH-EXISTS.md` (campaign, `7b60f15`/`69d77cb`/`64f7e73`).
Probe code: `vvdata/weyl-campaign/even-correction/probe-mod2.patch`.

### ✅ Three more models collected

    14_37   7 of 14 keys    ModelChecks 12036/0   neg ctl 62/5    9bb63ec
    6_73    8 of 15 keys    ModelChecks 12162/0   neg ctl 62/5    88634f9   (recovered: the
    6_107   8 of 15 keys                          neg ctl 62/6              SIGTERM'd 09-14 run)

None is a Guo-Yang base, so all three are at `10_61` evidence level -- no external oracle.
⚠ The five Guo-Yang re-runs (`95_1 119_1 159_1` at ~27 h, `115_1 123_1` at ~7 h) are ALL still in
`BorcherdsForms`. Their logs are 28 bytes and Magma buffers, so **CPU-vs-elapsed is the only
progress signal** -- all five read 99.6%, i.e. healthy, not stuck.

### ✅ The hatch's EXISTENCE question is settled, and it is not what the record assumed

`PROBE_MOD2` decides "does ANY even correction exist" by lattice membership -- `target in L + 2Z^nds`
-- in ~20 s per base. `PROBE_INTSWEEP` could never answer it: it tests ONE discriminant at a time
against a guessed amount list, and at an obstructed base the charge equation `a*phi_j = -phi(target)`
PINS the amount per discriminant, so the corrections that exist need TWO OR MORE discriminants and
were never searched for. Result: **TRUE at all 28 bases with annihilator data**, restricted to
genuine discriminant coordinates. ⇒ conditions 1+2+3 are jointly satisfiable everywhere surveyed;
**condition 3 is not a filter on existence.** The test is strictly stronger than the known parity
criterion at **16 of 28** bases (scaling-free unit-perturbation control), so this is real
information, not a restatement of the 28/28 parity survey.

⚠ **The `QUADCONSTRAINTS.md` cost verdict does not transfer to obstructed bases.** Its 24 (at
`34_3`) and 48 (at `35_1`) are measured at UNOBSTRUCTED controls, where `phi = 0` makes every amount
legal and integrality picks the expensive discriminants. At an obstructed base the charge equation
sets the price instead.

### ⚠⚠ THE HAZARD: the target coordinate is NOT the divisor coefficient

    generic disc             factor 1    target +2 -> divisor +2
    AL fixed-point disc      factor 2    target +2 -> divisor +1   (|d| = m or 4m, m | D*N)
    d = -3, -4               factor 4+   target +4 -> divisor +1   (compounds with the above)

An ODD divisor change breaks condition 1 -- the cover is not preserved and the "model" is a
DIFFERENT CURVE, undetectable at an obstructed base where there is no oracle. **The existing probe
printed `DIVISOR MISMATCH` and CONTINUED**, so this was latent in the machinery. ⚠ The recorded
`34_3` results are unaffected (disc -164 is generic; re-confirmed this session), but `PROBE_EVEN`'s
"prefer the LARGEST |disc|" heuristic does nothing to avoid the dangerous ones.
⇒ The check now computes the ACTUAL change (`div_f - ram`) and ERRORS unless every entry is an even
integer. **Keep the guard even if the rule looks complete** -- the first rule covered the factor-2
case and still let disc 4 through; the guard, not the rule, caught it.

### Corrections to claims made earlier the same day

* **"cost 10 at `38_5`" -- WRONG, it is 18.** Every cost-10 witness used disc 19 or 20 at amount 2,
  i.e. the factor-2 discriminants. Demand is `2g+23`, not `2g+15`.
* **"22 of 28 bases have 2-torsion" -- WRONG**, a scaling artifact (elementary divisors of `dM*A`).
* **"the lattice test is just the parity criterion" -- WRONG**, generalised from `38_5`, which is one
  of the 12 bases where they happen to coincide.
* **A two-amount elimination** `(y2_a)^2/y2_b` removes the degree cost EXACTLY at `34_3` (7/7 keys,
  23 independent checks, negative-controlled) and is **PROVABLY IMPOSSIBLE** at an obstructed base:
  legal perturbations all carry the same nonzero charge, so no two are proportional and the residue
  never cancels. Do not re-attempt it.

### In flight at the end of this session

* `38_5` auto-hatch pipeline run (per-key corrections, cost 18, all divisor changes verified even).
  First run to get past `BorcherdsForms` at an obstructed base -- it is what will finally measure
  `#rat` and condition 4 there, neither of which has ever been observed at an obstructed base.
* `34_3` **quadratic-CM-point validation**: fit at the TRUE perturbed degree (`DEGBUMP=24`) so the
  rational points UNDER-determine `f` instead of contradicting it, then let the existing quadratic
  machinery cut `P(B)` down. The rational supply is capped (10 at `34_3` however many are asked for)
  while the extra points all arrive QUADRATIC, so quadratic points are the only supply that grows.
  ⚠ `require not IsEmpty(B)` in `QuadraticConstraintsOnEquations` means quadratic constraints are
  reachable only AFTER the rational fit succeeds -- i.e. exactly when they are not needed. That gate
  is the hatch's real blocker, not the CM supply.
* `358_1` pole-ladder extension: alive, 99.6% CPU, **32 GB RSS**, ~7 h inside
  `WeaklyHolomorphicBasis`. ⚠ lovelace runs `earlyoom` with `--prefer ...|magma`.

## Handoff — 2026-09-16 (later still) — COLLECTION: `14_37` LANDS, `14_71` RESOLVES, NOTHING ELSE IS DONE YET

Picked the backlog-collection track back up. lovelace load 107/256 (shared, other users' jobs
still dominate — did not launch anything new).

### ✅ NEW MODEL: `14_37` — 7 of 14 keys, via Hauptmodul rebase

`~/shimura/bk3/14_37.log` finished after **25.5 h wall** (started under the 2026-09-15 batch).
7 of 14 cover keys came back empty from the direct construction and were filled by the
Hauptmodul-rebase sweep (`"sweeping 1 Hauptmodul root(s) on base 5699"` — the same mechanism as
`22_5`'s recovery). Collected, `ModelChecks` run against the WHOLE suite with it included
(**108 model files, 12036 checks, 0 failures**), and a negative control (corrupt one leading
coefficient in the `[1,7,74,518]` entry) correctly goes red (62/5). Committed, `9bb63ec`.
⚠ No external oracle at `14_37` — `10_61`/`34_11`/`74_5` evidence level, not a Guo-Yang base.

### ✅ `14_71` RESOLVED via `defhi` (P=700 extension): OBSTRUCTED, deficit 2

`~/shimura/defhi/14_71.log` extended the pole ladder past `deficit.m`'s `P<=266` cap. Reads
`3 3 2 2 2` at `P = 498, 550, 574→realigned, 550, 700` — **flat at deficit 2 across P=550→700**,
150 of pole-order headroom past the point it last moved. Same invariance signature as the
`158_1` positive control from the prior collection pass. Corrects the earlier `deficit 3` reading,
which was mid-descent, not the true value. Ordinary (1-dimensional), not a new 2-dim candidate.
Folded into `vvdata/weyl-campaign/obstructed-rerun-2026-09-10/screened-2026-09-14.txt` on the
campaign branch (`cdcf801`, not yet pushed). ⇒ **known obstructed 132 → 133**, still a lower bound.
`358_1` is the only remaining unresolved verdict — its own `defhi` extension is still running on
lovelace as of this collection (`~/shimura/defhi/358_1.log`, header only so far, ~4 h elapsed).

### Everything else on lovelace: still in flight, nothing new to collect

    bk3 pipeline    115_1 123_1 (relaunched, both at "Computing Borcherds forms")
                    62_7 (at "equations of covers"), 6_107 6_113 6_73 6_89 6_137 (all at
                    "candidate discriminants"/"CM points"), 314_1 (deferred one cover on an
                    ambiguous sign, W={1,157}, continuing -- not a failure), 95_1 119_1 159_1
                    (all still at "Computing Borcherds forms")
    defic5 screens  ~25 large-D N=1 / D=6 large-N targets still header-only (never returned;
                    see the campaign log for the list) -- check again before re-screening any of
                    that range

Two uncommitted scratch files sit at the repo root from a prior session, `tests/_probe_gy.m` /
`tests/_probe_gy2.m` (quick Guo-Yang model-completeness listings, not registered tests since they
lack the runner's naming convention). Left alone — harmless, not blocking anything, and per
`CLAUDE.md` scratch scripts belong on the campaign branch under `vvdata/weyl-campaign/`, not here;
worth moving or deleting next time this file is touched.

## Handoff — 2026-09-16 (later still, second) — THE KRY CHAPTER 7 LEAD IS CHASED, AND RETIRED

Fetched the KRY book (PDF), `pdftotext -layout`'d it, and read Chapter 7 — §7.1 (statement of
Theorem C, the height-pairing/Fourier-coefficient identity, and Theorem 7.1.1) and §7.6
(the explicit local formula `ν̃_p(T)` at ramified primes, Props 7.6.2–7.6.4) in full; the rest of
the chapter only by section-header structure. Re-read `vvdata/weyl-campaign/deficit.m` alongside it
to pin down exactly what the target computation is.

**Verdict: retired as a deficit predictor.** Theorem C's height pairing is an Arakelov-theoretic
real number on the *integral model*, tied by the book's own Ch. 9 to central derivatives of
`L`-functions — the Gross–Zagier analogue for Shimura curves. `deficit.m`'s number is a plain
finite-dimensional linear-algebra rank over q-expansion coefficients, with no scheme, height, or
archimedean data anywhere in it. The §7.6 local formula, despite having the right "local density at
a ramified prime" flavor (same family as this project's own `κ_p`/`SchoferFormula.m`), computes an
intersection multiplicity for ONE fixed pair `(t1,t2)`, never a rank over a basis against a target
SET — that question doesn't appear in Ch. 7 at all. If anything this REINFORCES the earlier
`S_{3/2}`-dimension finding rather than circumventing it: the closest global quantity KRY computes
is exactly the "hard", central-L-value-flavored kind the deficit rank was already diagnosed as.
Full argument and the exact citations: `vvdata/weyl-campaign/kry-ch7-notes.md` (campaign, `dba02a3`,
pushed). `PLAN.md`'s "⚠ THE KRY LEAD" section has the recorded verdict.

⇒ **The theory-arc parallel track is now exhausted of concrete leads.** Nothing currently on record
suggests a closed-form deficit predictor exists; the model backlog (collection track) remains the
only track producing results. Not pursued further this session — reverting to backlog collection.

## Handoff — 2026-09-16 (later) — THE `S_{3/2}` DIMENSION IS COMPUTABLE, AND IT ANSWERS THE WRONG QUESTION

`paper/DRAFT-borcherds-obstruction.md` §5 asked for `dim S_{3/2}(ρ_L^*)` in closed form as a
predictor for `deficit.m`'s measured number. Both halves are now settled.

**Route B (naive) refuted immediately**: `deficit = genus(X_0^D(N))` fails on the first two bases
tried — `146_1` (genus 7, deficit 0) and `194_1` (genus 9, deficit 0) don't even share a genus with
`38_5` (genus 9, deficit **1**), so genus alone cannot determine the deficit.

**Route A (the real Riemann–Roch formula) is now sourced, implemented, and validated.** The formula
is Borcherds' own, from the GKZ paper itself (Duke 97 (1999), p. 9) — not the Bruinier/Kuss citation
the draft originally guessed at. Implementing it hit a real wall: this repo's `WeilRepresentationST`
builds the honest `|L'/L| × |L'/L|` matrix, and `|L'/L|` runs 72,200 to 1,299,272 on the calibration
bases — far too large to diagonalize (the naive attempt was killed after 7+ minutes on the smallest
one). The fix: every quantity the formula needs reduces to `O(n)` Gauss sums via two algebraic
tricks (a projector isolating one eigenspace of the negation involution, and a shift-bijection that
factors `tr((ST)²)` into a product of two simpler sums) — **all six resulting trace identities were
checked bit-for-bit against the real matrices on `6_1`/`10_1` before being trusted on anything
larger.** Banked as `vvdata/weyl-campaign/gksz-dim-formula.m` (campaign, `6565a95`).

**The result is not what was hoped for, and that is itself the finding**: `dim M_{3/2}(ρ_L^*)` comes
out in the **thousands** (`38_5` → 1594, `146_1` → 888) against measured deficits of 1 and 0. This
is not a bug — Serre duality (§1) says the full obstruction space is dual to *arbitrary* principal
parts, but `deficit.m` only ever asks to hit a small, fixed target set of tracked CM-divisor classes.
The honest statement is `deficit = rank(S_{3/2}(ρ_L^*) → target*) ≤ dim(target)`, bounded by the
*target's* dimension (small, a dozen or so classes), not the ambient one — so the Riemann–Roch number
is essentially irrelevant to the actual deficit. Whether that rank degenerates depends on whether
*specific Fourier coefficients* vanish at the *specific* discriminants in the target set — a
Waldspurger-type coefficient/central-L-value question, which is "hard" arithmetic, not a "soft"
dimension count, and structurally cannot be answered by any Riemann–Roch or trace-formula argument.
⇒ **No closed-form deficit predictor is expected to exist along this route.** What survives is a
real, validated, fast upper bound on the deficit — useful, but not the predictor the draft wanted.

Full argument, with all the algebra: `paper/DRAFT-borcherds-obstruction.md` §5.

⚠ **Paused here, deliberately, for a fresh session to pick up** — the live lead (Kudla–Rapoport–Yang,
*Modular Forms and Special Cycles on Shimura Curves*, Ch. 7's inner product formula, motivated by
`T` reducing to the supersingular locus at ramified primes) is real but was not chased into an
implementation; the conversation that found it was already long. Full plan, with the exact chapter,
what's already confirmed, and the calibration bases to validate against: `PLAN.md`, "PARALLEL TRACK
— added 2026-09-16".

## Handoff — 2026-09-16 — POST-RESET COLLECTION: 57 NEW OBSTRUCTED BASES, TWO NEW 2-DIM SPACES

The local Mac was reset mid-session. **Nothing was lost**: nothing was running locally (no Magma
process), and every job in flight was already on `lovelace`, detached (`ppid 1`) and unaffected.
Full detail of everything below: `vvdata/weyl-campaign/obstructed-rerun-2026-09-10/screened-2026-09-14.txt`
(campaign, `c834599`).

### ✅ Branch hygiene: campaign was one merge behind `main`, now fixed

`git diff origin/main origin/m0-theta-campaign --name-only -- ':!vvdata/weyl-campaign/*'` showed a
real (non-doc) divergence on `tests/InternalBorcherds.m` — campaign's last merge point (`928bf22`)
predated `main`'s `714e831`/`5b57a2d`. Merged and pushed (`d55b5c0`); the invariant is clean again.

### ✅ The screening collection: known obstructed 75 → 132, 2-dim spaces 5 → 7

Two batches finished on lovelace since the last collection pass:

* **`defic4`'s remaining 7 bases** (its first 3 were already recorded): `278_1`, `326_1`, `346_1`
  obstructed; **`314_1` CLEAR** (pipeline launched, deferred one cover on an ambiguous CM sign —
  not a failure); `302_1` and `334_1` obstructed at **deficit 2**; `358_1` still descending at
  `P=266`, no verdict.
* **`defhi.sh`** (new: extends the pole ladder from `deficit.m`'s `P<=266` cap out to `P=700`) was
  built to settle the three still-moving cases, using `158_1` — a known obstructed base — as a
  positive control: it stays flat at deficit 1 for 630 of pole-order headroom past the old cap,
  which is the evidence that a flat `P<=266` tail is the real invariant and not truncation.
  `302_1`/`334_1` hold at deficit 2 across the same range. ⇒ **2-dimensional obstruction spaces:
  FIVE → SEVEN** — add `302_1`, `334_1` to `166_3 22_19 74_7 10_67 58_13`; `326_1` resolves to
  ordinary (dim 1), not the fourth candidate it looked like at `P<=266`.
* **`defic5`**, a systematic odd-prime-`N` sweep (`D = 6,10,14,22,26,34,38,46,58,62,74,82,94,118,
  122,134,142,146,194,206,326,362,386,394`) plus a batch of large-`D` `N=1` targets: **52 more
  obstructed verdicts**, checked against `bases49.txt` and every prior screened/stumbled list —
  zero overlap, script-checked not assumed. One more unresolved case, `14_71`.

⇒ **Known obstructed: 49 recorded + 5 stumbled into + 78 screened = 132**, still a lower bound.
⚠ **Read the jump as "we screened a systematic sweep", not as a new density statement** — these
bases were chosen, not sampled. ⚠ **Unresolved, not obstructed and not clear**: `358_1`, `14_71`
(both still descending at their last rung, neither got the `P=700` extension). 25 more targets in
the `defic5` batch never finished (still running as of collection, header-only logs).

### `bk3` pipeline: nothing new to commit yet, all still in flight

`62_7 6_107 6_113 6_73 6_89 6_137 314_1` are all still computing (mostly at "equations of covers").
**`95_1` is now actually running** — `cc8beb8` (campaign, already landed before the reset) dropped
the stale `vx_skip` guard in `vvdata/weyl-campaign/genmodels.m` once its justifying defect (fixed
2026-09-05, `d9b52d0`) was ten days gone; `95_1` is a Guo-Yang target with a published equation and
was being `quit`-skipped for no remaining reason. `115_1`/`123_1` — the other two bases the guard
was fencing off besides `95_1`/`129_1` — have **not** been relaunched since.

### Next

* Collect `bk3` once any of `62_7 6_107 6_113 6_73 6_89 6_137 314_1 95_1` finishes; verify with
  `VerifyModelSet` + a negative control before committing, per the usual discipline.
* Relaunch `115_1`/`123_1` now that `vx_skip` is gone.
* Extend `358_1` and `14_71` with `defhi.sh` (`P=700`) to get real verdicts.
* The ~25 not-yet-finished `defic5` bases will keep returning verdicts; collect and fold into the
  same log file rather than starting a new one.

## Handoff — 2026-09-15 (night, 7th) — RETRACTION: "20 OF 23 KEYS HOLD DIFFERENT CURVES" IS WRONG

The section below this one (night, 5th) claims 20 of 23 multi-entry genus-1 keys hold different
curves, and the one below that (night, 4th) calls two `GonzalezRotger.m` entries "the wrong torsor".
**Both claims are retracted.** The arithmetic was right; the criterion was wrong.

### WHAT WAS WRONG

`Genus1Classes.m` asserted that entries under one key must be **GL2-equivalent** as binary quartics.
That is not what "the same curve" means. A degree-2 map from a genus-1 curve `C` to `P^1` is a
`Q`-rational degree-2 divisor class, and those form a torsor under `E(Q) = Pic^0(C)`. So when `E(Q)`
is nontrivial, **one curve carries several inequivalent quartic models** — inequivalent as binary
quartics, identical as curves.

`E(Q)` is nontrivial for every Jacobian in play. Measured:

    30a6 [2,2]   78a2 [2,2]   30a2 [2,6]   30a3 [2]   102b3 [2,2]   130b2 [2,2]
    154a1 rank 1 [2]     138a1 rank 1 [2]     426b1 rank 1 [2]

Not one trivial group. Two bases picking different degree-2 classes produce exactly the observed
signature — inequivalent quartics, same Jacobian, 14 of 14. ⇒ **"14 of 14 share a Jacobian" is not a
smoking gun, it is the EXPECTED outcome**, and the "one mechanism, not twenty mistakes" paragraph in
the night-5th section argues for the wrong conclusion from it.

Corroborating, and the check that should have been run first: the entries of all eight keys probed
have **identical everywhere-local solubility profiles** (real place and every prime up to 47).

⚠ This is the failure `CLAUDE.md` opens with — correct arithmetic about the wrong object — and it
happened while building a test to catch that very class. The tell was available and ignored: a
check that calls 87% of its population defective is far more likely to have the wrong criterion than
to have found an 87% defect rate.

### WHAT REPLACES IT

`tests/Genus1Classes.m` now compares the **Jacobian**, which is an invariant of the curve and blind
to the choice of degree-2 class — the honest genus-1 analogue of `ConicClasses.m`'s Brauer class:

    23 multi-entry genus-1 key(s), 23 pairwise Jacobian comparison(s)
    all 23 multi-entry key(s) agree on the Jacobian

⚠ Necessary, NOT sufficient: inequivalent torsors of one `E` share a Jacobian. At genus 0 the Brauer
class happens to be a COMPLETE invariant of a conic; at genus 1 the Jacobian is not, and that gap is
real and unclosed. The file's header now carries the whole argument, so the GL2 criterion is not
re-invented.

Controls RUN: quadratic-twist one entry (changes its Jacobian) -> RED, naming `6_17 W=[1,2]` and
printing both `aInvariants`; read a single model file -> RED on the count guard; restored -> GREEN.

### `GonzalezRotger.m`, corrected

`KNOWN_TORSOR_DRIFT` is renamed **`NO_EXHIBITED_ISO`** and no longer claims the entries are wrong.
The asymmetry note now gives BOTH reasons a failure proves nothing — `IsGL2Equivalent` need not
return the full orbit, and our model need not be GL2-equivalent to theirs at all. The error message
says "could not be proved isomorphic", and is flagged in-file as a "something changed, look at it"
guard rather than a defect claim.

⚠ **Nothing should be deleted from `models_6_5.m` or `models_6_13.m`.** The night-4th section says
those entries "should be removed from the data"; that recommendation is withdrawn.

### A second, smaller retraction in the same class: `test_bp_KY`

The `InternalBorcherds` revival said its 118 mismatches meant "the discrepancy is confined to the
`Wpoly2` branch". **Also wrong, also the wrong object.** Arbitrated with the brute-force density
oracle at RANK 1 — which is what `test_bp_KY` builds, where production and `Whittaker2.m` both work
at rank 2:

    rank 1, Q = [2 kappa], kappa in {1,2,3,5,6,7}, m in {1,2,3,4}
       -> 24 comparisons, library vs oracle: 0 disagreements   (k = 14, re-confirmed at k = 18)

⇒ **`Wpoly2` is correct at rank 1 as well as rank 2.** What fails is the test file's own right-hand
side — its transcription of [KY, Prop 5.1] at `p = 2`. The file's guess at "a sqrtp factor I am
missing" is in the right half of the identity but the wrong place; a missing `sqrtp` would have moved
the odd primes too. Corrected in-file.

⇒ **Tally for the day: three blind spots opened, and the library was clean in all three.** Every
defect found today was in a test, a harness, or one of my own readings.

### What survives, unchanged

The positive half: where an isomorphism IS exhibited it is a proof, certified as an exact identity
in `Q[x]`, and 11 of 13 genus-one entries now carry one where previously all 11 rested on a Cremona
label. `NPROOF` is asserted separately so a run that degrades back to matching invariants goes red.
The oracle stands at 49 comparisons. `ConicClasses.m` is untouched and still passes.

## Handoff — 2026-09-15 (night, 6th) — THE p=2 BLIND SPOT IS CLOSED, AND THERE WERE FIVE SHAPES

`tests/Whittaker2.m` validated `Wpoly2` on exactly two 2-adic Jordan shapes. Production feeds it
**five**. All five now have expected values pinned, and **the library is correct on every one**.

    H0      unimodular hyperbolic        [[0,1],[1,0]]        was covered
    H1      2-modular hyperbolic         [[0,2],[2,0]]        was covered
    d1d1    odd type, v_2(det) = 2                            NEW
    d1d2    odd type, v_2(det) = 3                            NEW
    even1A  2*[[2,1],[1,2]]                                   NEW -- and I had MISSED it

### ⚠ MY OWN SHAPE CENSUS WAS WRONG, AND WRONG IN THE REPO'S SIGNATURE WAY

I reported four shapes. My classifier recorded only the 2-VALUATION of the off-diagonal entry after
`pAdicDiagonalization`, so it labelled every 2x2 block "H" — conflating the two inequivalent even
binary `Z_2` lattices. There are exactly two up to scaling: `H = [[0,1],[1,0]]` (det -1) and
`A = [[2,1],[1,2]]` (det 3). Checked independently: `det A / det H = -3`, and `-3 = 5 mod 8` is not
a square in `Z_2^*`, so they are NOT interchangeable — and `Wpoly2` routes them down different
branches of Yang's formula. Correct arithmetic, wrong object, again.

`even1A` is what the ODD discriminants give: `6_1` at -3/-19, `14_1` at -11, `10_1` at -3, `34_1` at
-3/-11, `38_1` at -11/-19, `26_1` at -11/-19, `6_5` at -19, `10_3` at -3. All twelve real Grams
return the same value row and the entry `3` appears nowhere else in the file, so it is a distinct
branch and not a relabelling.

### THE VERDICT: THE LIBRARY IS CORRECT

Brute-force representation-density counting, independent of Yang's formula and of the library:

    calibration on H0/H1        32/32 exact, ratio 1        (reproduce a KNOWN value first)
    mu = 0, all five shapes     72/72 agree
    nonzero cosets             108/108 agree
    Magma re-derivation        786 comparisons, 0 mismatches (k=12, 1602 s)
                               966 comparisons, 0 mismatches (k=10, 122 s) with even1A

### ⚠ THE PLATEAU IS REAL AND LONG — NO REPEATS RULE, ANYWHERE

Measured across 300 `mu = 0` comparisons, the last `k` at which any approximant moved was **k = 8**;
across 666 coset comparisons, **k = 2**. `6_1`/-24 at `m = 32` reads

    1 1 1 1 1 1 1 2 2 2 2 2 2 2      (k = 3 .. 16)

— seven identical values before the true one. My own "two repeats" stopping rule produced NINE false
mismatches earlier this session. The committed table was taken at a fixed **k = 14, re-confirmed at
k = 16**, 228/228 agreeing at both.

### ⚠ TWO HARNESS FACTS THAT COST TIME

* **A command-line `kmax:=N` DOES NOT REACH A TEST FILE.** `run_tests.m` does `Read()` + `eval`, and
  Magma's eval scope cannot see top-level command-line assignments — the variable is simply
  unassigned and the default silently wins, which looks exactly like the flag being honoured, only
  slower. It cost a 27-minute run at the wrong `k`. `tests/_offline/Whittaker2Oracle.m` therefore
  reads `W2O_KMAX` from the ENVIRONMENT (`GetEnv`, as `Y2TWIST`/`M0PROGRESS` already do).
* **`ElementOfNorm` is ORDER-OF-CALL dependent, not merely seed dependent.** Inserting a
  `pAdicDiagonalization` call between iterations changed which `lambda` came back at `10_3`/-3 and
  `6_5`/-4 — to a 2-adically equivalent lattice, but a different Gram. ⇒ **Anything that pins a
  `lambda^perp` by re-deriving it is not reproducible.** The 25 Gram matrices are committed as
  LITERALS, and `two_adic_shape` re-checks that each row still is the shape it claims, so a typo
  fails before any value is compared.

### Counts and controls

**1067 comparisons** asserted (90 + 8 + 3 from the old parts, 300 new at `mu = 0`, 666 on cosets),
plus a shape census asserting 10 `d1d1` / 9 `d1d2` / 6 `even1A`. Runtime `0.15 s -> 1.1 s`; no brute
force in the committed file. Controls RUN:

    perturb one expected value            -> RED  (the value assert)
    empty the coset loop silently         -> RED  (ONLY the counter catches this)
    corrupt a Gram so its shape changes   -> RED  (the which-object guard, before any value)
    restore                               -> GREEN, "Done!  1067 comparisons."

⚠ One weak layer, flagged in-file: for `even1A` the coset values are all 1 (Yang's `K_mu` vanishes),
so that part of the coset sweep asserts little.

## Handoff — 2026-09-15 (night, 5th) — 20 OF 23 MULTI-ENTRY GENUS-1 KEYS HOLD DIFFERENT CURVES

The `6_5`/`6_13` torsor drift is not an anomaly. It is a CLASS, and `tests/Genus1Classes.m` now
records it.

### THE MEASUREMENT

Entries under one `(D,N,W)` key are the same quotient `X/W` computed over DIFFERENT BASES
(`all_eqns[k][base]`, `EquationsCovers.m`), so they must be isomorphic over `Q`. Across every model
file:

    890 populated keys ; 23 multi-entry keys of genus >= 1 (all genus 1) ; 23 pairwise comparisons
    ->  3 proved isomorphic,  20 NOT EVEN GL2-EQUIVALENT

Violating bases: `10_13 10_3 10_7 14_3 14_5 15_2 22_7 26_5 6_13 6_17 6_23 6_5 6_71`.
The three that pass: `10_7 W=[1,7]`, `21_2 W=[1,2]`, `6_13 W=[1,6]` — so `6_13` has both, and this
is per-KEY, not per-base.

### ⚠ ONE MECHANISM, NOT TWENTY MISTAKES — THEY ALL SHARE A JACOBIAN

Cross-checked against an invariant that does not use `IsGL2Equivalent` at all:

    same-Jacobian pairs: 14        different-Jacobian pairs: 0

Every multi-entry genus-1 key sampled holds inequivalent quartics whose Jacobians carry the SAME
Cremona label. Inequivalent quartics with one Jacobian are **different torsors** of it. Independent
transcription errors do not land on one `H^1` class fourteen times out of fourteen.

⚠ The suspect is the `y^2` sign/scale resolution (`find_y2_signs`, the per-disc signs): a choice
made differently per base would move the torsor while leaving the Jacobian alone. **UNTESTED.** Four
single-cause stories have been refuted by controls in this repo already; treat it as a lead.

### ⚠ THE GENUS-0 ANALOGUE PASSES — WHICH LOCALISES IT

    ConicClasses.m: 236 genus-0 conic(s); 41 multi-entry key(s); all internally consistent

Same data model, same generation path, same multi-base structure — and across bases the conics agree
every time. So this is **genus-1-specific**, not a general cross-base inconsistency.

### WHY NOTHING WAS DELETED

Two of the twenty are arbitrated by Gonzalez-Rotger (`6_5 W=[1]` entry 1, `6_13 W=[1]` entry 2) and
are pinned in `GonzalezRotger.m`'s `KNOWN_TORSOR_DRIFT`. **For the other eighteen there is no
oracle**, so removing an entry would be guessing which base got it right. Recording them makes the
class visible in CI and stops it growing; if the mechanism is found and fixed, the entries regenerate
correctly and the data question answers itself.

### THE TEST

`tests/Genus1Classes.m`, 0.14 s. Two entries are the same curve over `Q` iff there is a GL2(Q) map
`[a,b,c,d]` and a RATIONAL `lambda` with `f_2(x) = lambda^2 (cx+d)^(2g+2) f_1((ax+b)/(cx+d))`, and
the identity is certified in `Q[x]`.

⚠ The `lambda^2` is the mechanism, not decoration: `IsGL2Equivalent` decides equivalence MODULO ANY
SCALAR, and `y^2 = f` curves are isomorphic only when that scalar is a SQUARE.
⚠ Asymmetry, deliberate: a square constant PROVES isomorphism; finding none proves nothing, since
`IsGL2Equivalent` does not promise the full orbit. Hence a recorded list, not an assertion that these
curves differ.
⚠ It errors when a RECORDED violation starts passing. That firing is GOOD NEWS — it is what fixing
the mechanism looks like — and it forces the record to stay accurate instead of going stale.

Negative controls RUN: drop one violation from the record -> RED ("NEW key(s)"); record a violation
for a key that passes -> RED ("now PASS"); read one model file -> RED (count guard, "only 0
multi-entry key(s)"); restored -> GREEN.

## Handoff — 2026-09-15 (night, 4th) — EXHIBIT THE MAP: TWO COMMITTED ENTRIES ARE THE WRONG TORSOR

`tests/GonzalezRotger.m` compared the genus-one full curves by an INVARIANT — the Jacobian's Cremona
label. That is necessary and not sufficient: quartics with the same Jacobian can be **inequivalent
torsors** of it. Upgrading the check to exhibit the isomorphism found that this is not hypothetical.

### THE DEFECT

At **`6_5`** and **`6_13`** the `W=[1]` key holds TWO entries which are **not GL2-equivalent to each
other** — genuinely different curves — and **both carry the Jacobian label the paper states**:

    6_5   entry 1  Jacobian 30a6  NO Q-isomorphism to GR's curve     <-- spurious
    6_5   entry 2  Jacobian 30a6  isomorphic, T = [1,8,-1,0], lambda = 1/64
    6_13  entry 1  Jacobian 78a2  isomorphic, T = [0,1,-1/8,-3], lambda = 176
    6_13  entry 2  Jacobian 78a2  NO Q-isomorphism to GR's curve     <-- spurious

⚠ And the old check read only `models[key][1]`, so **at `6_5` it was certifying the entry that is
NOT the published curve** — and reporting a match. Same shape as the `10_3` `[1,2]` drift: internally
consistent, externally wrong, invisible to an invariant.

⚠ **THE DATA IS NOT YET FIXED.** `data/models/models_6_5.m` and `models_6_13.m` still carry the
spurious entry, recorded in the test as `KNOWN_TORSOR_DRIFT` so a NEW one goes red. Removing them
touches entry counts that `ModelChecks`/`VerifyModelSet` read, so it is a separate change.

### THE CERTIFICATE

Gonzalez-Rotger's own relation (Section 2, p.3), checked as an exact identity in `Q[x]`:

    f_GR(x) = lambda^2 * (c x + d)^4 * f_ours((a x + b)/(c x + d)),   lambda in Q

Given it, `(X,Y) |-> ((aX+b)/(cX+d), Y/(lambda (cX+d)^2))` is an isomorphism over `Q`. So nothing
calls `IsIsomorphic` or `Jacobian()` — these curves have no rational point by construction, so that
route returns `ERR` on both sides and prints a vacuous MATCH, and `IsIsomorphic` on a genus-0
`CrvHyp` is wrong on 2.29-10 (Magma#125).

⚠ **The `lambda^2` is the whole point.** `IsGL2Equivalent` decides equivalence of binary quartics
**modulo any scalar**; `y^2 = f` curves are isomorphic only when that scalar is a **SQUARE**. A
non-square constant is a different torsor. All 11 bases are GL2-equivalent to GR's quartic; the
square-class test is what separates them.

⚠ **Asymmetry, deliberate.** A transformation with a square constant PROVES isomorphism. Finding
none does NOT prove non-isomorphism — `IsGL2Equivalent` does not promise the full orbit. So proofs
are asserted; failures to prove are reported, never asserted upon.

### THE ORACLE NOW

    ok (13 genus-one entr(ies) checked, 11 of them by an EXHIBITED isomorphism;
        + 15 AL-quotient(s) + 21 splitness check(s); 0 base(s) without a usable W=[1] model)

**47 -> 49 comparisons**, and 11 of them are now proofs rather than invariant matches. `NPROOF` is
asserted separately from `NCMP`, so a run that silently degraded back to matching invariants goes
red. Negative controls RUN:

    empty KNOWN_TORSOR_DRIFT   -> RED  ("NEW wrong-torsor entr(ies)")  -- it SEES 6_5/6_13
    perturb GR's 14_1 quartic  -> RED  (the paper's own self-check fires first)
    cripple exhibit_iso        -> RED  ("no entry could be proved isomorphic...")
    restored                   -> GREEN

### ⚠ HOW MUCH OF THESE BASES ANY ORACLE ACTUALLY TOUCHES: 30%

Measured over the 11 Gonzalez-Rotger bases:

    98 populated model keys  ->  29 touched (30%)
    149 entries              ->  49 touched (33%)

Two-thirds of the model data on the bases where we HAVE an oracle is checked by internal consistency
alone. Concentrated in the large ones: `6_13` is 15 keys / 27 entries with 3 keys touched, `10_7` 15
keys with 2. ⚠ And **10 of the 11 have no `X0_*.m` test at all** — only `15_1` does.

⚠ **Gonzalez-Rotger cannot close that gap**, so do not write `X0_D_N.m` files for it: an `X0_*` test
is driven by hand-transcribed PUBLISHED COVER EQUATIONS (`cover_data`), and GR publish only the
genus-one full curve and the involutions. Such files could carry `cover_data[{1}]` and nothing else,
restating the key this section already proves. The 70% has no published source.

## Handoff — 2026-09-15 (night, later still) — A SUITE FILE THAT RAN NOTHING, AND WHAT IT HID

`tests/InternalBorcherds.m` reported `Success! 0.000 s` in every suite run on record. It defines
`test_kronecker_sigma`, `test_bp_KY` and `test_W` and **called none of them**, and nothing outside
called them either — the four apparent references are `test_Whittaker2`/`test_WeilRepresentation`,
prefix collisions, checked. A file of definitions asserts nothing. Compare `tests/Whittaker2.m` and
`tests/WeilRepresentation.m`, which invoke their procedure on the last line.

### What the silence hid: all three were broken, none of it mathematical

    the import named tests/BorcherdsProducts.m, but Wpoly/Wpoly2/Wpoly_scaled live in the LIBRARY,
      SchoferFormula.m -- so the symbols could never resolve at all
    ShimuraCurveLattice returns a QuaternionLatticeData record where it used to return 5 values
    ElementOfNorm takes the order and the basis, and returns ONE value where it returned two

⚠ And a fourth, which is the repo's signature trap: `Ldata`Q` holds the Gram matrix over the
**rationals**, while the original line built it over the **integers** (`ChangeRing(Qinv^-1, Z)`).
The two compare EQUAL — `Ldata`Q eq Qint` is `true` — and `lambda_v*Q` then fails with "incompatible
coefficient rings". Equal as values, different as objects. I made that substitution on the strength
of the equality test and had to undo it.

### Now green, and it checks something

    InternalBorcherds: 104 sigma identities, 6 published Wpoly values...Success! 0.270 s

`test_kronecker_sigma` and `test_W` now **return their assertion counts**, and the file asserts the
counts (104 and 6). A caller that only knows "it did not throw" cannot tell a thorough run from an
empty one. Negative controls RUN, not assumed:

    perturb a published Wpoly value (w22)        -> RED
    silently empty the sigma loop (kappas := []) -> RED   (the COUNT catches this, not the asserts)
    restored                                     -> GREEN

The six `Wpoly_scaled` values are Yang's published ones at `d = -4` and `d = -3`; they had not been
checked by any run this repo has a record of.

### ⚠ `test_bp_KY` IS DELIBERATELY NOT WIRED IN — and its failure is now LOCALISED

It is a probe, not a test: its assertion is commented out in the body and it returns a list of
mismatches. Once the import was repaired so it could run at all:

    test_bp_KY(10)  ->  28 mismatches
    test_bp_KY(20)  -> 118 mismatches

⚠ **All 118 sit at `p = 2` with `mu = 0`.** Every odd prime agrees, and so does `p = 2` at
`mu = 1/2`. So the discrepancy is confined to the **`Wpoly2` branch** — worth recording because the
file's own header guesses at "a sqrtp factor that I am missing", and a missing `sqrtp` would have
moved the odd primes too. Left as a probe until that is understood; do not assert on it.

### Import paths in this repo are not what they look like

`import "X.m"` from a file that `run_tests.m` eval's resolves against the **working directory**, but
from a file reached through another `import` it resolves against the **importing file's directory**.
So `tests/InternalBorcherds.m` needs `"SchoferFormula.m"` when run as a test and would need
`"../SchoferFormula.m"` if anything ever imported it. Nothing does, and top-level statements are
illegal in a file used as a package anyway — which is a second reason those three calls could not
simply have been added while the file was being imported somewhere.

## Handoff — 2026-09-15 (night, later) — `find_t` DID NOT PROVE ANYTHING, AND NOW IT SAYS SO

`6_131` and `6_137` were recorded as screen failures on `assert success eq 0` inside `find_t`
(`BorcherdsForms.m`). **Neither is obstructed.** The assert was reading "Magma's integer LP gave up"
as "the problem is infeasible", and those are not the same statement.

Full account, with every measurement: `vvdata/weyl-campaign/find-t-lp-solver-gives-up.md`.

### The witness

The family is `M = 12N` (`D = 6`). The three solved neighbours return the SAME eta-exponent vector
with only the pole order scaling — `N = 73, 89, 107` give `k = 144, 176, 212 = 2N-2`. Extrapolating
to `N = 131, 137` and substituting into `find_t`'s own constraint blocks:

    N 131  M 1572  k 260 : eq true  ge true  le true  ge2 true  ==> FEASIBLE
    N 137  M 1644  k 272 : eq true  ge true  le true  ge2 true  ==> FEASIBLE

⚠ **The hypothesis I started from is REFUTED by its own witness.** I expected the hard-coded
`SetLowerBound(LP, n, -1000)` to be too tight at large `M`. The witness's smallest entry is **-260**
— the bound was never binding. Drafting the check before touching the code is what caught it.

### The failure is sporadic in BOTH directions

    M = 1572  bound -1000                            -> gives up (success 25)
    M = 1572  bound -261 -300 -500 -800 -1500 -3000  -> success, k = 260 every time
    M = 1284  bound -1024                            -> gives up, while -1000 -2000 -5000 succeed
    M = 732, 948, 1068, 1308, 1356   bound -5000     -> gives up, while -1000 -1024 -2000 succeed

Tighter works, looser works, looser-still fails, and nothing tracks `M`. ⇒ **`success != 0` carries
no mathematical information**, and no single bound is safe.

### The fix and its guard

`find_t` retries over `[-1000, -1024, -2000, -800, -5000, -20000]` and errors only if every one
gives up, saying in the message that this is not a proof of infeasibility. `-1000` is tried FIRST,
so the production path for every base that already worked is bit-for-bit unchanged — verified, not
assumed (`M = 876/1068/1284` reproduce `k = 144/176/212`). A solution whose minimum entry **equals**
the bound is rejected and the next bound tried, since a binding bound may have truncated the search.

`tests/FindT.m` (51 s, 5 bases) checks the optimum AND verifies each `t` against the constraint
blocks — a `k` assertion alone would only say the solver returned what it returned last time.
Negative controls RUN, not assumed: single bound `[-1000]` → red; `ETA` perturbed → red; restored →
green. Suite: 11 files touching `BorcherdsForms` re-run green, incl. `X0_15_1 X0_35_1 X0_6_11`.

⚠ **Not a two-base footnote.** Of the 78 never-screened even-`D` targets at `#div(M) <= 20`, **24
are `6_N` with `M = 12N >= 1788`** — the same family, past where the old code first aborted.

### ✅ THE TIMEOUT DIAGNOSIS IS CONFIRMED BY MEASUREMENT, and it cost a CLEAR base

The two re-screens that have finished both spent longer in `WeaklyHolomorphicBasis` than the old
`timeout 1800` cap allowed, which is why their first attempt left a header and nothing else:

    6_89    WHB 2090 s, wall 2206 s   ladder 0 0 0   CLEAR       -> pipeline launched
    178_3   WHB 2340 s, wall 2500 s   ladder 1 1 1   OBSTRUCTED

So the cap was not a marginal call: **`6_89` is a clear base that the harness threw away.** Eight
re-screens still running.

### ⚠ TWO SUITE FILES VERIFY NOTHING — found in passing, NOT fixed

`tests/BorcherdsProducts.m` (0.010 s) and `tests/InternalBorcherds.m` (0.000 s) have **zero
top-level statements**. Both are pure helper libraries that lack the `_` prefix the runner uses to
exclude helpers, so they run as tests and assert nothing.

`BorcherdsProducts.m` is a genuine library (35 test files import it) and is harmless apart from
inflating the file count. **`InternalBorcherds.m` is not**: it defines `test_kronecker_sigma`,
`test_bp_KY` and `test_W` and **calls none of them**, and nothing outside the file calls them either
(the four apparent hits are `test_Whittaker2`/`test_WeilRepresentation`, prefix collisions —
checked). Compare `tests/Whittaker2.m` and `tests/WeilRepresentation.m`, which define a procedure
and then invoke it on the last line. This looks like a plain omission, but adding the three calls
may turn the file red, so it is a separate piece of work.

### ⚠ A TRAP I WALKED INTO — don't repeat it

I edited `BorcherdsForms.m` while 8 test jobs were running from the same tree. `AttachSpec` loads
packages on demand, so those runs could have mixed code versions; their results were discarded and
the whole set re-run on a frozen tree. This is the local-checkout form of "never `git pull` a clone
that has jobs running from it". **There is also another Claude session running Magma on this Mac**,
so `pgrep magma` counts are not yours alone.

## Handoff — 2026-09-15 (night) — THE COLLECTION, AND TWO BATCHES LOST TO HARNESS, NOT MATHS

Collection pass over lovelace. Nothing here is a new mathematical fact; the value is that **two
batches that read as negative results were not results at all**, and both are back in flight.

### ✅ THE SUITE IS GREEN ON 2.29-10, `X0_15_1` INCLUDED

`~/shimura/suiteout`: 84 logs, **75 report `Success!`, 0 failures**. `X0_15_1` passes on lovelace
now — so the genus-0 conic-class fix (`f3d98a5`) is confirmed on the very version that exposed
Magma#125, not just on the Mac. Supersedes "81 of 84 green, the only failure is `X0_15_1`".

⚠ The tree's `HEAD` is 9 commits behind `origin/main` and does **not** contain the fix — it carries
it as working-tree modifications. Checking `HEAD` alone would have given the wrong answer about
which code ran; `git diff origin/main -- <the fixed files>` is the check that means something.

The 9 logs without `Success!` are all accounted for and none is a failure of the library:

    X0_10_23 X0_6_29 X0_6_31 X0_6_37   still running
    run_filters                        still running (FilterByGeneralizedComplicatedFixedPoints
                                       3103 s, FilterByTrace 2789 s -- slow, but progressing)
    _gyinvol _gyinvol_crv              diagnostics, they print and never assert
    _basesweep _rebaselever            ⚠ PARAMETRISED HELPERS.  They need `Dd`/`Nn` and a blanket
                                       per-file sweep runs them without arguments, so they fail
                                       with "Identifier 'Dd' has not been declared".  A harness
                                       artefact of the sweep, NOT a defect -- do not chase it.

### ⚠ `6_73` NEVER FAILED — it was SIGTERM'd, and it is screen-CLEAR

`bk3/DRIVER.log` reads `EXIT 143 6_73`, seven minutes after its own screen finished, i.e. at the
moment the batch was relaunched. `143` is SIGTERM: the process was killed from outside, and its log
stops at "Computing Borcherds forms...". Its ladder is `1 0 0 0` — **clear**. **Relaunched.**

### ⚠ TEN EVEN SCREENS PRODUCED NO VERDICT — a 30-minute timeout, not an obstruction

`6_89 178_3 278_1 298_1 302_1 314_1 326_1 334_1 346_1 358_1` each left a **19-byte log** holding
only its `DEFICIT BASE D N` header. The screen driver caps the run at `timeout 1800`, and the
comparable bases that DID finish spent longer than that in `WeaklyHolomorphicBasis` alone —
`254_1` 2091 s, `262_1` 2402 s, `274_1` 1383 s. So the cap, not the maths, is the likely cause.

**All ten relaunched with no timeout** (`~/shimura/def4.sh` → `~/shimura/defic4/`). ⚠ Recall that a
killed Magma run loses its buffer, so a truncated log is indistinguishable from "never started" —
which is exactly how these read. **`screened-2026-09-14.txt` must not gain an entry for any of them
until a real ladder comes back.**

### ✅ A NEGATIVE CONTROL FOR THE ODD LADDER, 6 FOR 6

`~/shimura/defic/ODD_*.log` ran the **even** screen's criterion (`Ncols - Rank`, swept over `P`) on
six bases that **all build and have committed models**:

    15_1  2 4 7 11 15 20      21_2  2 3 4 5 8 10      39_1  4 5 7 10 14 19
    51_1  6 9 9 15 20         55_1  7 8 10 16 22      57_1  5 7 10 15 19

Every one reads a large and **rising** deficit — 6 of 6 would be called OBSTRUCTED, and 6 of 6 are
wrong. This is the sharpest statement yet of why odd `D` needs `wdef` and its own `m`-ladder: it is
not that the even statistic is noisy on odd `D`, it is that the even statistic is **anti**-correlated
with the truth there.

### Jobs stopped (approved, not unilateral)

* **The legacy `34_11` run** — 10 days, **131 GB RSS**, Magma 2.29-9, frozen code. `models_34_11.m`
  was rebuilt on current code and committed in `b80c3c5`, so the job could only produce a
  superseded answer. Killing it returned ~151 GB.
* **The six odd screens** (`141_1 145_1 91_1 55_2 65_2 143_1`) at 13 h with no `wdef 0`. Recorded as
  **not cleared — no verdict**, which per the `21_2` refutation is all an odd non-clear ever means.
  Ladders as far as they got: `141_1` 7→4, `91_1` 5→2, `145_1` flat 8, `65_2` flat 3, `55_2` flat 1,
  `143_1` printed nothing but its header in 13 h.
* **`_gyinvol_crv`** — 10 h 49 m stuck at `26_3`, where the by-construction route declined and it
  fell back to `IsIsomorphic` on a paired CRV presentation. That is the documented hang, not a slow
  test; it will not finish.

### ⚠ `21_4` SHOULD NEVER HAVE BEEN SCREENED — `N = 4` is not squarefree

It died at `BorcherdsForms.m:55` on the squarefree-`N` assertion, the known method boundary. Drop it
from the odd list; the filter is in `PLAN.md` step 3 and was not applied when that list was built.

### Provenance check that DID come back clean

`deficit.m` and `genmodels.m` live at the ROOT of every lovelace tree but are committed only under
`vvdata/weyl-campaign/` on the campaign branch — precisely the shape that drifts silently. Hashed
all four trees against the committed versions: `scq-current`, `scq-0914b` and `scq-suite` are
**byte-identical** to `origin/m0-theta-campaign`. Only the old `ShimuraCurveALQuotients` tree has a
stale `deficit.m` — and that is the tree `def.sh` points at, so **do not screen with `def.sh`**;
`def4.sh` runs from `scq-current`.

### In flight at the end of this pass

    pipeline (bk3)   14_37  62_7  6_107  6_113   (~13 h 45 m)  +  6_73  (fresh)
    screens (defic4) 6_89 178_3 278_1 298_1 302_1 314_1 326_1 334_1 346_1 358_1
    suite            4 X0_* files and run_filters

## Handoff — 2026-09-15 (later) — MAGMA #125 FILED; `X0_15_1` EXPLAINED AND FIXED

### ⚠ A RED `X0_*` MAY BE THE MAGMA VERSION. CHECK `GetVersion()` FIRST.

`tests/X0_15_1.m` was red on lovelace and green on the Mac with a **byte-identical model**. Cause
found and filed: **[Magma-Maths/Magma#125](https://github.com/Magma-Maths/Magma/issues/125)**.

For a genus-0 `CrvHyp` given by a **degree-1** model, `IsIsomorphic` returns `false` when the leading
coefficients differ by a NON-SQUARE. Witness, verified verbatim on both versions:

    C1 := HyperellipticCurve(x);
    C2, phi := Transformation(C1, [2,0,0,1], 1, P!0);   //  y^2 = 1/2*x
    IsIsomorphism(phi);      // true   on BOTH
    IsIsomorphic(C1, C2);    // V2.29-7 true,  V2.29-10 FALSE

Magma CONSTRUCTS the isomorphism, certifies it, then denies one exists. `true` is correct: the
handbook defines the notion as "a matrix T and a scalar e ... that induce `y^2 = f1 |-> y^2 = f2`",
and `T = diag(a,1)` with scalar `e` sends `y^2 = x` to `y^2 = (a/e^2)x`, covering all of `Q^*`.
Nothing in the docs excludes genus 0 from `CrvHyp` (`HyperellipticCurve(f,h)` promises "the
nonsingular hyperelliptic curve", and `HyperellipticCurve(C::CrvCon)` converts a conic by design) --
that was the main way this could have turned out NOT to be a bug, so it was checked before filing.
Scope, measured: **degree 1 only**; degree 2 (conic class) and degree 3 (genus-1 twist) agree across
versions. Machines: **lovelace runs 2.29-10, the Mac 2.29-7.**

**FIXED IN-REPO** by never asking `IsIsomorphic` about a genus-0 pair: `tests/BorcherdsProducts.m`
decides those by CONIC CLASS -- `RamifiedPrimes` of the quaternion algebra, the same invariant
`tests/ConicClasses.m` uses; a degree-<=1 model is split, so mixed degrees compare correctly. A
genus-0 cover carrying `ws_data` now ERRORS rather than silently skipping, since that branch
exhibits no map to conjugate involutions by.

⚠ Two dead ends worth not repeating: the first guess, "isomorphic as curves but not as hyperelliptic
curves", is WRONG here -- these are isomorphic in BOTH senses. And on the genus-1 key of the same
test, where that distinction WOULD have bitten, matching `(I,J)` invariants does not settle it
(necessary, not sufficient -- quartics with the same Jacobian can be inequivalent torsors); the
honest check is to exhibit the map, which is `x -> 9x, y -> 108y`, exactly the `scales` the test
file already records.

### WHAT IS STILL RUNNING ON LOVELACE — collect these first

Trees: `scq-current` (legacy jobs), `scq-0914b` (odd screens), `scq-suite` (tests; HAS the
`ScaleForSchofer` fix and the genus-0 test fix). ⚠ Never `git pull` a tree with jobs running from it.

* **`bk3` pipeline** — `34_11` and `74_5` DONE and committed; `134_3` died on CM supply; `6_109`
  FAILED (obstructed, see the false-clear section). Still running: `14_37 6_107 6_113 62_7 6_73`,
  with `62_7` already at "Computing equations of covers".
* **`oddscr` odd screens** — `33_1` clear, `69_1` clear (both BUILT since), `33_2` "obstructed"
  (= NOT CLEARED, nothing more). Still running: `141_1 143_1 145_1 21_4 55_2 65_2 91_1`.
* **`defic` even screens** — 6 still running; ladders go in
  `obstructed-rerun-2026-09-10/screened-2026-09-14.txt` on the campaign branch.
* **`suiteout`** — the parallel test run. ⚠ `onetest.sh` hardcodes that output directory, so a
  second sweep OVERWRITES the first; the `X0_*` re-run after the genus-0 fix landed there.

### Scoreboard for this stretch

    new models      21_1  33_1  69_1  34_11  74_5      (33_1 and 21_1 match published equations)
    oracle          Gonzalez-Rotger 43 -> 47 comparisons, 0 bases without a usable W=[1] model
    new coverage    the second hauptmodul row, 251 assertions, previously unchecked
    obstructed      73 known, a LOWER BOUND
    upstream        Magma#125 filed

## Handoff — 2026-09-15 — THE SCALE FIX CORRUPTED NOTHING, AND TWO MORE MODELS LANDED

### ✅ AUDIT CLOSED: no committed model was corrupted by the `ScaleForSchofer` bug

The fix changes CM values at `d = -4` for ODD `D*N`. `33_1`/`69_1` announced themselves by crashing,
but a base where the doubled scale produced a merely WRONG value would have built silently and
passed every internal check — this repo's canonical failure mode. So every odd-`D*N` committed model
was checked, and the screen is free: `d = -4` can only enter a table if the base HAS a disc `-4` CM
point, and `NumberOfOptimalEmbeddings` decides that in closed form.

    14 odd-D*N models
     8  n(-4) = 0  -- PROVABLY unaffected, no disc -4 CM point exists
        (111_1 15_1 35_1 39_1 51_1 55_1 65_1 87_1)
     3  21_1 33_1 69_1  -- built today WITH the fix, each validated
     1  57_1            -- its d = -4 is a POLE (published table), so no Schofer value there
     2  77_1 93_1       -- REGENERATED under the fixed code: BYTE-IDENTICAL to what is committed

⇒ **Measured, not assumed, for every odd-`D*N` base.** Nothing to re-derive.

### ✅ TWO MORE MODELS FROM THE SCREEN-CLEARED BATCH

    34_11   10 keys   VerifyModelSet 44/0   neg ctl 5 failures   ModelChecks 44/0
    74_5    12 keys   VerifyModelSet 45/0   neg ctl 5 failures   ModelChecks 45/0

`34_11` is the base that had a 10-day legacy run grinding on frozen code; screened clear, then built
on current code. Neither has an external oracle — `10_61` evidence level.

### The second-row check, final tally

**33 tables, 251 checks, 0 mismatches** (the six expensive tables landed after the earlier commit
said 27/195). Figure corrected in `GuoYangCheck.m`, here and in `PLAN.md`.

### Still open / running

* `bk3`: `14_37 62_7 6_113 6_107` still going; `6_109` FAILED (see the false-clear section) and
  `134_3` died on CM supply.
* Odd screens: `33_2` reads "obstructed", which per the refutation means **not cleared**, nothing
  more. `141_1 143_1 145_1 21_4 55_2 65_2 91_1` still running.
* Suite 81 of 84 green; the only failure is the Magma-version-dependent `X0_15_1`.

## Handoff — 2026-09-14 (night, later) — THE RUNAWAY CLASS IS FIXED, NOT JUST ROOT-CAUSED

The section below this one says the cause was "a huge principal part x a column whose scale
differs", and treats the scale difference as legitimate. **It is not legitimate — the scale itself
was wrong**, and correcting it removes the runaway at the source. `4bfb859`.

### THE BUG: `w_1` counted as an Atkin-Lehner involution

`ScaleForSchofer`'s second Ogg clause was `(d mod 4 eq 0) and ((D*N mod (d div 4)) eq 0)`, with no
lower bound on `m`. At `d = -4` that is `m = 1` — the IDENTITY — and `d div 4 = -1` divides
everything, so **the clause fired at `d = -4` on every base**. Even `D*N` hides it (the first
clause, `(d eq -4) and IsEven(D*N)`, gives the same answer and is right to). Odd `D*N` gets a
Schofer scale at `d = -4` that is **too large by a factor of 2**.

⇒ That is the whole runaway class. The huge common factor `C` cancels in `ReduceTable` (which
subtracts the per-row minimum) everywhere except the one column carrying `2C`.

**Confirmed against a separate code path.** `NumFixedPointsByCMOrder` — the repo's own Ogg
implementation, which REQUIRES `m > 1` — reports disc `-4` fixed by **NOTHING** at odd `D*N`
(`33_1 69_1 21_1 57_1`) and by `w_2` at even `D*N` (`6_1 38_1`).

⚠ **The published tables cannot arbitrate this cell.** Every offline Guo-Yang table with a FINITE
`d = -4` value has even `D*N`; the one odd-`D*N` table with `d = -4` (`57_1`) has a pole there.

### ✅ TWO NEW MODELS, both on the PLAIN RECIPE

    33_1   GonzalezRotger MATCH, Jacobian 33a1 (published).  Oracle 45 -> 47 comparisons, and
           "0 bases without a usable W=[1] model" for the first time.
           neg ctl: the -1 twist gives 528h2 and is REJECTED.  VerifyModelSet 44/0 (neg ctl 4).
           ⚠ Built twice, with and without HMFIT=1: BYTE-IDENTICAL.  It does not need the flag.
    69_1   VerifyModelSet 44/0 (neg ctl 5).  ⚠ No external oracle -- 10_61 evidence level.

`21_1` is unaffected (its `d = -4` is a pole) and still needs `HMFIT`.

### ⚠ `X0_15_1` FAILS ON LOVELACE AND PASSES ON THE MAC — Magma 2.29-10 vs 2.29-7

Controlled three ways: it passes locally WITH the fix, it fails on lovelace on the UNMODIFIED tree
run SERIALLY, and the `15_1` model regenerates **byte-identically** on both machines and matches
what is committed. So the mathematics is reproducible and the divergence is in the test's
isomorphism check. **A red `X0_15_1` is not evidence of a regression — check the Magma version
first.** Regression coverage for `4bfb859` was 76 of 84 files green on lovelace.

### ✅ THE SECOND HAUPTMODUL ROW IS NO LONGER UNCHECKED — 195 new checks

Guo-Yang publish the PRIMARY column only, so the `s~` row — an independently computed Borcherds
form — had never been compared against anything. It does not need new data: `s~ = 1 - s` is Mobius
and the cross-ratio is Mobius-invariant, so the published column pins what that row's cross-ratios
must be (`want` is literally unchanged). Now asserted in `tests/_offline/GuoYangCheck.m`.
Measured before it became an assertion: **33 tables, 251 checks, 0 mismatches**. `NEGCTL=1`
perturbs one non-frame `s~` value and the check catches it, naming the disc.

⚠ **Selection bias, and it limits the conclusion**: a table exists only where Guo-Yang published
one, i.e. a base that builds. A base whose `s~` row is wrong may fail to build and so have no
table — `21_1`, which motivated all this, is exactly that case and is NOT in the set. "All pass"
means **no systematic defect among bases that build**, not "no defect".

### `21_1`'s wrong CM value, pinned exactly

Its `s~(-7) = 9`; the other five discriminants agree it must be **36** (`-1/4 + 5/4`, `-9/16 +
25/16`, `-1/16 + 17/16`, `-25/144 + 169/144`, `1 + 0`, each exactly 1). The pipeline reads
`scale_tilde` FROM that cell, so the bad datum becomes the yardstick, satisfies the relation by
construction, and cannot be flagged — the error surfaces as "the four others are inconsistent".
`HMFIT=1` refits and is externally confirmed there (Jacobian `21a2`).
Fingerprint for whoever chases it: `9 = 3^2`, `36 = 2^2*3^2`, so the LogSum is short by exactly
**`2*Log2`** — and 2 is UNRAMIFIED here (`2 | 21` is false), in fact SPLIT in `Q(sqrt -7)`
(`-7 = 1 mod 8`). Not today's ramified-prime story. Ruled out: another Ogg mis-classification —
on `X_0^21(1)` only `-7` and `-28` are AL-fixed, each by `w_7` alone.

## Handoff — 2026-09-14 (late night) — ⚠ `6_109` IS A CONFIRMED FALSE CLEAR

**A real pipeline run contradicts the screen, and it takes a documented claim with it.**

`6_109` was screened CLEAR (`1 0 0`) and launched in today's `bk3` batch. It ran 58 min on current
code and died with **"Failed to find all Borcherds forms"** — the Borcherds obstruction itself.

⚠⚠ **The claim "`6_109` reads `1 0 0` and BUILDS" was NEVER BACKED BY A RUN.** There is no
`models_6_109.m` anywhere — not committed, not on lovelace. It was an inference (low rung reads 1,
higher rungs read 0, therefore truncation, therefore it builds) that got written down as a "live
proof" and then cited in both `PLAN.md` and `HANDOFF.md` as the justification for the `>= 2 rungs`
rule. Delete that citation wherever it appears.

⇒ **Known obstructed is 73**, and the even screen now has a **false CLEAR** on record — the clear
direction is not airtight either, in either parity.

**The `>= 2 rungs` rule itself still stands, but on a different example.** `146_1` reads
`1 0 0 0 0`, has a committed model AND a published Guo-Yang table, and genuinely builds. Re-ground
the rule on `146_1`; `6_109` is now a counterexample to it, not evidence for it.

### A HYPOTHESIS for the mechanism — NOT established, do not quote as fact

The pipeline floors `min_m` at `-(n_oo + k - 1)` and otherwise takes it from the divisor's own `m`s,
so the pole order it actually uses is `max(floor_pole, max|m| over the divisor)`. Rungs DEEPER than
that are pole orders the pipeline never visits. Measured floors:

    6_109   floor 325   ladder at [325, 357, 401] = 1 0 0     FAILS
    6_107   floor 330   ladder at [330, 362, 406] = 0 0 0     still running, screen clear
    146_1   floor  55   ladder 1 0 0 0 0                      BUILDS

If this is right, `6_109`'s only reachable rung is the floor, which reads **1**, and its two zeros
are unreachable — while at `146_1` the divisor pushes the pole order past the floor into the zeros.
That would mean the rule as written ("a low rung reading 1 is truncation; trust the invariant value
above") **discards the one rung the pipeline is guaranteed to use**. ⚠ Checking this needs the `m`s
the divisor actually supplies at each base; it has NOT been done.

### The `vx_skip` retest is dead again, still with no verdict

`115_1` and `123_1` are gone from `ps` with no models written; `123_1`'s log ends mid-search. No
`vx ge 0` assert was ever hit, so the vx fix is not implicated — they simply did not survive.
`vx_skip` stays in `genmodels.m`.

## Handoff — 2026-09-14 (night) — THE RUNAWAY CLASS IS ROOT-CAUSED

`RationalNumber`'s guard has said **"Cause OPEN"** since 2026-09-13. It is no longer open. The chain
below is measured end to end at `33_1`, with `21_1` as the contrast base, under `RUNAWAY=1` (a new
env-gated instrumentation, inert by default, left in place because it is what found this).

### The chain

1. **The Borcherds form at `33_1` has a 19-DIGIT PRINCIPAL PART.** `c(-m)` reaches
   `3133789104529289709` at `m = 2, 6, 7, 8, 10`, while at every other `m` it is a single digit.
   The same dump at `21_1`, which builds, tops out at **10**. So every Schofer value at `33_1`
   carries a gigantic `Log11` component — call it `C`.
2. **`kappa_11(m)` is IDENTICAL at every discriminant** for exactly those `m` (`m=6: -2Log11`,
   `m=7: -4Log11`, `m=10: -4Log11`, and `m=2,8: 0`). So `C` is a pure COMMON factor, not a
   per-point quantity — it should cancel.
3. **It does not cancel at one column, because `ScaleForSchofer` is not constant.** Measured:

       d = -4    n_d 4  W_size 2  scale -1/2      <-- twice everyone else
       d = -12   n_d 2  W_size 2  scale -1/4
       d = -15, -67, -88, -163    n_d 4  W_size 4  scale -1/4

   `d = -4` gets `W_size` halved by Ogg's condition without `n_d` falling with it, so its value is
   **2C/4** where every other column is `C/4`.
4. **`ReduceTable` subtracts the per-row MINIMUM valuation**, so it removes `C/4` from the whole row
   — clearing every column except `d = -4`, which is left holding `C/4` exactly. Measured before /
   after, row 1: `[0, -2277330272783157992, -1138665136391578996, ...]` becomes
   `[0, -1138665136391578997, -1, -1, -1, 0, -1]`.
5. `RationalNumber` then meets a 19-digit exponent and the guard fires.

⇒ **The runaway needs BOTH a huge principal part AND a discriminant whose Schofer scale differs.**
It is not a precision failure (already refuted, and now explained: exact linear algebra reproduces
byte-identically), and it is not the ramified prime being mis-handled — `11` appears only because
that is where this form's large coefficients happen to pair.

### TWO STORIES REFUTED ALONG THE WAY — both by measurement, before they were written up

* **"`Solution` picks a bad representative from `sol + Kernel`."** Plausible — the kernel is
  19-dimensional at `33_1` — and WRONG: `maxsol` is **2**. The solution vector is tiny; it is the
  ECHELON BASIS that carries the digits. `21_1` has kernels of dimension 10-13 and builds fine.
* **"A degree mismatch at `d = -4`" (the factor 2 looked like a degree-2 point read as degree 1).**
  Refuted: every point is in `pt_list_rat` and `find_degs` returns `1` for all seven. The factor 2
  is `ScaleForSchofer`, not the field of definition.

### ⇒ THE NEXT EXPERIMENT, AND WHY IT IS WELL-POSED

`C` is an ARTIFACT, and we know it is because adding a kernel element changes the form without
changing its divisor — so `C` is not an invariant of the problem. The concrete move is to
**LLL-reduce `sol` against the kernel so as to minimise the resulting FORM's coefficients** (not
`sol`'s, which are already small). If a representative with a small principal part exists, `C`
collapses and step 3 has nothing to leave behind.

⚠ And the judge is already in place: `33_1` has a Gonzalez-Rotger target,
`y^2 = -3x^4 - 10x^2 - 243`, Jac `33a1`. Same for `69_1`/`Log23`, the other member of the class,
which today's odd screen also clears.
## Handoff — 2026-09-14 (evening) — `21_1` BUILDS, AND THE ODD SCREEN EARNS ITS CORRECTION

### ✅ X_0^21(1): a NEW EQUATION, confirmed by the external oracle

The cheapest target in the backlog (`M = 84`) has a model — 4 cover keys, **50 s** — and
`tests/GonzalezRotger.m` arbitrates it: **the `W=[1]` full curve has Jacobian `21a2`, matching the
published equation.** Ours is `-343x^4 + 94x^2 - 7`, theirs `-7x^4 + 94x^2 - 343`: the same curve by
`x -> 1/x`. The oracle now makes **45 comparisons** (was 43): 10 genus-one curves + 15 AL quotients
+ 20 splitness checks.

Evidence, each piece negative-controlled:

    GonzalezRotger   MATCH 21a2      neg ctl: the -1 twist gives 336e4 and the oracle REJECTS it
    VerifyModelSet   44 checks / 0   neg ctl: twisting the [1,3] entry gives 3 failures
    ModelChecks      44 checks / 0

⚠ **Built under `HMFIT=1`, which is still OFF by default and is now validated at ONE base.** Row
added to `data/models/PROVENANCE.md`. This is the first EXTERNAL confirmation that the fit picks the
right datum, which is more than the flag had before — but it is one base.

### THE BUG UNDER IT: a fitted scale's SIGN is gauge, and one consumer was not

`HMFIT=1` advanced `21_1` two stages and then died on
**"y^2 and s have poles in different places"**. The probe (`scratchpad/probe_21_1.m`) says exactly
which object disagreed: the three covers and the SECOND hauptmodul all have their pole at `d = -4`;
the FIRST hauptmodul's row has no `Infinity()` entry at all, because its value there is
**`-Infinity`**. Cause: HMFIT fitted `scale = -4/9` where the default reads `+4/9`, and
`s_new := s/scale` turns `Infinity()` into `-Infinity()`.

**That sign is pure gauge** — the sign criterion quantifies over `eps1, eps2 in {+-1}`, so `u -> -u`
maps `eps1 -> -eps1` and leaves the satisfied set unchanged — and the default path can never produce
a negative scale, because it reads one off the ABSOLUTE-value table. Two consumers depend on that:
`RationalConstraintsOnEquations` finds the pole with `Index(table, Infinity())`, and `s_new[i_st0]`
is FORCED to `+1` (at that index `stilde` vanishes, so the relation reads `eps1*s/scale = 1`).
Fix: keep the positive representative inside the HMFIT block. One line, inside an env-gated
experiment; `15_1` was checked first to confirm the convention (`+Infinity`, signs per-discriminant).

⚠ `HMFIT=1` does **not** unblock `33_1`, the other Gonzalez-Rotger target: it still dies on the
`Log11` runaway. Different cause, still open.

### ✅ THE `wdef` CORRECTION PAYS OFF IMMEDIATELY — and the 72-base census SURVIVES

**False-positive check on the even screen: clean.** `10_71 22_31 38_13 6_101` re-read with the
sharper statistic all hold at `wdef >= 1`. The verdicts called on `Ncols - Rank` stand.

**First odd-`D` screens ever run** (10 launched, `~/shimura/oddscr` on lovelace):

    33_1   clear   wdef 0 at m = -15   ⚠ deficit 1 -- the OLD statistic would have called it OBSTRUCTED
    69_1   clear   wdef 0 at m = -6    ⚠ deficit 3 -- likewise

Both are independently known **not** to be Borcherds-obstructed: each fails downstream in the
runaway class (`Log11`, `Log23`), and this session reproduced `33_1`'s failure directly. So the two
bases where the plain deficit and `wdef` disagree are exactly the two where the answer is already
known — and `wdef` is the one that gets them right. That is the correction's first live test.

### Still running

`bk3` pipeline (7), `defic` screens, `oddscr` (8 more), the `115_1`/`123_1` `vx_skip` retest.
## Handoff — 2026-09-14 (later) — THE ODD-D SCREEN: BUILT, AND ITS PREMISE CORRECTED

Read this before `PLAN.md`'s "one substantial piece of deferred work" — that item is now DONE, but
not in the way it was written, and the difference is the point.

### ⇒ WHAT CHANGED MOST: `Ncols - Rank` IS THE WRONG STATISTIC, AND `wdef` IS THE RIGHT ONE

`BorcherdsForms` takes a new `DeficitScreen := false` parameter (`b78149b`). Setting it reports the
obstruction deficit from INSIDE the intrinsic, which reaches the odd-`D` 0-side block without
copying it and without refactoring the hot path — everything new is behind the flag, so the
production path is untouched (`tests/X0_15_1.m` passes, 15.0 s).

**But supplying the missing 0-side rows does NOT make the even-`D` criterion work on odd `D`.** Two
separate errors in the old framing, both measured:

* **The ladder is over `m`, not over `P`.** On even `D` the deficit does not depend on `m` at all,
  and the diagnostic is its INVARIANCE as the pole order grows. On odd `D` the 0-side contributes a
  row block fixed by `m_choice` while a deeper `P` keeps adding columns, so **the deficit GROWS with
  `P` at fixed `m`** — `15_1` at `m = -3` reads `1 3 6 10 14 19` across `P = 10..266`, and `15_1`
  builds. That growth is what the old "overestimate ~20" was really measuring. The odd ladder runs
  over `m`, each rung read at the shallowest `P` legal for that `m`.
* **Full column rank is sufficient, NOT necessary.** `Ncols - Rank` asks whether *every* vector is
  in the image; the search only ever asks it of a target supported on the CM points' coordinates.
  The screen now also reports **`wdef`** — the deficit restricted to the span of achievable targets
  — and decides on that. The gap is real and decisive: **`55_1` reaches `deficit 3` / `wdef 0` and
  builds**; `21_2` sits at `deficit 2` across its whole ladder and builds.

⚠ `P` is chosen from the 0-side discriminants so `relevant_ds` stays a superset of
`relevant_ds_0_oo`. Without that the `Index()` fill returns 0 — exactly the failure `95_1` hit.

### KNOWN VALUES REPRODUCED BEFORE ANY NEW NUMBER WAS TRUSTED

    38_5    obstructed, deficit 1 at every rung     recorded: deficit 1 at poleord 190   ✓
    34_3    clear, 0                                recorded: 0                          ✓
    146_1   1 then 0                                recorded: 1 0 0 0 0                  ✓
    142_1   wdef = deficit = 1 at every rung        CONFIRMED obstructed by a real run   ✓
    158_1   wdef = deficit = 1 at every rung        CONFIRMED obstructed by a real run   ✓
    15_1    clear (wdef 0 at m = -7),  16 s         odd, builds                          ✓
    55_1    clear (wdef 0 at m = -15)               odd, builds                          ✓

`142_1`/`158_1` matter most: they show the sharper statistic does **not** wrongly clear a base that
genuinely is obstructed.

### ⚠ LIMITS — quoting this outside them is the failure mode

* **On odd `D` the screen is fast only when it CLEARS.** An obstructed verdict needs the whole `m`
  ladder, and the 0-side basis at deep `m` costs (`pole_order = -D0*m`, e.g. **4005** at `15_1`).
  That is the opposite of what a screen wants, so odd-`D` triage should read a clear and stop.
* ⚠⚠ **THE ODD "OBSTRUCTED" VERDICT IS REFUTED — `21_2` IS A CONFIRMED FALSE POSITIVE.** Corrected
  later the same day, after this section first said only that the verdict "has no positive control".
  It is worse than uncontrolled: `21_2` is a **Guo-Yang base with a committed model**, and its screen
  exhausts all 8 rungs of `all_ms` at `wdef >= 2` and prints **obstructed**. ⇒ **On odd `D` the
  screen yields CLEAR verdicts ONLY.** Anything else means "this screen learned nothing".
* The mechanism is the one already written down: **`wdef` over-approximates** the achievable targets
  — it takes the whole span of the CM coordinates, not the specific `div_coeffs` combinations the
  search actually forms — so `wdef >= 1` is one-sided in BOTH parities. Even `D` gets away with it
  empirically (7/7, and `142_1`/`158_1` re-confirm); odd `D` demonstrably does not.
* Combined with the cost asymmetry above, this makes the expensive path the worthless one: the full
  `m` ladder costs the most and yields the only verdict you cannot use. **Time-box odd screens, take
  a clear, stop.**
* Odd builders scored so far: `15_1` clear, `55_1` clear, `39_1` clear, **`21_2` FALSE POSITIVE**;
  `51_1`, `57_1` still running.

### THE SECOND COLLECTION: 12 more screens, 11 new obstructed, known obstructed 61 -> 72

Ladders appended to `obstructed-rerun-2026-09-10/screened-2026-09-14.txt` on the campaign branch.
NEW obstructed: `10_71 10_73 10_79 10_89 14_53 146_3 254_1 262_1 58_13 6_127 6_139`. Cleared: `6_73`.

* **`58_13` reads `2 2 2` — a FIFTH 2-dimensional obstruction space** (with `166_3 22_19 74_7
  10_67`), and the second found by screening rather than by a failed run.
* `262_1` (`4 2 1 1`) has only two invariant rungs — the weakest call in the batch.
* **`6_131` and `6_137` get NO VERDICT**: the screen itself died on `assert success eq 0` at
  `BorcherdsForms.m:114` inside `find_t` (`M = 1572, 1644`). That is the polytope/t-ladder stage —
  an infrastructure limit at large `M`, not a rank fact. Do not record them either way.
* The **5 stumbled-into** obstructed bases are identified as `142_1 158_1 166_1 214_1 6_97` (grep of
  the `bk2` run logs); checked disjoint from the 18 screened and from `bases49.txt`, so
  `49 + 5 + 18 = 72` has no double count. Still a LOWER BOUND.

### ⚠ A CLEARED BASE CAN STILL FAIL — ON CM SUPPLY, WHICH THE SCREEN DOES NOT SEE

8 pipeline runs were launched on the screen-cleared bases (`~/shimura/bk3` on lovelace, default
recipe). **`134_3` died in 68 s: "Could not find enough rational CM points!"** — the screen cleared
it on rank and it failed on the *other* triage axis. ⇒ **The deficit screen predicts the Borcherds
obstruction and nothing else.** Route a screen-cleared base through the CM-supply check too.

### WHAT IS RUNNING (lovelace, `~/shimura/scq-current`; do NOT `git pull` it)

* **`bk3`: 7 pipeline jobs** on cleared bases — `14_37 62_7 34_11 6_107 6_109 6_113 74_5`; four are
  already past Borcherds forms and into the CM-point stage. `134_3` failed (above).
* **`defic`: 11 deficit screens** still going (`178_3 278_1 298_1 302_1 314_1 326_1 334_1 346_1
  358_1 6_89` and friends).
* **The `vx_skip` retest is STILL ALIVE** — `115_1` and `123_1` at ~11 h under `genmodels_novx2.m`,
  no `vx ge 0` assert. ⚠ A previous reading in this session that they had died was wrong: the `ps`
  output had been truncated by a `head`. If either completes, delete `vx_skip` from `genmodels.m`.
* 4 legacy jobs on the FROZEN clone (`34_11` at 10 d, `95_1`, `159_1`, `119_1`). `34_11` now also
  has a fresh run in `bk3` on current code; the legacy one is on `magma-2.29-9` and old source.

## Handoff — 2026-09-14 — THE DEFICIT SCREEN WORKS. Read this first.

**Everything below is committed and pushed on both branches; the invariant prints nothing.**

### ⇒ WHAT CHANGED MOST: obstruction is now CHEAP to detect, so STOP ATTEMPTING BASES BLIND

`vvdata/weyl-campaign/deficit.m` computes `deficit = Ncols(mat) - Rank(ech_basis*mat)` — a rank
comparison depending on NO divisor choice, so it skips the CM points, the field-of-definition work
and the 96-triple search. **Validated 7/7** on bases it was never calibrated on, at **seconds to
minutes** against HOURS for a pipeline run.

    magma -b DD:=<D> NN:=<N> deficit.m        # from a tree with ShimuraQuotients.spec

⚠⚠ **A SINGLE RUNG CARRIES NO INFORMATION. The INVARIANCE across rungs IS the diagnostic.** At a
low pole order the basis has not caught up and the deficit reads high from pure truncation.
Require `>= 2` rungs and a stable value. Live proof: `6_109` reads `1 0 0` and **builds fine**;
`146_1` reads `1` at its floor and `0` above. Both would have been condemned off one number.
`f6b8b55` now guarantees `>= 3` rungs (the old fixed list `[51,102,134,190,266]` gave exactly ONE
when `floor_pole = n0+k-1` exceeded 266).

⚠ **EVEN `D` ONLY, and the tempting shortcut is REFUTED.** The omitted odd-`D` block only
`VerticalJoin`s rows, so the reported deficit is an UPPER bound and a reported 0 *would* be valid —
but measured on six odd-`D` bases that all BUILD, the overestimate is **~20** (`15_1`→20, `39_1`→19,
`51_1`→20, `55_1`→22, `57_1`→19, `21_2`→10). No odd-`D` base will ever read 0. ⇒ The 0-side block
must really be implemented: `BorcherdsForms.m:876-978`, 102 lines, self-contained. **EXTRACT it as a
file-local function and `import` it — do NOT copy it** (its own comment warns the lines "must move
TOGETHER"). ⚠ It refactors the hottest path inside a memoised loop, so it needs a FULL SUITE run
(~4 h). Deliberately deferred.

### FIRST PRODUCTION SCREEN: 7 new obstructed, 8 cleared

`obstructed-rerun-2026-09-10/screened-2026-09-14.txt` has the ladders.
NEW obstructed: `22_31 6_101 274_1 38_13 218_1 226_1 10_67`. Cleared for running:
`134_3 14_37 62_7 34_11 6_107 6_109 6_113 74_5`.
⚠ **`10_67` has deficit 2 — a 2-DIMENSIONAL obstruction space**, the fourth known (with `166_3`,
`22_19`, `74_7`) and the first found by screening rather than a failed run.

⇒ **Known obstructed is 61** (49 recorded + 5 stumbled into + 7 screened) and STILL A LOWER BOUND —
only ~20 of 105 reachable even-`D` targets are screened. **Never quote 49.**

### 15 NEW/CORRECTED MODELS (102 model files, ModelChecks 0 failures)

`38_3 46_3 35_2 51_2 57_2 77_1 58_3 26_7 62_3 46_5 82_3 74_3 86_3 22_17` new, `10_3` CORRECTED.
Every one verified with `VerifyModelSet` **and an individual negative control** (twist one
genus>=1 entry by the non-square `-1`, confirm the check fails).

### A SECOND EXTERNAL ORACLE — `tests/GonzalezRotger.m`, 43 comparisons

Gonzalez-Rotger, arXiv:math/0612732v2 Table 1. **It resolved the `10_3` drift** that `ModelChecks`,
`ConicClasses` and `VerifyModelSet` had all passed for a week: the three committed entries were
CONSISTENTLY wrong, and an internal-consistency test cannot arbitrate consistency. Supplies target
equations for `21_1` (`y^2 = -7x^4+94x^2-343`, Jac `21a2`) and `33_1` (`y^2=-3x^4-10x^2-243`, `33a1`).

### ⚠ WHAT IS STILL RUNNING ON LOVELACE — COLLECT THIS

Tree: **`~/shimura/scq-current`** (a COPY reset to `main`; never `git pull` the legacy clone).
Outputs: models in `~/shimura/bk2/`, deficit logs in `~/shimura/defic/`.

* **13 pipeline jobs** from batches 1-2 — collect any `models_*.m`, verify + negative-control,
  commit. Already-handled bases are listed above.
* **16 deficit screens** — read with the `>= 2 rungs + invariant` rule, append to
  `screened-2026-09-14.txt`.
* **`115_1` and `123_1`: the `vx_skip` RETEST, and it looks promising.** `genmodels.m` hardcodes
  `vx_skip = {95_1,115_1,123_1,129_1}` and `quit`s; the vx defect was FIXED on 09-05 and those four
  were "gated on it", never retested. Run with the skip stripped (`genmodels_novx2.m` on lovelace)
  they are deep into Borcherds forms with **no `vx ge 0` assert** — `123_1` logged
  `BFPOOL pole_order=1845 Zero=true pool=1935 rank=1846 cols=1846`, i.e. it spanned the Zero side
  FULLY, which is exactly what used to blow up. **If either completes, delete `vx_skip`.**
* **4 legacy jobs** (7-9 d, silent logs) on the FROZEN `f87b0ae` clone. ⚠ `34_11` was screened
  **clear** (deficit `0 0 0`) so it is slow, not futile — unlike `69_1`, which was killed after the
  screen showed its failure was reproducible locally in 8 min.

### Other open items

* **`21_1`** (cheapest target, `M=84`): `find_signs_hauptmodul` reads its normalisation off the two
  discs where a value VANISHES, so those satisfy the relation by construction and can never be
  flagged. Five discriminants agree on `scale_tilde = 36`; the pipeline takes **9** from `d = -7`,
  the disc it reads FROM. `HMFIT=1` (env-gated, off) fits it and advances `21_1` TWO stages, to
  "y^2 and s have poles in different places". The GR oracle can now judge the result.
* **The runaway class** (`33_1` Log11, `69_1` Log23): precision REFUTED (byte-identical at 3x
  `Prec`); the prime is RAMIFIED both times. Cause open.
* **Class-constancy** (`M0MultiplierExact`): `55_2` dev `2.95e-15` is precision-INDEPENDENT (same
  at `Prec` 100/200/300), so NOT roundoff; `87_2` dev `2.2e-5` is ~0.2% of scale. Probably two
  different problems sharing one message.
* **CI**: the permission gate now runs ONCE (`acbc5b1`) instead of in ~40 matrix jobs — that flake
  reddened a run on 09-13. Runner loss on the ~90 min jobs (`X0_10_23`, `X0_6_37`,
  `ExternalCMValues`) is NOT fixable from inside a job; rerun it.

## Handoff — 2026-09-13 (evening) — the BACKLOG restarts: 5 new models and a second oracle

The even-correction line was closed (below); the session then switched to the dormant model
backlog. Everything here is committed and pushed on both branches.

### ✅ FIVE NEW MODELS, all VerifyModelSet-clean and each NEGATIVE-CONTROLLED

    base   M     where      keys        VerifyModelSet   negctl (twist one entry by -1)
    38_3   228   local      10 / 15     116 / 0          --
    46_3   276   local      10 / 15     116 / 0          116 / 3
    35_2   280   lovelace   12 / 15     146 / 0          146 / 7
    51_2   408   lovelace   10 / 15     116 / 0          116 / 6
    57_2   456   lovelace   10 / 15     116 / 0          116 / 5

`46_3` and `35_2` were both Tier 0 "RationalNumber crashers" — the crash IS the m0 signature — so
the 2026-08-23 plan's falsifiable prediction now has **two** confirmations, three weeks on.
⚠ None of the five has a Guo-Yang oracle: this is the `10_61` evidence level.

### ✅✅ A SECOND EXTERNAL ORACLE — Gonzalez-Rotger, genus one (`tests/GonzalezRotger.m`)

*"Non-elliptic Shimura curves of genus one"*, JMSJ 58 (2006); arXiv:math/0612732v2, **Table 1 p.8**.
Their eleven genus-one `(D,N)` are **exactly** the eleven our own curve data gives.

⇒ **8 of 8 comparable committed models MATCH the published equation** — `14_1 15_1 34_1 46_1 6_5
6_7 6_13 10_7`, none of which had any external corroboration before. Not comparable: `21_1`/`33_1`
(no model) and **`10_3`, whose `W=[1]` key is EMPTY** — filling it would bring an oracle to the base
with the unresolved `[1,2]` drift, which currently has none. Its Jacobian should be `30a2`.

⚠ Compare the JACOBIAN, not coefficients — the models need only be `Q`-equivalent. ⚠ And NOT via
Magma's `Jacobian()`/`EllipticCurve()`: both want a rational point and these curves have **none** by
construction. A first attempt returned `ERR` on both sides and printed a vacuous `MATCH` for every
base. The test uses the paper's own `I,J` invariants, **self-checks them against the labels the
paper itself states** before comparing ours, and both guards are negative-controlled.

    X_0(21,1) : y^2 = -7x^4 + 94x^2 - 343   Jac 21a2   w_21 = (x,-y), w_7 = (-x,y)
    X_0(33,1) : y^2 = -3x^4 - 10x^2 - 243   Jac 33a1

### THE BACKLOG, RECOUNTED — and the tiers are not ordered the way the plan assumes

⚠ The recorded 377/73/304 does NOT reproduce. Under the natural filter (`D>1`, star curve with
covers) the current data gives **798 targets / 93 done / 705 missing**; the plan's missing-demand
histogram peaks at 9 and stops at 21, ours peaks at 21 and runs to 31. Either a further filter is
undocumented or the curve data was regenerated. **Do not quote a % complete until this is pinned.**
Of the missing, **192 are FRESH and reachable** (squarefree `N`, `#div <= 20`, not obstructed, not
already tried).

* ✅ **`#div >= 24` IS a real wall, CONFIRMED on Normaliz** — see `normaliz-wall-probe.md`. Probe the
  SOLVER, not the pipeline: at a matched `n/M`, `#div=12` solves in 1 s and `#div=24` TIMES OUT at
  600 s. The plan's routing stands, for a new reason (time, not polymake's OOM).
* ✅ **The `#div = 16-20` "marginal" tier is REACHABLE** — 3 of 4 tried built (`35_2 51_2 57_2`);
  `69_2` failed for an unrelated reason. That tier had never been tested on Normaliz.
* ⚠ **The `#div <= 12` "reliable" tier is NOT uniformly reliable** — `85_1` (M=340) died with a
  **SEGFAULT at 122.7 GB** after 98 min. `#div(M)` predicts POLYTOPE cost, not downstream memory.
  (Also corrects the old "Magma dies ~11 GB" figure again: that was a machine, not Magma.)
* ⚠ Non-squarefree `N` (`6_25 14_9 6_49`) hits the known assertion-failed METHOD BOUNDARY. Filter
  it out of any batch; I wasted three launches on it.

### THREE CODE FIXES (`a9f33d5`), full suite 77/77

1. **`SchoferFormula.m`: guard `f2 = 1`.** `IsPrimePower(1)` ERRORS in Magma, so the `Yang_tt`
   branch crashed at any FUNDAMENTAL discriminant. Cannot regress (the path previously hard-errored).
2. **`LogSum.m`: the runaway guard ASSERTED A WRONG CAUSE.** "the LogSum did not converge upstream"
   is REFUTED — `Prec 300` reproduces `69_1`'s coefficient **byte-identically**. Cause is OPEN; the
   `Prec` experiment is spent, do not repeat it. The prime is RAMIFIED (`23 | 69`).
3. **`tests/GuoYangQuotientOracle.m`: count the unbuildable quotients.** It printed "0 skipped"
   three lines below three `CurveQuotient failed` lines (Magma #123). Now reports them.

### `21_1`: the error names the wrong discriminants (`2416ad0`, HMFIT)

`find_signs_hauptmodul` reads its normalisation off the two discs where a value vanishes, so those
two satisfy the relation BY CONSTRUCTION and can never be flagged. At `21_1` the error blames
`[-15,-43,-51,-67]` — but those four, plus `-91` (from `gtsweep`'s wider CM selection) and `-28` up
to sign, all agree on `scale_tilde = 36`, while the pipeline takes **9** from `d = -7`, the disc it
reads FROM. Discrepancy exactly **4**, a perfect square, so not absorbable by the `+-1` signs.

`HMFIT=1` (env-gated, OFF by default) fits the normalisation over all rational CM points and
independently reproduces that: `scale = -4/9, scale_tilde = 36, satisfied at 5 of 6`. `21_1` then
advances TWO stages to a NEW blocker: **"y^2 and s have poles in different places"**.
⚠ HMFIT trusts the majority; it does not prove it. Settle it with the GR oracle, not the fit.

### lovelace

`69_1` KILLED as futile (we reproduced its runaway locally in 8 min). Four legacy jobs (7-9 d,
logs silent 8-9 d) still run against the FROZEN `f87b0ae` clone. New work runs from
**`~/shimura/scq-current`**, a COPIED tree reset to current `main` — never `git pull` the clone the
legacy jobs run from. 17 jobs, load 17/256. Normaliz there is 3.10.2 vs 3.11.1 locally; harmless,
lattice points are canonically sorted on cache read and write.

## Handoff — 2026-09-13 (later) — the post-condition-4 blocker, examined

**Supersedes the "blocker then MOVES to `QuadraticConstraintsOnEquations`" claim in the section
below.** Full account with every count: `vvdata/weyl-campaign/even-correction/QUADCONSTRAINTS.md`
(campaign, `9a949b4`), with `ratfit_compare.py` and `ratfit-cmextra.patch` beside it.

### ⚠ THE STAGE ATTRIBUTION WAS WRONG — read the `require`, not the traceback

    Runtime error in 'QuadraticConstraintsOnEquations':
    Error in Schofer table values at rational points - no solution found!

That `require` (`EquationsCovers.m:68`) tests `kernels[j]`, an ARGUMENT, built one intrinsic earlier
in `RationalConstraintsOnEquations` (`EquationsCovers.m:31`). When it fires the quadratic stage has
done no arithmetic. It also **cannot** be the blocker in principle: its relations are built from
`B[1]` and solved inside `P(B)`, so they select within an existing solution space and can never
create one. At `34_3` it is moot anyway — **`#quad = 0`**, so the stage is a no-op in both runs.

⇒ The real failure is the **rational linear fit**: no `f` of degree `<= 2g+2` matches the rational
CM values. Empty at all 7 cover keys, at FULL COLUMN RANK (baseline: `dimB = 1` at all 7).

### ✅ THE HATCH'S FOUNDING PREMISE HOLDS — and more strongly than "branch locus"

* branch locus preserved: **0 of 42** (rational CM point, key) cells change zero-ness;
* `h = y2_pert / y2_base` is a perfect **6th power at 58 of 58 cells** once the bad primes of
  `D*N = 102` are divided out, and its 6th root depends **only on `s`, not on the cover key**.

⇒ `h = c_key * G(s)^6`, one shared `G`, `c_key` a bad-prime constant (the re-chosen
`find_y2_scales` row scale). `amt = 6` is the exponent — exactly what adding `6*Z(164)` to the
divisor predicts. Predicted first, then measured.

⇒⇒ **`G^6 = (G^3)^2` is a perfect SQUARE, so `y^2 = c*f_old*(G^3)^2` is the SAME CURVE as
`y^2 = c*f_old`.** The even correction preserves the whole cover, not just its branch divisor, and
**the pipeline is discarding a correct answer** — its ansatz is simply too short.

### ⇒ THE BINDING QUANTITY IS A COST, NOT A FIFTH CONDITION

    degree gain    = amt * deg Z(disc)
    the fit needs    2g+5 + amt*deg Z(disc)   RATIONAL CM points

`deg Z(164) >= 2` at `34_3` (`= 1` refuted on an exactly-determined system at the `g=0` keys), so
the demand is `>= 19` against a supply of **10** — and 10 is genuine: `#ds = 19` still yields
`#rat = 10`, the six extra points all quadratic. **The hatch is not obstructed at `34_3`; it is
PRICED OUT.**

⚠ No selection rule so far considers `deg Z(disc)`, and `PROBE_EVEN` actively prefers the
**largest** `|disc|` — hence large class number, hence large `deg Z(disc)`. **That heuristic works
directly against the cost.** The cheapest legal correction is `amt = M` at a discriminant with
`deg Z(disc) = 1`.

### ✅ THE MECHANISM REPLICATES AT `35_1`, AND CONDITION 4 IS DISC-DEPENDENT

`amt 12` / disc 32 — `35_1`'s recorded clearing config — **reproduces on current code**: condition 4
clears, then the run dies in exactly the same place as `34_3`, with `dimB = 0` at both fitted keys
and **full COLUMN rank** at every non-vacuous degree bound. Same error naming the same wrong stage.
⇒ The rational-fit analysis is NOT `34_3`-specific.

⇒⇒ **CONDITION 4 IS DISC-DEPENDENT — the control `PLAN.md` asked for, answered:**

    35_1, amt = 12:   disc  32  ->  condition 4 CLEARS
                      disc 112  ->  condition 4 FAILS (RationalNumber)

Same base, same amount, different discriminant, different outcome. **"The modulus at base X" is
not a well-defined object**; every recorded modulus (`6` at `34_3`, `12` at `35_1`, `| 6` at
`21_2`) is really a `(base, disc)` measurement. ⇒ **RETIRE the modulus-formula hunt as ILL-POSED**,
do not resume it — the six refuted formulas were fitting a quantity that does not exist.

⚠ `35_1`'s `g=2` cover (`9044`, `W=[1,5]`) is **not fitted through this path**; only `9045` (`g=0`)
and `9046` (`g=1`) reach `RationalConstraintsOnEquations`. A demand quoted "at the `g=2` cover" is
a wrong-object number.

⚠ **Supply is only meaningful with the ask attached**: `34_3` returns 10 rational points whether
asked for 7 or 19 (a genuine cap); `35_1` returns 8 asked for 9, and **12** asked for 21.

Affordability with the correction that actually WORKS on each base:

    base   working correction    cost       g=1 cols needed   supply   short by
    34_3   -164, amt 6           6*4 = 24        30            10       ~3x
    35_1    -32, amt 12         12*4 = 48        54            12       ~4.5x

### ⚠⚠ RETRACTED SAME DAY — "minimum cost 24" IS REFUTED AT `35_1`

The section immediately below reports a cost floor of 24 at `34_3` and explicitly warns it is one
base. **`35_1` refuted it within the hour**, which is the fifth time in three days a law fitted
from one base has been REPLACED rather than refined by the second.

    base    cheapest INTEGRALLY SOLVABLE discriminant        cost
    34_3    -164,  deg Z = 4                                  24     nothing legal below 24 (0 of 11)
    35_1    -112,  deg Z = 1, non-ram                          6     legal at 2 of 3 cover keys

⇒ **Legality does NOT require an expensive discriminant.** "Cost >= 24" is a property of `34_3`.

**What SURVIVES is the framework, not the number.** `cost = amt * deg Z(disc)`, and the fit
needing `2g+5 + cost` RATIONAL CM points, comes from the 58/58 sixth-power structure and is
untouched. Also untouched: `deg Z(164) = 4`, so `34_3` itself is still priced out at 31 points
against a supply of 10.

⇒ ⚠⚠ **"THIS REOPENS THE HATCH AT `35_1`" WAS WRONG AND IS WITHDRAWN** (see the section above).
Two errors: `amt = 6` is ILLEGAL at `35_1` (the modulus is 12), so "cost 6" never existed; and at
the legal `amt = 12`, `-112` **fails condition 4** anyway. The cheap discriminant that appeared to
break the floor is not usable at all. ⇒ The replacement for the refuted floor is not a number but a
**COUPLING**: on both bases the discriminant satisfying every condition is an expensive one, and
the cheap ones are filtered — by condition 3 at `34_3`, by condition 4 at `35_1`. ⚠ Thin (2 discs
tested for condition 4 at `35_1`, 1 at `34_3`) — do not promote this to a law either.

⚠ **TWO CAVEATS ON THE REFUTATION ITSELF, both open:**
* **the `35_1` deg Z sweep is UNVALIDATED.** `degz.m`'s self-check is hardcoded to `34_3`'s five
  discriminants, most of which are not valid CM discs for `D = 35`, so it printed `0 of 5` and the
  guard (which only fires at `34_3`) let the sweep run. `deg Z(-112) = 1` rests on the METHOD,
  validated only at the other base. `degz_known.py` (drafted, not yet applied) derives the
  per-base known values from the divisor-degree identity and is what closes this.
* **no discriminant is integrally solvable at ALL THREE `35_1` cover keys** at `amt` 6 or 12 — the
  legal ones cover `9045`/`9046` but not `9044`. `34_3`'s `-164` was legal at all 7. So `35_1`
  trades a cost problem for a COVERAGE problem. ⚠ Note `IntegralSolution` is OFF by default, so
  condition 3 is MEASURED, not enforced — which sits oddly with the record that `amt 12`/disc 32
  cleared `35_1` cleanly, and that tension is unresolved.

⚠ **A counting trap at `35_1` that does not exist at `34_3`:** each key emits MANY `INTSWEEP` /
`PROBEDS` lines (successive pole-order iterations), so line-count != key-count. `costtab.py` counts
lines and is therefore WRONG at this base; dedupe by key before reading it.

### ✅ THE COST IS MEASURED — and `34_3` cannot afford the cheapest LEGAL correction

`degz.m` computes `deg Z(d)` from `FieldsOfDefinitionOfCMPointFast` (the quantity
`replace_column` already uses) and **refuses to print a sweep unless it first reproduces the five
values the divisor-degree identity pins independently** — 5 of 5. That identity is itself new and
worth keeping: `deg f = (multiplicity at -3) · deg Z(3)` reproduces **all 7** baseline polynomial
degrees at `34_3` exactly, the disc `-3` CM point being `s = ∞`.

    deg Z(164) = 4     =>  cost 6*4 = 24, i.e. 31 rational CM points needed, against 10

⚠ `deg Z` is **not** `h(d)/2` or any similar formula: `h(-408) = 4` with `deg Z(408) = 1`, while
`h(-68) = 4` with `deg Z(68) = 2`. Do not fit one; call `degz.m`.

`costtab.py` crosses that against a hoisted `PROBE_INTSWEEP` over all 21 relevant discriminants.
Restricting to those available at **all 9 keys** (a correction must apply at every cover key) and
non-ramified:

    amt 6:   degZ 1 (11,20,24,27) 0/9 | degZ 2 (56,68) 0/9 | degZ 3 (116) 0/9
             degZ 4:  164 -> 9/9,  180 -> 0/9
    amt 12:  56 -> 9/9,  164 -> 9/9,  180 -> 9/9;  every degZ 1 still 0/9

⇒ **`-164` is the UNIQUE usable discriminant at `amt = 6`, and the joint most expensive available.**
The minimum cost over ALL legal `(disc, amt)` pairs is **24**, reached twice by different routes —
`164/6` (`6·4`) and `56/12` (`12·2`). **Nothing below cost 24 is integrally solvable, 0 of 11.**
The failures are condition 3, not condition 2 — every candidate reports `inimage true`.

⇒⇒ **The hatch is STRUCTURALLY UNAFFORDABLE at `34_3`, ~3× over supply** — not a bad choice of
discriminant. A clean negative result with a stated mechanism, replacing an unexplained failure.

⚠⚠ **ONE BASE. Do not promote "cost >= 24" to a law.** `24` coincides with several quantities at
`34_3` and nothing distinguishes them — this is exactly the shape of the four laws refuted on
09-12/13. The honest statement is the measurement. `degz.m` and `costtab.py` run unchanged at
`35_1` and `21_2`, and that is the next experiment.

### ⚠⚠ A THIRD NULL-RUN TRAP, same family as the two below — a knob at a dead call site

The first `CMEXTRA` knob went into `EquationsOfCovers` (the `[4/6]` path). **`genmodels.m` never
calls it** — `AllEquationsAboveCovers` has its own copy of the `num_vals` computation. Both runs
came back byte-identical in the same 207 s, which reads as "more CM points don't help": a clean
refutation of the degree story **from a run that added no points**. Caught because `[4/6]` appeared
in no log, baseline included. ⇒ The knob now prints `CMEXTRA num_vals` unconditionally.

⚠ Also: **the degree sweep is VACUOUS whenever `ncols > nrows`**, which it is at the default
`MaxNum = 7` (`#rat = 6` = `2g+4` for `g=1`, a square system). A kernel appearing above the true
bound there means nothing.

### ⚠ `genmodels.m` output is NOT byte-comparable to `data/models/`

The unperturbed run gives 15 cover keys against the committed 10, and shared entries differ — but
they are **the same curves in a different normalisation** (`W=[1,2,17,34]` fresh = committed/9;
`W=[1,17]` fresh(x) = committed(17x/3); `W=[1,102]` both). `genmodels.m` does not pin `base_label`.
A byte-diff is not a valid reproduction check for this driver.

## Handoff — 2026-09-13 (newest)

Continues 09-12. Everything below is committed AND PUSHED on both branches; the branch-divergence
invariant prints nothing against `origin`.

### ✅ CONDITION 4 — a fourth requirement on any even correction, and it is SATISFIABLE

The instrumentation asked for on 09-12 was built and run. `RationalNumber` fails **iff some prime
carries a NON-INTEGRAL exponent** (`LogSum.m:137`) -- a `LogSm` is a formal `sum_p coeff_p log p`,
so the failure is a fractional EXPONENT, not an irrationality, and naming the prime is the whole
diagnostic. Measured chain:

* the prime is always a RAMIFIED prime of `D` (17 at `D=34`, 5 at `D=35`), never the level prime;
* it is present BEFORE the final rescaling (`scale = -1/4`, denominator 4, wrong prime);
* the principal-part coefficients are INTEGERS, **362 of 362**;
* `Kappa0`'s OWN log-`p` coefficients are fractional -- 235 at denominator 3, 113 at 9, natural at
  `p=17` where `p+1=18`.

⇒ Those fractions are INTRINSIC and are supposed to cancel: a legitimate divisor makes
`sum_m c(-m) kappa_p(m)` an INTEGER, and the perturbation breaks that. So a usable perturbation
needs FOUR conditions, not three:

    1. EVEN                        -- cover unchanged                 (parity survey 28/28)
    2. phi(target) = 0             -- Borcherds' criterion            (always solvable, gcd(phi)=1)
    3. integral solution           -- a form exists at all            (170/834 candidates at 34_3)
    4. sum_m c(-m) kappa_p(m) in Z at every ramified p | D, every CM d   <-- NEW, and not implied by 3

**Condition 4 IS satisfiable** -- at `34_3`, `35_1` and `21_2` there are amounts giving ZERO
non-rational cells, clearing `ValuesAtCMPoints` for the first time since 2026-08-30. ⚠ The blocker
then MOVES to `QuadraticConstraintsOnEquations` ("Schofer table values at rational points -- no
solution found"), a NEW and UNEXAMINED stage that may be CM supply, a known rescue axis.

### ⚠ THE MODULUS LAW IS UNRESOLVED — SIX FORMULAS FITTED AND REFUTED

    D=34 (2*17)   M = 6  EXACTLY    (N=3 and N=7 AGREE -- so it is NOT N-dependent)
    D=35 (5*7)    M = 12 EXACTLY
    D=21 (3*7)    M | 6             (2,4 untestable there -- they fail condition 3)

Refuted: `2N` (at `34_7`), `2*oddpart(p+1)` for `p=17` (at `34_3`), `2(pmin+1)` (at `21_2`),
`2(pmax+1)` (at `34_3`), `2*gcd` (at `35_1`), `2*lcm` (at `34_3`). **Do not propose a seventh.**

⚠ **UNCONTROLLED CONFOUND, possibly making the question mis-posed**: the perturbation DISCRIMINANT
differs across the three bases (164, 32, 16), and nothing shows the modulus is a property of the
BASE rather than of the DISCRIMINANT. ⇒ The next experiment is the disc-dependence control -- vary
the disc at FIXED base -- because it decides whether "the modulus at base X" is even well defined.

⚠ **THE 28-BASE DIVISIBILITY SCREEN IS VOID** (it used `2N`). **The hatch's REACH IS UNKNOWN** --
neither "3-4 of 28" nor the original "49" is supported. Do not quote either.

### ⚠⚠ TWO NULL-RUN TRAPS, SAME CAUSE — always print the perturbed-key count

`PROBE_EVEN_COPRIME` silently excludes the chosen discriminant and perturbs NOTHING while returning
0 cells. Hit twice: `180/+4` at `34_3`, and THREE runs at `21_2` (`N=2`, disc 16 even) that read as
"all amounts clear" -- which would have refuted the `p_min` formula **on no computation at all**.
That filter rests on the non-coprimality hypothesis, which the even-correction README itself records
as REFUTED. ⇒ **Run with `PROBE_EVEN_COPRIME=0`, and count `perturb disc` lines before believing any
verdict.**

### ✅ tests/ConicClasses.m — the first check on a genus-0 twist class (`4f4667a`)

`ModelChecks`' four tests are STRUCTURALLY BLIND to a quadratic twist at genus 0: a conic and its
non-square twist share genus, genus formula, the trivial Weil polynomial AND the point count over
every `F_p` (every smooth conic over a finite field is isotropic, so both are `P^1` with `p+1`
points). **282 of 822 committed entries are genus 0 across 75 files** -- the largest exempted class
in the repo, previously unvalidated, and why `10_3`'s `[1,2]` drift was invisible to CI.

The new test needs no theory: entries under one `(D,N,W)` key are the same quotient over different
bases, so they must share a class in `Br(Q)[2]`, computed as the quaternion algebra `(a, disc)`.
**215 conics, 38 multi-entry keys, all consistent, 0.05 s.** NEGATIVE-CONTROLLED in situ: twisting
one `10_3` entry by `-2` makes it fail and name the key. Suite **76/76**.
⚠ It catches INTERNAL inconsistency only -- it does NOT resolve `10_3`, whose committed entries
agree with each other. Predicting WHICH conic is right is the open arbiter question; `W=[1]` entries
come out RAMIFIED, consistent with `X^D(R) = empty`, which points at Ogg's real-points criterion.

### Runtime on lovelace: where the time goes

Measured scaling (at `51_1`): pool grows LINEARLY in pole order, but `qexps ~ PO^3.3`,
`EchelonForm ~ PO^4.4`, `ech_etas ~ PO^2.5`. Extrapolated to the pole orders the live jobs reach
(`119_1` hit 1309, `111_1` 1665): ~0.6 h per pool build at `PO=1300`, ~1.7 h at `PO=1700`, and
`EchelonForm` OVERTAKES `qexps` around `PO~1200-1600` -- so these runs sit exactly at the crossover.
⚠ The recorded "EchelonForm is only 6 s of 222 s" is from a small pole order and does NOT
extrapolate.

Levers, ranked: (1) build the DEFICIT PREDICTOR -- the deficit is a rank comparison independent of
any divisor choice, so obstruction is computable from the Weil representation WITHOUT the pipeline;
nobody has built it, and it is the only lever that changes the campaign's complexity rather than one
run's. (2) Attack COEFFICIENT GROWTH in the elimination (multi-modular + CRT): the pool is already
100% triangular when sorted by valuation, so `EchelonForm` never searches for pivots -- but "skip
the RREF" is REFUTED, it only relocates 33-digit coefficients downstream. (3) Reuse q-expansions
across the t-ladder (`qexp(t^j f) = qexp(t)^j qexp(f)`), untested, needs higher absolute precision.
(4) Audit `Prec` per base -- precision is `M^2`.
⚠ The scaling table is ONE base; verify the exponents elsewhere before investing.

⚠ Checked: `bfp` (`95_1`) and `vxfix` (`159_1`) DO contain the vx fix `d9b52d0`, so those multi-day
runs are on valid code. All five lovelace jobs still alive; `34_11` at **81.6 GB** with 1450 GB free.

## Handoff — 2026-09-12

### ⇒ LIVE JOBS: `X0_111_1` IS DONE AND PASSED; the five on lovelace are still running

**`X0_111_1` SUCCEEDED on lava in 58595 s (16.3 h).** It survived the 1678-vector pool that the
last handoff flagged as near the ~2000-vector / ~11 GB wall -- no OOM. Result line:

    X0^111(1): 1 curve comparison(s), 0 involution comparison(s), 1/1 expected covers matched;
               4 committed model cover(s) re-derived (0 CRV skipped)

It anchors on the FULL genus-7 curve, which is a stronger anchor than `93_1`'s quotient-only one.
⇒ **Both previously un-re-derived Guo-Yang bases now have passing re-derivation tests.**
⚠ Nothing has been copied off lava; the log is `$HOME/x0_111_1.log` there and its clone is still at
`8dac84c`, i.e. behind `main`.

**lovelace: all five still alive** (`34_11` INTSOL at 8 d 2 h, `95_1`, `159_1`, `69_1`, `119_1`), all
elapsed ~ CPU so all progressing. ⚠ **lovelace is USABLE AGAIN** -- load 40 on 256 cores, down from
324. The last handoff's "do not launch there" no longer holds.

### ⚠⚠ THE `A_m` MAIN LINE IS WITHDRAWN — the hatch is blocked on INTEGRALITY, not on a theorem

Full argument and data: `vvdata/weyl-campaign/even-correction/AM-REASSESSMENT.md` on campaign
(`bb3700e`), which also carries a superseding header on that directory's README.

`PLAN.md` named the `A_m` theorem as the main line because it "unblocks 49 obstructed bases and is
the only item that does". **That justification does not survive its own evidence.**

* `A_m` is DEFINED by `sum_m c(-m) A_m = mult(f)` -- that identity is how the values were solved,
  not a property proved of them. So inserting `A_m` into `Kappa0` adds exactly `mult(f) log N` per
  firing CM point.
* `SchoferFormula.m:1024` ALREADY adds precisely that, and it is `prop:kappa0`'s conclusion verbatim.
* That code was LIVE at `619051a`, the commit both hatch branches were cut from -- so the recorded
  17 non-rational cells were measured WITH the correction applied.
* Re-run on current code: **12 of 18 bad cells are at NON-FIRING discriminants** (`-24 -51 -228
  -408`), where no level term is owed and the `A_m` defect is outside its own stated scope; at the
  6 firing cells the outer term DID fire with nonzero `m0mult` and they are non-rational anyway.
* `rem:gauge` already says the `N | m` support rule is a GAUGE, which explains without any new
  theorem both why `N | m` had to be imposed by hand and why `prop:closedcoef`'s `-a_E` "reproduces
  1 of 13".

**⚠ A TRAP THIS CREATES:** implementing `A_m` inside `Kappa0` WITHOUT removing the outer m=0 term
would DOUBLE-COUNT. `SchoferFormula.m:609` says the correction "actually belongs" there -- true as
bookkeeping, but it is a MOVE, not an ADDITION.

**What actually blocks it:** the perturbed form is NON-INTEGRAL. `m0mult = (1/2)c_eta(0)` goes from
integers (baseline) to quarter-integers (perturbed), so `c_eta(0)` is half-integral; a fractional
multiplier puts a fractional exponent on a prime, which is exactly the `RationalNumber` failure --
and unlike the log-`N` story it predicts failures at firing AND non-firing discriminants, which is
what both runs measure. `IntegralSolution := true` does not rescue it (the perturbed divisor admits
no integral form), **and that verdict is safe only because its control was run: the unperturbed
baseline passes cleanly under the same flag** (0 cells, 12 keys).

⚠⚠ **CORRECTED LATER THE SAME DAY — "blocked on integrality" is an OVER-CLAIM.** The sweep meant
to confirm it refuted it. Integrality is real but PARTIAL: baseline 0 cells, the original heuristic
run 18, and **any integral perturbation exactly 11** (164/+2, 164/+4, 56/+4 all give 11), with
`m0mult` integral again. The residual 11 is INVARIANT in discriminant and amount, so it does not
depend on which even divisor is added. The remaining cause is NOT identified, and the next step is
INSTRUMENTATION (print what `RationalNumber` is handed at one stable bad cell), not a fifth
single-cause story — four have now been refuted by controls.
⚠ A genuine defect WAS found in the old instrumentation: its "prefer the largest |disc|" heuristic
sent 4 of 7 keys to disc 296, which is integrally solvable at NO key, while 164 is solvable at all.
⚠ Useful accident: the `180/+4` run applied no perturbation at all (180 is divisible by N=3 and the
probe requires coprimality), giving a NULL CONTROL that returns 0 cells — so the 11s are caused by
the perturbation, not by the harness.

⇒ The hatch is a SEARCH for an even perturbation that is ALSO integral -- two conditions, not
one -- rather than a wait for an open theorem; but both conditions together are still not enough. ⚠ `A_m`/`b` remains a genuine open question in its
own right (no product of local densities reproduces `b`); it is simply not what blocks the 49.

⚠ Not a byte-reproduction of 2026-08-31: 18 cells vs 17, because the CM evaluation set differs
(`-56/-68` then, `-228/-408` now). Same phenomenon, different sample -- do not read the two counts
as a change in the effect.

⚠ Recorded, not chased: **`mult(f)` is NOT determined by `div(f)`.** The INTSOL and default
baselines pick forms differing by a trivial-divisor kernel element and report different `m0mult`
vectors, while both give 0 bad cells and the same 12 keys.

### WHAT THE OBSTRUCTION ACTUALLY IS — and it is NOT the eta quotients failing to span

Asked directly this session, so it is written down here. **The eta-quotient basis is not the
problem, and "the space is one pole order too small" is REFUTED BY MEASUREMENT**
([[borcherds-obstruction-is-real]], probe at `38_5`):

    bump 0:  poleord 190  rows 164  cols 36  rank 35      deficit 1
    bump 8:  poleord 198  rows 172  cols 38  rank 37      deficit 1

Enlarging the weakly holomorphic space adds forms AND divisor-columns at the same rate, so the
deficit is invariant. Decisively, the annihilator `phi` is **stable under enlargement** -- equal
entry-for-entry on shared discriminants and merely extending to the new ones -- which is the
signature of reading coefficients off a FIXED modular form, not of a truncated basis.

**The actual reason is Borcherds' criterion.** A divisor is the divisor of a Borcherds product iff
it pairs to zero against every form of the obstruction space -- the weight-3/2 cusp forms of the
lattice's dual Weil representation. At `38_5` that space is 1-dimensional with generator `phi`, and
`phi(target) = -22 != 0`, so the requested ramification divisor **is not the divisor of any
Borcherds product**. That is why all 96 triples fail identically: the search is futile by
construction, not unlucky.

Two corollaries worth keeping:
* **The working bases have NO obstruction space at all** (`34_3`, `38_7`: rank = cols, deficit 0 at
  every key, target found on triple 1). The obstruction is absent there, not dodged.
* **It is not a size threshold.** `38_7` is strictly LARGER than `38_5` in every dimension and is
  fully surjective. The cokernel dimension is an ARITHMETIC INVARIANT of the discriminant form, not
  a monotone function of `DN` -- which is why `14_19`/`14_29` work while `14_17`/`14_23`/`14_31`
  fail.

⚠ **DO NOT CONFUSE THIS WITH THE ETA-QUOTIENT EXPLOSION**, which is a different failure mode
(TIMEOUT, not form-failure) and was root-caused to a `D0` bug and FIXED
([[odd-d-etaquotient-explosion]]).

⚠ And note the two axes are independent: the deficit persists with a fully integral basis and under
the integral solve. So the obstruction (is the target in the image at all?) and integrality (does
the preimage contain an integral point?) are DIFFERENT questions -- which is exactly why the even
perturbation can fix the first and break the second.

### ✅ `EquationsByRebase` now runs under a pinned `base_label`; both `model_drift_ok` flags are off

`main` `47ea828`. PLAN item B. The `base_label eq 0` gate is gone and the pin is threaded into the
stage's inner `EquationsAbovePointlessConics`, which had been silently reverting to the default base
(propagation adds NEW bases to `re_eqns`, so that was a real hole).

⚠ **`STAR` IS NOT FORCED TO THE PINNED BASE, AND THE FIRST ATTEMPT THAT DID SO WAS WRONG.** Forcing
`STAR := base_label` is what PLAN item B literally specifies; at `26_3` it runs the stage and fills
NOTHING. Measured `base_count` there: `<8092,1> <8098,1> <8103,3> <8104,3> <8105,7>` -- the pinned
base carries 3 first-level equations, the heuristic picks `8105` with 7, and only `8105` admits a
usable Hauptmodul root.

⚠⚠ **AND THE COMMITTED DATA IS A MIXTURE.** `models_26_3.m`'s `[1,2]` and `[1,13]` were filled by
`7a923ae` on a DEFAULT run -- the old gate was `base_label eq 0`, so that run cannot have been
pinned. So the file carries 13 keys in the `base_label := 8103` presentation plus 2 from an unpinned
rebase, and reproducing it REQUIRES the unpinned `STAR`.

**The check that caught it is the one worth keeping**: `model_drift_ok` tolerates MISSING keys, so
"Success!" with the flag on proved nothing. The flagless run is what failed, with exactly
`[1,13], [1,2] NOT PRODUCED AT ALL`. Both flags are now off and both tests pass on their own merits
(`X0_26_3` 189 s, `X0_10_13` 773 s), with both Guo-Yang oracles green -- and that oracle holds
Guo-Yang's own curves for exactly the two keys the rebase fills. Full suite **75/75, 0 failures**.

### ✅ ModelRegen retargeted from a MEASURED sweep — and it immediately found a real drift

All 38 bases with no `X0_*.m` test and no recorded cost were measured, one magma process each
(so an OOM costs one base, not the batch), capped at 600 s. **Only 12 of 38 finished.**
`CHEAP_BASES` 8 -> 14: `+191 s` for `+40` comparisons, against the old list's 802 s for 114.

    added:      34_1 5.8s->4   46_1 13.3s->4   6_5 22.3s->23
                106_1 38.6s->3  122_1 49.6s->3  118_1 61.4s->3
    left out:   34_3 204s->10   178_1 269s->3   202_1 374s->2        (poor value, not failure)
    capped:     the ENTIRE D=6 large-prime-N family (6_23 .. 6_83, 11 bases), plus
                10_17 10_29 10_31 10_37 10_41 10_53 10_61 14_13 14_19 14_29 22_13 34_5 34_7 38_7 58_5

⚠ **`10_3` DRIFTS, and it is REAL AND PRE-EXISTING** -- 3 non-isomorphic entries at `W=[1,2]`,
0 missing. Confirmed NOT caused by this session's change: it drifts identically with
`EquationsCovers.m` reverted to HEAD. The entries are genus-0 conics and two differ from fresh
output by exactly **-1/2, a NON-SQUARE** -- a quadratic twist, i.e. the unpinned-y2-scale class
([[committed-models-can-be-unreproducible]]), not a lost cover. **Which side is correct is NOT
settled: there is no Guo-Yang oracle at `10_3`.** It is deliberately left OUT of the default list
rather than given an undiagnosed `MR_KNOWN_DRIFT` row -- an undiagnosed row is how `26_3`'s went
inert.

⚠ **`15_4` cannot run here at all**: `N = 4` is not squarefree and `BorcherdsForms.m:55` asserts.
A METHOD BOUNDARY, not drift, and it fails in 0.7 s.

### A near-miss worth recording: I nearly reported a clean suite as truncated

`run_tests.m` prints `Tests failed:` **only when `#failed gt 0`**, so its ABSENCE means zero
failures -- it is NOT the truncation signature that `CLAUDE.md` warns about. My file count also used
a pattern that only matches when `Success!` lands on the same line as the filename, which is false
for every test that prints output first, giving 40 instead of 75. Counting the right thing:
**75 expected, 75 started, 75 `Success!`, 0 `Fail!`.** ⇒ When checking for truncation, count
`Success!` occurrences against the suite's OWN file filter, and read the summary's print condition
before treating its absence as evidence.

## Handoff — 2026-09-10

### ⇒ LIVE JOBS AT HANDOFF TIME — collect these before starting anything

**`X0_111_1` is RUNNING ON LAVA** and will outlive this session. To collect it:

    ssh -J lovelace lava
    tail -50 $HOME/x0_111_1.log            # look for "Success!" / "Fail!"
    pgrep -u $USER magma                   # empty => finished (or died)

It was at `m_idx=3 of 7` after 5.6 h. ⚠ Its pool reached **1678 vectors**; the recorded wall is
Magma dying around **~2000 vectors / ~11 GB**, so if the process is gone with no verdict in the log,
suspect OOM rather than a code fault, and re-run with a smaller `Prec` or on a bigger box.
⚠ Its clone is `$HOME/ShimuraCurveALQuotients` on lava at `8dac84c`; **it is now behind `main`** —
`git fetch && git reset --hard origin/main` BEFORE any new run there, but **NOT while that job is
alive** (`AttachSpec` compiles on demand; see [[never-update-a-clone-with-jobs-running]]).

**FIVE jobs on lovelace** (`ps -u $USER -o pid,etime,time,cmd | grep magma`), all ~100% CPU:
`34_11` (INTSOL=1, 5 d 16 h — PLAN's old item 1), `95_1`, `159_1`, `69_1`, `119_1`.
⚠ Do not `git pull` those checkouts while they run. ⚠ lovelace itself is SATURATED by other users
(load 324/256) — launch new work on **lava**, not there.

Everything else from this session is committed and pushed on both branches; the evidence for the
49/49 refresh is at `vvdata/weyl-campaign/obstructed-rerun-2026-09-10/` on campaign.


### ✅ THE OBSTRUCTED CLASS RE-RUN AGAINST CURRENT CODE: 49 of 49, ZERO FLIPS

Every OBSTRUCTED verdict on record was taken **2026-09-01/02**, and `BorcherdsForms.m` has had six
commits since — including **`d9b52d0` (09-05), "shift the oo-side basis by its own valuation, not
the 0-side n0"**, the vx fix, which is a CORRECTNESS fix to the very stage that raises "Failed to
find all Borcherds forms". So the 49-base figure justifying `A_m`'s priority rested on pre-fix
verdicts. Re-run 2026-09-10 with `spanprobe.m` at `PROBE_BUMP=0`:

    49 bases re-run    49 still obstructed    0 flipped    0 failing for another reason
    runtimes 18 s (38_5) to 1349 s (34_19)

⇒ **The obstruction is not an artifact of the pre-vx-fix code**, and `A_m`'s justification is now
refreshed evidence rather than a stale tally. A prediction recorded before the first run ("still
obstructed, ~60/40") held.
⚠ **The 49 was recovered, not assumed**: harvesting every obstructed verdict across `sweep122`, the
triage waves and the span probes yields EXACTLY 49 distinct bases, independently confirming the
"known 28 + 21 new" figure as the union of recorded verdicts.
⚠ **`38_5` returned in 18 s against 901 s recorded** (~6x, from the q-expansion bootstrap), and its
`pole_order=190 pool=164` reproduces the recorded `poleord 190 rows 164` — the same computation, not
merely another failure. ⚠ **Level does NOT predict cost** here either (18 s to 1349 s, uncorrelated
with M) — the third time that lesson recurred in one day.

### ✅ `X0_93_1` PASSES — 13389 s (3.7 h)

`tests/_offline/X0_93_1.m` (new): 1 external comparison against Guo-Yang's typo-corrected `[1,93]`
plus **3 committed model covers re-derived** (1 CRV skipped by design). This base mattered most
because `models_93_1.m` regenerates ONLY since the vx fix, so a silent regression there would have
left every committed artifact looking fine. Pre-flighted before the run, not after: their
`(3s^3-7s^2-3s-1)(3s^3+s^2-3s-9)` is isomorphic to the committed entry and all three refuted typo
repairs still fail.
`tests/_offline/X0_111_1.m` (new) is running on **lava** — it anchors on the FULL CURVE (genus 7,
hyperelliptic, published), which is stronger than 93_1's quotient-only anchor.
⚠ At m 3 of 7 its pool is **1678 vectors**, near the recorded ~2000-vector / ~11 GB wall where
Magma dies. If it disappears, that is the likely cause, not a code fault.

### ⚠ REMOTE MACHINES: lovelace is saturated, and PLAN's "four blockers" is FIVE

`lovelace` load **324 on 256 cores**, dominated by other users (`xw132`'s `k3rank` since Sep 06) —
the memory entry's warning that "idle is not a durable fact" holds. **Do not launch there.**
`lava` (`ssh -J lovelace lava`) was load 0.04 on 32 cores and is where `111_1` runs; it needed its
own clone, and the committed `polymake/` cache came with it.
⚠ **`PLAN.md` says four blockers; there are FIVE Magma jobs**, and the fifth is
**`34_11` with `INTSOL=1`, 5 d 16 h elapsed at ~100% CPU** — PLAN's old item 1, "the best-value
thing here". All five show elapsed ~ CPU, so they are progressing, not wedged.

### The math: two hypotheses formed, two retracted

Both concerned `A_m`; neither survived contact with the sources, and the record is worth more than
the hypotheses were.

1. **RETRACTED: "I derived the level-prime factor at general m."** Both "results" are already in
   `paper/level-prime-kappa.tex` — Result 1 IS `thm:closed` (`W_{m,N}(1) = (N-1) ord_N(m)`, with
   `cor:support` for the `N | m` vanishing, verified there over 180 checks against my 18), and
   Result 2 is in `sec:open`, which carries the same `alpha_k`/`G(X)` recipe AND the counts. Cause:
   I read the memory's "the next theorem is general `m` at a nonzero isotropic coset" as meaning the
   LEVEL PRIME was open at general `m`; it is not — the sentence means the intersection with the
   `D`-part. **Every number was right; I was wrong about which object was already known.**
   ⇒ **READ THE PAPER BEFORE DERIVING.** Memory entries and code are not a substitute for the
   30-page document in the repo.
2. **CHECKED AND DROPPED BEFORE REPORTING: "the `prop:closedcoef` refutation is a wrong-object
   comparison."** `rem:gauge` does say `-a_E` and `A_m` are two representatives disagreeing
   pointwise while both reproducing the multipliers — but (i) the memory POSTDATES `rem:gauge` by
   two days, (ii) its literal claim "`A_m` does not follow from `prop:closedcoef`" is TRUE, and
   (iii) decisively, `SchoferFormula.m:589` specifies the code needs the log-`N` coefficient of
   `Kappa0`, "nonzero exactly when `N | m`" — the LEVEL-supported object, whose support `cor:support`
   governs, not `-a_E`'s embedding support. **The memory is correct; the hatch is genuinely blocked.**

**What survives of the math:** `prop:closedcoef`, transcribed and evaluated against the repo's own
`Hurwitz`, reproduces `rem:gauge`'s stated values EXACTLY (`0,0,1,2,1,2` at `X_0^15(2)`) — a small
reusable confirmation that the closed form and its implementation agree.


    X0_*.m cover comparisons:   126 hand-written + 337 model-derived over ALL 34 bases
                                (was 126, and NOTHING else); 34 of 34 tests pass
    committed cover keys:       863 across 88 model files
      with a re-derivation test:  343 on 37 bases
      with none:                  520 on 51 bases   <- 476 of them validated ONLY by ModelChecks
    X0_* census:                32 of 34 pass; the 2 failures are MISSING-KEYS-ONLY and diagnosed

### The `X0_*` tests re-derived 41% of the covers. They now re-derive all of them, for free.

`test_AllEquationsAboveCoversSingleCurve` compared ONLY the keys hand-written into `cover_data`,
and `if not is_def then continue` dropped the rest IN SILENCE: 128 hand-written cover_data KEYS against
309 populated model keys, nine bases checking 1 of 15. (⚠ KEYS, not comparisons -- a key holds one
entry per base, so the comparison counts above are the larger multiset figures. Different objects;
do not quote one for the other.)

⇒ **The fix was not transcription.** `AllEquationsAboveCovers` is ALREADY PAID FOR by each test,
and `tests/_offline/ModelRegen.m` already had the right comparison -- it was offline only because
it paid for a SECOND pipeline run per base. So that comparison now runs as a second pass inside
the helper, reusing the run it already did:

* `tests/_modelfile.m` (NEW) -- `ReadModelSet(D,N)`. Isolated in its own file because an `eval`
  inside a procedure that closes over an outer variable segfaults Magma 2.29 (the trap that forces
  ModelChecks.m and ModelRegen.m into top-level form). PROBED in isolation before being built on.
* `tests/BorcherdsProducts.m` -- ModelRegen's MULTISET matching (so a 3 -> 2 loss cannot hide
  behind two committed entries matching one survivor), `<genus, f, h>` handling, CRV entries
  skipped AND COUNTED, a zero-comparison guard, and `model_drift_ok`.

**MEASURED at 6_11: 1 comparison -> 1 + 17, in 119.8 s against a 121.5 s baseline.** The check is
free; the pipeline run was the cost all along.
**NEGATIVE-CONTROLLED:** perturbing one entry to a same-genus DIFFERENT curve and adding a key the
AL group cannot produce makes it fail, naming both causes separately. It could have failed.

⚠ **IT IS A DRIFT CHECK, NOT A VALIDATION.** It says "current code still produces this", not "this
is correct" -- the committed file is what the pipeline itself wrote. Correctness still comes from
ModelChecks (Eichler-Selberg point counts) and the Guo-Yang oracles. The hand-written `cover_data`
entries must NOT be deleted in favour of it: those are Guo-Yang's PUBLISHED equations, and they
are the only entries carrying labelled involutions.

### ⚠ A pinned `base_label` loses EXACTLY the keys `EquationsByRebase` filled

Two tests fail, both MISSING-KEYS-ONLY, zero non-isomorphic anywhere in 34 bases:
`10_13` (`[1,2] [1,5] [1,26]`) and `26_3` (`[1,2] [1,13]`).

`AllEquationsAboveCovers` gates `EquationsByRebase` on `base_label eq 0` (`EquationsCovers.m:1061`),
so a test pinning a non-zero `base_label` cannot reproduce a key the rebase FILLED on a default run.
Both model files say so in their own headers -- `models_26_3.m` even names `[1,2]` and `[1,13]` as
the two that were empty and were "filled, unlocked by EquationsByRebase".

⚠ **THE CONTROL GROUP is what makes this a diagnosis and not an excuse.** `14_3`, `21_2` and `6_17`
also pin a `base_label` and ALL THREE PASS: `14_3`'s empties were fixed by the COPRIME FILTER FLIP,
not the rebase, and the other two never had any. The gate costs the rebase-filled keys and nothing
else.
⚠ **A PREDICTION WRITTEN DOWN BEFORE THE RUN WAS HALF WRONG, AND THAT IS WHY THE RULE IS NOW EXACT.**
It predicted drift at `10_13 14_3 21_2 6_17` and a PASS at `26_3`; the opposite happened for four of
the five. Had the flag been set from the prediction, three tests would have been needlessly
weakened and `26_3`'s real cause never found.

⇒ `model_drift_ok` therefore tolerates **MISSING keys only**. A key the pipeline DOES produce must
still be the committed curve, whatever `base_label` was pinned -- silencing both with one flag would
hide the failure that actually matters.

⇒ **OPEN, and well-evidenced: relax the `base_label eq 0` gate.** `EquationsByRebase` only ever
fills keys that are ALREADY EMPTY, so running it under a pinned `base_label` should not disturb the
pinned presentation's other covers -- and it would make both these tests reproduce their full model.
A pipeline change, so it needs oracle validation, not just a green test.

### The re-derivation gap is bigger than "128 of 309" -- that counted only the tested bases

    88  model files, 863 cover keys
    37  bases have a re-derivation test (CI or offline)  ->  343 keys
    51  bases have NONE                                 ->  520 keys (60%)
    44  bases are validated ONLY by ModelChecks          ->  476 keys (55%)

"Only ModelChecks" is not nothing -- genus, Weil divisibility and Eichler-Selberg point counts,
none of which touch the Borcherds machinery. But it NEVER RUNS THE PIPELINE, so drift there was
invisible to everything in the repo.

⚠ **And ModelRegen's default `CHEAP_BASES` had become PURE DUPLICATION: all nine had an `X0_*`
test.** Retargeted at bases with none. MEASURED PER BASE, because a batch total cannot tell a
3-minute base from a 26-minute one:

    6_1 10_1 14_1 22_1 6_7 6_13   72 comparisons, ~5 min for all six
    10_7            15 keys       26 comparisons, 185 s
    26_5  804 s | 14_11 1475 s | 22_7 1591 s | 65_1 813 s     <- measured, LEFT OUT for cost

⚠ **KEY COUNT DOES NOT PREDICT COST**: `10_7` has 15 keys and costs 185 s; `65_1` has 4 and costs
813 s. `14_43` was killed at 7 h 44 m unfinished.

⚠ **TWO CLAIMS I FIRST WROTE HERE WERE WRONG, both caught by being challenged rather than by a test.**
1. *"The new list is all even `D`, which is a hole."* **It is not a hole**, and "both D parities" from
   the old comment is itself the stale part. **10 of the 14 odd-`D` model bases have an `X0_*` test**
   (`15_1 15_2 21_2 35_1 39_1 51_1 55_1 57_1` in CI, `39_2 87_1` offline) and every such test now
   re-derives EVERY cover key, so odd-`D` model building is well exercised without ModelRegen. And
   **there is no D-parity branch in the code ModelRegen drives**: the only live `IsEven(D)` uses are
   in AL fixed-point code (`ShimuraQuotients.m:842`, `GeneralizedComplicatedFixedPoints.m:125,186`)
   reached from the FILTER/triage pipeline, never from `AllEquationsAboveCovers`;
   `BorcherdsForms.m:9`'s `assert IsEven(D)` is commented out. Parity mattered when ModelRegen was
   the only re-derivation for those bases; it is not any more.
2. *"`65_1` is the only odd `D` among the 51."* It is the only odd `D` among the **44** with neither
   a test nor an oracle mention. Among the **51** without a re-derivation test there are **four**:
   `111_1`, `15_4`, `65_1`, `93_1`. I quoted a figure for one set while naming the other.

⇒ Both were inherited framing rather than measured claims — the first copied from the comment being
replaced, the second a set I had computed earlier for a different purpose. Spend the ModelRegen
budget on COST, not parity.

**RUN END TO END with the new default: 8 of 8 reproduce, 0 drifted, 114 comparisons, 802 s** — inside
the ~20 min the old list cost, and none of the 7 new bases is re-derived anywhere else.

⚠ **AND THAT RUN FALSIFIED MY OWN REASON FOR ONE ENTRY.** `26_3` was included "to keep the
known-drift path exercised"; it reports `OK (16 compared, 1 CRV skipped)` — it REPRODUCES. The drift
`MR_KNOWN_DRIFT` records for `26_3` is entirely in its `W={1}` entry, that entry is a `"CRV"` entry,
and ModelRegen SKIPS every CRV entry. So **the `26_3` row of `MR_KNOWN_DRIFT` is INERT** — it cannot
fire under this code and has probably been inert since CRV skipping was added. With `14_43` out of
the default list, **the known-drift tolerance is now exercised by nothing**. A tolerance that cannot
fire is the mirror image of a check that cannot fail, and it was only caught by running the default
list instead of trusting the reasoning that chose it.

### The oracle's genus-0 branch was a one-bit check

`GuoYangQuotientOracle.m` compared genus-0 quotients by `HasRationalPoint` alone, so ANY two
POINTLESS conics MATCHED -- and 57 of its 170 comparisons (34%) take that branch. Now a real
`IsIsomorphic` on the conics, negative-controlled on `y^2 = -x^2-1` vs `y^2 = -x^2-3` (both
pointless, correctly distinguished).
⚠ **BE HONEST: it changed no verdict.** Every genus-0 quotient at all 20 oracle bases is a POINTED
conic, and pointed conics over Q are all isomorphic to P^1, so the old check was ACCIDENTALLY
equivalent; all 57 still match in the same 4.8 s. It is not equivalent in general -- 73 of the 281
genus-0 entries in `data/models/` ARE pointless (`6_5 6_7 6_83 82_1 93_1`) -- so this guards the
first such oracle base rather than discovering anything.

### ⚠ A LATENT SILENT-CORRUPTION RACE IN CONCURRENT RUNS (found, checked, NOT yet fixed)

`BorcherdsForms.m:180` writes every Normaliz solution to a DETERMINISTIC SHARED PATH
`polymake/polymake_solution_<M>_<n>_<m>`, and reads it back as `FileExists` -> `eval Read`.
`nmzsolve.py` writes that file NON-ATOMICALLY. So two concurrent Magma processes needing the same
UNCACHED triple can have one read a PARTIALLY WRITTEN file -- a valid-looking but truncated point
list, i.e. exactly the "a partially-cached base returns a wrong answer rather than an error" mode
`CLAUDE.md` flags as critical.

⚠ **This session ran up to five Magma processes at once, so it was exposed.** Checked rather than
assumed: NO `polymake_solution_*` was written during any of it and `polymake/nmzsolve.err` does not
exist, so every solve hit the committed cache and no race occurred. The results stand.
✅ **FIXED**: `nmzsolve.py` now writes via a pid-suffixed temp file and `os.replace` (atomic on
POSIX), at both write sites. Validated: byte-IDENTICAL output to the old writer on the same point
list, no temp file left behind, and exercised END TO END by a real `nmzsolve.py` invocation (not
just the cache-read path, which is all a passing test would have touched).
⚠ SHARED-PATH FILE -- **merge it down to the campaign branch.**

### ⚠⚠ AND A SECOND, WORSE ONE FOUND WHILE TESTING THAT: THE SOLUTION CACHE KEY IS INCOMPLETE

The cache key is `(M, n, m)` ONLY. It omits `k`, `sq_disc` and `cuspidal` -- **and all three change
the answer.** Measured at `(M,n,m) = (8,1,0)`, varying only the omitted parameters:

    k24=12 sq_disc=1 cuspidal=0  ->   4 points
    k24=12 sq_disc=0 cuspidal=0  ->   4 DIFFERENT points
    k24=24 sq_disc=1 cuspidal=0  ->  10 points
    k24=12 sq_disc=1 cuspidal=1  ->   0 points

So two call paths asking for the same `(M,n,m)` with different parameters means the second silently
gets the FIRST one's point set -- correct arithmetic about the wrong object, no error anywhere. And
the two call sites DO differ: `HolomorphicEtaQuotients` (`BorcherdsForms.m:194`, live, reached from
line 289) passes `sq_disc := true` pinned at `(M,0,0)`, while the Borcherds path (line 425) takes
the `sq_disc := false` default.

⚠ **LATENT, NOT MATERIALISED** -- measured, not hoped: of the 503 committed solution files **NONE is
a `*_0_0`**, so the `(M,0,0)` site has never cached anything and nothing can be mis-served today.
⚠ **DO NOT WIDEN THE KEY WITHOUT MIGRATING THE CACHE IN THE SAME COMMIT.** All 503 files are named
under the narrow key; widening makes them all invisible, and above the cached frontier a fresh solve
fails SILENTLY. That would turn a latent collision into a guaranteed silent regression everywhere.
Documented at the read site in `BorcherdsForms.m`.

⇒ **HOW IT WAS FOUND, because the method generalises:** regenerating a committed cache file to check
the atomic-write change gave a DIFFERENT point set. The tempting read was "my change broke it". The
actual cause was that the filename does not record the parameters, so I could not reconstruct the
original constraint system -- and that *is* the bug. A mismatch was evidence about the CACHE KEY,
not about the edit, and the writer had already been proven byte-identical independently.

### ⚠ THE BRANCH-DIVERGENCE INVARIANT IS RED, and not in the harmless direction

`CLAUDE.md`'s check -- `git diff origin/main origin/m0-theta-campaign --name-only --
':!vvdata/weyl-campaign/*'` -- "should print nothing but doc files". It prints **53**, including
`EquationsCovers.m`, `SchoferFormula.m`, `run_tests.m`, 15 model files and 30+ test files.

Direction checked, not assumed: **main-only 51 commits, campaign-only 180, and campaign is NOT an
ancestor of main.** So campaign carries real independent work (`rankcheck_gauge.py` on the
`rem:gauge` ambiguity, `cusp7.m`) AND is missing all 51 of main's recent commits -- which include
`EquationsByRebase`, the quotient oracle and the model fills.
✅ **RESOLVED the same day: `main` merged into `m0-theta-campaign`, no conflicts, both pushed.**
The invariant now prints **NOTHING AT ALL** (not even doc files), and main-only commits are **0** —
campaign contains everything on main. Sanity-checked by running from the campaign worktree itself,
which is the only thing that proves the point: `X0_38_1` passes in 8.5 s and the quotient oracle
makes its 170 comparisons there. Campaign keeps its own 182 commits of research work.
⚠ It will drift again the moment `main` moves. **Run the invariant, do not rely on discipline.**
⚠ NOT affected: `tools/regen-model.sh` runs campaign's `genmodels.m` but from the main checkout's
cwd, so `AttachSpec` loads MAIN's packages. Model regeneration is fine.

### A grep that read a fragment and generalised (again)

Building "which bases have an external oracle" by matching `<D, N,` tuples MISSED `93_1`, whose
Guo-Yang check is a bespoke `gy93_*` block at the END of `GuoYangEquations.m`. The count was rebuilt
searching for the model FILENAME and the `D_N` tag too. Caught only because the number contradicted
what the project already knew -- the same failure mode `CLAUDE.md` opens with.

## Handoff — 2026-09-09

    Guo-Yang published equations:   42 reproducible bases
    full curve stored:              38      <- 10_19 and 22_5 regenerated today
    remaining blockers:              4      95_1  119_1  159_1  69_1  -- still running on lovelace
    X0_*.m tests checking involutions:  34 of 34   (was 23 on 09-07)
    Guo-Yang quotient comparisons:     188 over 24 bases, 0 skipped, 0 mismatches

### The big change: the quotient oracle

⚠ **`GuoYangEquations.m` compares only the equations Guo-Yang PRINT** — usually the full curve
alone — so a base with fifteen cover keys got ONE external comparison. But they also print the
INVOLUTIONS, and every quotient follows from those:
`CurveQuotient(AutomorphismGroup(C,[w]))` is `X/W`. That turns one comparison per base into one per
cover key. `tests/GuoYangQuotientOracle.m` (generic, 20 bases) plus four hand-derived files for the
CRV bases now make **188 comparisons in a few seconds**.

**Three errors in Guo-Yang's tables are now determined**, each by evidence rather than preference:
* `93_1`: `-3t` is a typo for `-3s` (confirmed by the journal version).
* `14_5`: the table's `w_35` sign is wrong; **their own Example 36** has it right.
* `10_13`: the table **SWAPS `w_10` and `w_13`** — settled by Ogg's fixed-point rule, with the
  clincher internal to their paper: **their own CM table** puts disc `-52` at Hauptmodul `0` and
  `-40` at infinity, contradicting their involution table and agreeing with our pipeline.

### The pipeline change: `EquationsByRebase`

An EMPTY cover key is often a **Hauptmodul normalisation artefact, not an obstruction**. The
pipeline builds a genus-`g` curve only as a FIBRE PRODUCT, needing degree exactly `g+1` over a
shared base; which degree a quotient has depends on whether infinity is a branch point, which is
ours to choose. `t -> r + 1/u` at a RATIONAL ROOT fixes the degree profile. At `22_5` this
reproduces Guo-Yang's degree-12 polynomial VERBATIM.

Wired in as the last stage of `AllEquationsAboveCovers`; it is a **no-op unless some cover is
empty**, and only ADOPTS keys that were empty. The `ws` transport works because the rebase is
LINEAR on the weighted ambient: `psi = (r*x + z, y, x)`.

⚠ **Cost**: bases WITH empty covers get slower (`X0_10_11` 471 s, `X0_10_13` 872 s). Bases without
pay nothing.

### Traps that cost real time today, all now guarded

* **`run_tests.m` was SILENTLY TRUNCATING the suite.** It globbed every `tests/*.m`, including
  helpers; several end with `exit;`, which kills Magma. 72 of 79 files reported and no summary was
  printed — a truncated run looks like a clean one. Now mirrors the CI matrix. **Check the file
  count and that `Tests failed:` is present.**
* **A model entry may be `<genus, f, h>`, meaning `y^2 + h*y = f`.** Dropping `h` gives a DIFFERENT
  curve of the same genus; 9 entries across 7 files have one. This produced a false "defect" report
  against `models_87_1`, retracted.
* **An incomplete oracle cannot refute anything.** `GuoYangQuotients_10_19.m` was missing `w_10`
  and `w_95`, so "matches no quotient" really meant "matches none of the ones I computed" — that
  produced a wrong retraction of the rebase lever, since re-corrected.
* ⇒ **All three of the day's wrong verdicts were REFUTATIONS**, each correct about its arithmetic
  and wrong about its object. **A failing check needs its object verified as much as a passing one.**
* **A SKIP is a silent gap**: 13 oracle comparisons were being skipped behind a green summary line.
  All recovered; the cause was on the Guo-Yang side (`CurveQuotient` returns a plain `Crv`).
* **`tools/regen-model.sh`'s flag table had gone stale on both rows** — `CMNONCOPRIME` is a dead
  name (the code reads `CMCOPRIME`), and `Y2TWIST` would have produced models differing from the
  committed files and been read as drift. Table now empty; `51_1` and `22_5` verify IDENTICAL.

### Filed upstream

**[Magma-Maths/Magma#123](https://github.com/Magma-Maths/Magma/issues/123)**: `AutomorphismGroup`/
`CurveQuotient` fail for curves in weighted projective (toric) ambients — `IdentityMap` returns a
`TorMap`, not a `MapAutSch`. It blocks the oracle on exactly the `CRV` paired presentations, which
is why `10_19`, `22_5`, `10_13` and `26_3` each need a hand-derived oracle file.

### Later the same day — the fill, and an audit that found a systemic gap

* **ALL 18 remaining EMPTY cover keys filled**, across `6_29 6_31 6_37 10_11 10_13 10_23 14_5
  26_3`. **0 empty cover keys remain: 347 of 347 populated across 38 Guo-Yang bases.** Every filled
  key was checked against the quotient oracle **in a scratch directory BEFORE installing** — that
  ordering is what makes the data trustworthy, and it should not be inverted.
* ⚠ **`14_5` gained two cover keys that never existed in the file** (`[1,5,7,35]`, `[1,7,10,70]`).
  Its AL group has order 8, so there are 15 proper cover keys; the file had 13.
* ⚠ **Existing entries can come back RESCALED BY A SQUARE** (11 did at `10_19`). That is a
  re-presentation, not a regression. Verify entry-by-entry isomorphism; only a MISSING cover is a
  failure.
* **`10_13`'s labelling differs from Guo-Yang by a GROUP AUTOMORPHISM**, and only half is proven.
  Ours differs by `5 <-> 26` AND `10 <-> 13`, fixing `2, 65, 130`; the map is multiplicative so the
  swaps stand or fall together. `10 <-> 13` is PROVEN by fixed points with their own CM table as
  clincher. `5 <-> 26` CANNOT be: both quotients are genus 2 and Riemann-Hurwitz forces `r = 0`, so
  both involutions are FIXED-POINT FREE and Ogg's rule says nothing. We adopt our labelling for
  both; the second half is **inferred by consistency, not established**.
* ⚠ **A SECOND SILENT TRUNCATION, pre-existing**: `tests/test_weil_polynomial.m` ended with
  `quit;`, which kills Magma since `run_tests.m` evals every test in one process. It sorts
  second-to-last, so `trace_formula.m` never ran locally — which is probably why it was believed
  deliberately skipped. It is not slow: **2.4 s**. Fixed.
* ✅ **`tests/_offline/X0_87_1.m`'s long-standing failure DIAGNOSED AND FIXED**, and validated:
  **passes in 4081 s**. The cause was a **DROPPED h-TERM** — the model stores `[1,29]` as
  `<3, f, h>` with `h = x^3+x^2+1`, and the generator emitted only `f`, so the test compared a
  DIFFERENT curve of the SAME GENUS. It stays offline because it is slow, not broken.
* **`X0_206_1` went 1 -> 4 of 4 covers**, including its `h`-bearing `[1,103]`.

⚠⚠ **AND THE AUDIT THAT MATTERS MOST: the `X0_*` tests re-derive only 41% of the covers.**
**128 `cover_data` keys against 309 populated model keys.** Nine bases check 1 of 15
(`10_11 10_13 10_23 6_11 6_17 6_19 6_29 6_31 6_37`) and eleven check 1 of 4. The helper SILENTLY
SKIPS an absent key, so this is invisible unless counted.
⇒ The MODELS are well checked (~190 oracle comparisons over 25 bases against Guo-Yang). What is
thin is the **RE-DERIVATION** claim — that the pipeline reproduces them — which for most bases
rests on ONE cover. Closing it is mechanical but must handle `<genus, f, h>` entries and `CRV`
pairs, both of which have already caused defects, and it costs CI time.

### Still open

* **The 41% re-derivation gap above** is now the largest single opportunity: work the thin tests in
  order of missing covers, starting with the nine at 1-of-15.
* The four lovelace blockers are mid-FIRST-PHASE after ~3 days; weeks away, not days.
* `93_1` and `111_1` still have no `X0_*` re-derivation test (14-20 h per run).

## Older — Handoff 2026-09-07

**Supersedes** the 2026-07-17 handoff about producing cover models, archived as
`HANDOFF_2026-07-17.md`. That task is not dead, but it is gated on the blocker described below.

Everything here is committed and pushed. **`git pull` first — local `main` may be stale.**

**➡ For what to do next, see `PLAN.md`** — five tracks, a do-not list, and the recurring traps.
This file is the record of *what happened*; `PLAN.md` is the record of *what to do*. When the two
disagree about state, this file wins.

## Handoff — 2026-09-07 (the 09-06 section below is still accurate, just earlier)

    Guo-Yang published equations:  42 reproducible bases
    we now have a model for:       38      <- 111_1 recovered
    remaining blockers:             4      95_1  119_1  159_1  69_1   -- ALL RUNNING

* **`111_1` recovered** — 20.2 h, DEFAULT flags, another base the vx fix unblocked. Verified by
  **exact full-curve `IsIsomorphic`, true in 0.05 s**, which is cheap only because its `W={1}` is
  HYPERELLIPTIC. ⚠ Pinnable to one commit (`f87b0ae`, clean tree) — unlike `10_61`/`14_43`.
* **`93_1` and `26_3` upgraded to FULL-CURVE PROOFS** (were quotient-level). `IsIsomorphic` hangs on
  CRV pairs, so the isomorphism is CONSTRUCTED: Mobius map from the hyperelliptic quotient, both
  sides carried by constant squares, then `IsIsomorphism` certifies it. Hundredths of a second.
  ⇒ **The BASE chooses the `V_4`** (assaferan): `26_3` would not match until rebuilt with
  `base_label := 8103`, which is the `V_4` Guo-Yang use. When a CRV pair will not match, try
  another base before concluding anything about the curve.

### Three blind spots removed, each of which immediately exposed a real defect

* **`VerifyModelSet` skips every `CRV` entry** — so 21 paired presentations across 16 files had
  NEVER been checked. `tests/CRVStructure.m` found **5 storing their parent conic twice**
  (reducible schemes, not the genus-1 curves recorded). ROOT CAUSE: at `g = 1` the required degree
  `g+1 = 2` is also a conic's degree, so the conic could fill BOTH roles in the fibre product.
  Fixed; those covers now defer.
* **The `X0_*` helper silently skipped unmatched cover keys** — it could pass while verifying
  NOTHING. It now counts comparisons and errors on zero. That immediately turned CI red, correctly:
  **`X0_10_19` had been passing green in CI while making ZERO comparisons**, at 84 min a run.
* **CI never set `NORMALIZ_BIN`** — so polytope solves failed SILENTLY ("no solutions", not an
  error). Now installs `normaliz-bin` and exports the path. ⚠ Scope was MEASURED: every other
  `X0_*` job reported full coverage, so `10_19` was the only affected test.

### The coprime guard: no evidence it is needed

Full sweep of the 11 `N>1` `X0_*` tests (for `N=1` the filter is provably a no-op, which excludes
19 of 30 rigorously). **8 of 10 pass identically with `CMNONCOPRIME` on and off.** The 2 failures
(`10_13`, `6_17`) are both CRV tests whose PINNED COORDINATE MATRIX breaks under re-presentation —
not correctness. Removing that artifact is what `tests/_crviso.m` does.
⚠ A first sweep appeared to show `10_13` failing; that was MY OWN foreground timeout killing the
sweep's Magma, which then recorded a killed run as a failure. Retracted.

### Process notes

* **`nohup ... &` inside a background call reports completion for the WRAPPER**, not the job.
* **My own foreground timeout killed a background sweep** — a killed run and a failing run are
  indistinguishable in a one-line summary. Capture the error text before believing a regression.
* Bugs of mine caught only because a test could fail: an eager `AutomorphismGroup` (5x slowdown,
  870 s -> 71 min), an inverted conic scalar (`x/rg` for `rg*x`), hardcoded variable order, and
  image polynomials built in the wrong ring. Each was invisible in the first case tried.

## Handoff — 2026-09-07, later (test coverage; supersedes the earlier 09-07 block on these points)

**Re-derivation coverage went 4 → 11 Guo-Yang bases.** Passing an `X0_D_N.m` test IS reproduction,
the stronger claim than `GuoYangEquations.m`'s stored-model comparison. Now covered:
`51_1 55_1 57_1 14_5 14_3 26_3 21_2 15_2 22_3 22_5` in CI, `39_2` offline. 34 `X0_*` tests in CI.
* `X0_21_2` is the **first test that checks a CRV entry** — possible only because the helper now
  CONSTRUCTS those isomorphisms (`tests/_crviso.m`) instead of calling `IsIsomorphic`, which hangs.
* `MR_KNOWN_DRIFT` is down from 5 to **2**: only `14_43` (`INTSOL=1`) and `26_3` (deliberate
  `base_label := 8103`).

**⚠ `Y2TWIST` WAS THE WRONG SUSPECT, and I nearly flipped it on a confounded measurement.**
`PROVENANCE.md` had predicted for two days that making twist selection default was "the right
long-term fix" for `15_2`/`22_3`/`22_5`. Measuring `Y2TWIST=1` against the COMMITTED models showed
large gains — and those were the **coprime flip's**, from hours earlier the same day, because every
committed model predated it. The control is flag-on vs flag-off on the SAME code: run that way all
three bases are IDENTICAL and the deferral path logs zero messages. The selector never fires. The
flip was reverted; the mechanism is kept (it is sound: unique-or-defer) but not defaulted.
⇒ **Compare against a current baseline, never a committed artifact.**

**The real win was already sitting there.** Those three needed NO flag — their committed files were
simply STALE. Regenerated with the plain recipe: `22_5` 3 → 11 populated covers, `15_2` 12 → 15,
`22_3` 13 → 15, nothing lost, `GuoYangEquations` still passing.

**`X0_87_1` is the one known-broken test** and is under diagnosis. Established: the MODEL is fine
(`ModelRegen` reproduces it; `GuoYangEquations` matches its `W={1}`), and the test is well-formed
(expects exactly the model's 4 single-entry keys). So the failure is `assert is_isom`. Leading
hypothesis, which has bitten twice already: `ModelRegen` compares the AGGREGATED model while the
helper iterates EVERY BASE of every cover, so a second base with a different presentation fails only
the helper — fixed at `26_3` and `21_2` with a `base_label`.

## Handoff — 2026-09-06 (this session; supersedes the state notes below)

**Five models produced, and the Guo-Yang denominator was wrong.**

    Guo-Yang published equations:  42 reproducible bases (NOT 43 -- see below)
    we now have a model for:       37
    remaining blockers:             5   95_1  111_1  119_1  159_1  69_1   -- ALL RUNNING

* **`93_1`** — the vx fix unblocked it (default recipe, 14.1 h). It also **settles a typo in their
  table**: their `-3t` is `-3s`, determined by isomorphism from our own model against three refuted
  alternatives, and later **confirmed independently** by the journal version.
* **`26_3`** — recovered with `CMNONCOPRIME=1`, 189 s. Full `V_4` diagram matches; their conic
  `-8x^2-3` comes out coefficient for coefficient.
* **`15_4`** — a FOURTH provenance category: **literature-derived, not pipeline-produced**, and it
  never can be (see below). `a = -1` is confirmed by our point counts, `b = -1` by the full-curve
  trace-formula comparison.
* **`10_61`, `14_43`** — first two models out of the OBSTRUCTED class (41 h, 42 h). ⚠ Neither is a
  Guo-Yang base, so **no external oracle** — `ModelChecks` alone. Weaker evidence; quote it as such.

**⚠ 42, NOT 43.** `15_4` is outside the Guo-Yang method *by the authors' own statement* — their
published Remark 39 says the normalizer of the Eichler order strictly contains the Atkin-Lehner
group there, so the star quotient our pipeline forms is the wrong object. It is not a blocker; it
is out of scope.

**⚠⚠ WE HAD BEEN READING THE SUPERSEDED PAPER.** arXiv:1510.06193 has exactly ONE version (2015).
The paper of record is **Compositio Math. 153 (2017) 1-40**, substantially revised and NOT on
arXiv; our `ShimuraCurves-arxiv.tex` is the arXiv one. The journal fixes `93_1`'s equation and
`39_2`'s involutions, and adds Remarks 38 and 39. PDF is in the user's Dropbox. **Check the
journal, not just v1.** Tu (Pacific J. Math. 269 (2014) — also now in that Dropbox folder, free
from MSP, not on arXiv) confirms `15_4` and covers `26_3`, but supplies nothing for any other
non-squarefree base.

**Speedup shipped:** the q-expansion bootstrap, `qexp(t^j f) = qexp(t)^j qexp(f)`, both sides —
**up to 18.6x** on that step at `pole_order 800`. ⚠ NOT yet shown to help any base end to end.

**The coprime guard is under doubt.** Three bases (`39_2`, `14_3`, `26_3`) produce Guo-Yang-matching
models with it OFF, and `26_3` is the very base whose bad discriminants justify its existence. A
targeted sweep is running. ⚠ Only the 11 `N>1` `X0_*` tests can possibly show anything: for `N=1`,
`gcd(d,1)=1` makes the filter provably a no-op.

### Process lessons this session cost something to learn

* **Never `git pull` a clone that has jobs running from it.** I did, to lovelace, with eight jobs
  running from that directory. `AttachSpec` loads packages ON DEMAND, so a long run can compile
  source that changed under it. `10_61` and `14_43` cannot be pinned to a single commit because of
  it. Launch long runs from a COPIED tree.
* **My own foreground timeout killed a background sweep**, and the sweep recorded the killed run as
  a FAILURE. That produced a false "`X0_10_13` breaks under `CMNONCOPRIME=1`", since retracted. A
  killed run and a failing run are indistinguishable in a one-line summary.
* **`nohup ... &` inside a background call reports "completed" for the WRAPPER**, not the job. Check
  process state; do not trust the notification.
* **Corrections made:** `93_1` was first reported as "34 -> 35" (that is the CM-TABLE count, a
  different set); `EchelonForm` was called negligible when it has the STEEPEST growth (~`PO^4.4`);
  `15_4` was diagnosed as a squarefree-`N` code issue when the authors had stated the real reason
  in a version we had not read.

## ⇒ READ THIS FIRST — 2026-09-04, late

> ### Spend your effort on WHICH OBJECT the claim is about, not on whether the computation is right.

That is the single most useful thing this session produced, and it was learned the hard way. Nearly
every error made on 2026-09-04 had the same shape: **the arithmetic was correct and the object was
wrong.** A rank computed over MONOMIALS when the claim was about FORMS (twice, in two different
ways). A LaTeX parser emitting perfectly valid polynomials from silently truncated input. Three
`grep`s that read a fragment of a file and generalised from it. Every validation in place was of
the form "is this number computed correctly" — **none of them could catch "is this the right
number."**

Two habits did catch things, and are worth keeping:
* **Reproduce a KNOWN value before trusting a new one.** The rank result was only believable
  because the same script had to reproduce the paper's own rank-4 panel and a principal part known
  from `tests/M0Multiplier.m`. Both caught silent data corruption that had produced a
  plausible-looking right answer for the wrong reason.
* **Draft an edit instead of applying it.** The claim "the paper is wrong, fix `rem:gauge`" was
  retracted *while writing the diff*, because writing it forced a close enough read of `sec:exact`
  to notice it says FORMS where I had MONOMIALS. Applying directly would have degraded a correct
  argument in a paper heading for submission.

### State right now

**Three runs were left going on lovelace** (`~/shimura/models/*.genmodels.log`, `M0PROGRESS=1`).
**All three completed their a0 tables**, which is the result they were launched for:

    34_11   13 fallback points of 64      (~8.7 h elapsed)
    10_61   27 fallback points of 64      (~6.7 h)   <- previously DIED at gate 3
    14_43   22 fallback points of 64      (~6.7 h)   <- previously DIED at gate 3

**A completed a0 table means the two-point check never fired, i.e. GATE 3 IS CLEARED at `10_61`
and `14_43`** — the gate that killed both before. That retires the "924 gate-3 failures at 10_61,
zero near-misses" entry and, with it, the claim that `10_61` "is not runnable, it has a real
upstream defect". The fallback rates scale as expected: 1/64 at `15_2` (M=60), 8/64 at `58_5`
(M=580), 22-27/64 at M≈1200 — the threshold is M-scaled, so bigger bases trigger it more.
All three were still in FINAL ASSEMBLY (class-constancy, isotropic agreement, rational snap) when
this was written; **check those logs first — they may have finished, and `10_61`/`14_43` producing
models would change the sweep record's "two runnable candidates of 122" count.** lovelace is shared and busy again (~68 magma
processes, mostly other users) — check `uptime` before adding load. Both branches and lovelace's
clone are clean and in sync; housekeeping list is empty.

### What shipped

* **The per-coset `tau` fix** (`475e72b`) — the only thing that moved the mathematics. `34_11` went
  from failing to passing; validated three ways (`15_2` exact, `58_5` keeps its models, and both
  match the Prop 9.15 closed form 9/9 — `58_5`'s as a PRE-REGISTERED prediction). New
  `M0PROGRESS=1` diagnostic (WriteStderr, because buffered `printf` is lost when a run is killed).
* **`51_1` and `57_1`** — never blocked, just never run. Guo-Yang coverage 32 → **34 of 43**;
  81 model files, `ModelChecks` 8309 checks 0 failures.
* **34 Guo-Yang CM-value tables** as offline tests (`tests/_offline/`, 254 checks, 0 failures) and
  **`tests/GuoYangEquations.m`** comparing 7 committed models to the PUBLISHED equations — the
  first thing in CI that checks our output against the literature. Both negative-controlled.
* **Housekeeping**: branches 4 → 2, 10 archive tags, lovelace pulled up (was 52 behind), the
  t-shift fallback ported to `main` after 9 days missing, 109 stranded Normaliz solves harvested.

### ⚠ What was CORRECTED — do not re-litigate these

* **`rem:gauge` is CORRECT. Do not edit the paper.** It was claimed wrong twice, from measuring the
  wrong object. The `oo`-only model is valid for genuine Borcherds forms (39 of them) and FAILS on
  monomials (residual 2.08), so a rank over monomials says nothing about it.
* **`A_m` genuinely needs new mathematics.** `sec:determined` determines the CONSTANT TERM at
  isotropic cosets ("an indicator, not a phase"); the canonical representative is the SCALAR
  `-a_E`. All-`m` at a nonzero isotropic coset is absent. No shortcut by extraction.
* **`93_1`/`95_1`/`159_1` are the vx class**, not squarefree-`N` — all have `N=1`, which IS
  squarefree. Measured: Magma's own `GalFldFun.m:305 assert vx ge 0`. `genmodels.m`'s `vx_skip`
  list is INCOMPLETE.
* **`gtsweep`'s `FIRE` lever does nothing** (measured on all three bases it claimed to fix).
* **The `cusp7` "scoped implementation task" was falsified** — a monomial-pool coverage gap, not a
  dump bug. MAIN LINE's "provably resolvable" premise is dead with it.
* **Gate 4's "GENUINE 43% violation" was numerical**, not mathematics.

### Where to pick up

`PLAN.md` "Picking this up cold" is current. In short: wait on the three runs; then the 9 remaining
Guo-Yang blockers are correctly classified for the first time (1 structural, 3 vx, 1 nonintegral,
1 non-rational, 2 odd-`D` basis ceiling, 1 open anomaly — `26_3`'s exact `z -> z/(z-1)`
involution). `22_5` and `14_3` need full-curve models GENERATED, not transcribed.

## Update — 2026-09-05: Guo-Yang coverage re-measured, and two code changes

Six commits, `64d9316`..`36ac71e`, all pushed. **Two of them change code and only ONE of the two
is fully validated** — read the status column before building on either.

### ✅ TWO GUO-YANG EQUATIONS RECOVERED — 34 of 43 (`39_2` and `14_3`)

Both were blocked by the **coprime-to-level CM filter**, not by mathematics, and both were filed
under diagnoses that had gone STALE rather than been wrong. `CMNONCOPRIME=1` (env-gated, OFF by
default) unblocks them; the published equations are what make the results believable.

    39_2   filed NONINTEGRAL.   filter on: 3 CM points vs demand 19 -> "not enough points".
           filter off: 24 points, 15 keys 0 empty. W={1} genus-7 hyperelliptic,
           IsIsomorphic to Guo-Yang in 0.06 s. Pinned in tests/GuoYangEquations.m (9 bases).
    14_3   covers under-determined by default, W={1} EMPTY (6 keys, 3 populated).
           filter off: 16 keys, 0 empty. W={1} genus-3 CRV pair, IsIsomorphic in 6817 s.
           Pinned in tests/_offline/GuoYangCurve_14_3.m -- OFFLINE because ~2 h would wreck
           GuoYangEquations.m's ~97 s.

`ModelChecks` passes both independently (82 files, 8767 checks, 0 failures) via trace-formula
point counts rather than the path that produced them.

⚠ **IS EVERYTHING REPRODUCIBLE FROM COMMITTED CODE? NO — and that is partly deliberate.** Of the
34: **24** are verified BY re-derivation (the `X0_D_N.m` tests run the pipeline, so passing IS
reproduction); **10** are stored-model comparisons that never run it. Of those 10, four do NOT
regenerate by default, for two different reasons that should not be conflated:
* `22_3`, `15_2` — ACCIDENTAL: the y2-scale guard (`1768517`) postdates the files. `Y2TWIST=1`
  restores them. This is drift, and it went unseen because every test reads the stored model.
* `39_2`, `14_3` — DELIBERATE: they exist only under `CMNONCOPRIME=1`, which stays off because it
  has no theoretical guarantee. Their justification is the published equation, not regeneration.
  Enabling the flag to make them "reproducible" would trade a documented gap for an undocumented
  risk on every base.
⇒ The target is NOT "everything regenerates by default", it is "every non-reproducing file has a
recorded reason and an independent validation" — which is what `MR_KNOWN_DRIFT` and the model
headers now encode. (`14_5`, `55_1`, `21_2`, `87_1` were still unchecked when this was written.)

### ✅ X_0^39(2) RECOVERED — the first of the two

**The one coverage result of the session.** `39_2` was filed as the NONINTEGRAL malformed-form
base; that was wrong. It is starved by the **coprime-to-level CM filter**: with the filter on it
sees 3 CM points against demand 19 and dies with "Could not find enough points"; with
`CMNONCOPRIME=1` it sees 24 and builds cleanly (15 keys, 0 empty). Its `W={1}` genus-7 curve is
**`IsIsomorphic` to Guo-Yang's published equation** (0.06 s), now pinned in
`tests/GuoYangEquations.m` (9 bases); `ModelChecks` passes it independently (82 files, 8573 checks,
0 failures, via trace-formula point counts rather than the path that produced it).

⚠ **The flag is NOT safe by default and is not enabled.** The `p | gcd(d,N)` local factor has no
live implementation (`kappaminuszero` is dead code), and at `26_3` two non-coprime discriminants
give provably wrong values. What makes `models_39_2.m` trustworthy is the INDEPENDENT ORACLE, not
the flag — **any further base produced this way must clear the same bar before being committed.**
⚠ That file does not regenerate by default; recorded in its header and in `ModelRegen`'s
`MR_KNOWN_DRIFT`, with the reason distinguished from the three y2-guard entries.

`26_3` is the same story but NOT yet closed: a model is produced, its `[1,78]` cover is verified
against Guo-Yang, but the genus-5 `W={1}` comparison had not returned when this was written.

### The Guo-Yang picture, measured rather than inherited

    43   published equations                    (see the counting trap below)
    34   we reproduce today, with a test        (24 pipeline + 9 GuoYangEquations + 14_3 offline)
     9   the gap: 8 with no model, + 22_5

`PLAN.md`'s COVERAGE section had `24` tested and a `10`-base transcription gap; both were stale.
**TIER 1' is finished as a transcription task** — `57_1` was the last transcribable base
(`64d9316`, a paired presentation like `21_2`; 8 bases, 112 s, still dominated by `21_2`'s ~100 s).
`14_3` and `22_5` are NOT transcribable: we do not possess the object to compare, so they are
model-GENERATION items. `10_19` looks like a gap in the stored models but is not — its
`X0_10_19.m` re-derives the curve via `AllEquationsAboveCovers` instead of reading a model file.

⚠ **Counting trap: the obvious grep for the 43 bases returns 41.** Two rows write the label
without braces round `D` (`$X^6_0(17)$`, `$X^6_0(29)$`), so a pattern anchored on `X^{D}_0(N)`
drops exactly those two and yields a plausible 41. Cross-check on the equation cell instead:
`multirow{1}{*}{\text}` occurs 43 times. Also, `6_17`/`6_29` appear ONLY in CM-value captions
elsewhere — having a `tests/X0_6_17.m` does not imply a published equation — and `15_1` has a test
but is not a GY equation base at all.

### Two code changes

| commit | change | status |
|---|---|---|
| `d9b52d0` | `BorcherdsForms`: shift the oo-side basis by its own valuation, not the 0-side `n0` | fix landed, **NOT yet shown to unblock any base** |
| `36ac71e` | `EquationsCovers`: `Y2TWIST=1` prototype, decide the unpinned twist instead of dropping the cover | works, and **yields 0 new GY equations** |

**The vx fix (`d9b52d0`).** The "vx class" crash is Magma's own `assert vx ge 0`
(`GalFldFun.m:305`) reached from `AbsEltseq` on a deep Laurent pole (`93_1`: `q^-60`). Cause is a
wrong-object normalisation at `BorcherdsForms.m:771`: `ech_fs_oo` holds the **oo**-expansions of the
**ZERO**-side etas `ech_etas_0`, but the shift applied was `n0`, which comes back from
`WeaklyHolomorphicBasis(... : Zero, n0 := n0)` and bounds the 0-side, not the pole at oo. Shift by
`max(n0, -min valuation)` instead, and carry the same `n_oo` into BOTH places that must agree:
`coeffs_to_divisor_matrix(-n_oo, ...)` (the shift DEFINES the column↔exponent mapping) and
`min_m := Minimum(min_m, -(n_oo + k - 1))` (`relevant_ds` must stay a superset of
`relevant_ds_0_oo`).
⚠ **Both of those were learned by running it, not by reading it.** Missing the `min_m` one made
`95_1` clear the assert and then die at `:891` with `column index not in [1..37]`; and `n_oo` was
unassigned for even `D` (the block computing it is odd-`D` only), caught by the regression.
**Safety property:** where every oo-pole already fits within `n0` — every base that currently works
— the maximum IS `n0` and the change is a literal no-op.
✅ **SWEEP DONE, AND THE FIX IS EXONERATED — do not re-open this.** The 8-base sweep came back
**6 IDENTICAL / 2 DIFFERS** (`22_3`, `15_2`). ⚠ The `DIFFERS` does NOT falsify the no-op claim, and
the instruction that stood here ("if any base says DIFFERS the commit needs revisiting") would lead
you to exactly the wrong conclusion. The discriminating test is **fresh-vs-fresh**: regenerate at
`d9b52d0~1` and compare to regenerating at HEAD. Both came back `IDENTICAL`, so the fix changed
nothing; the two bases differ because their COMMITTED models are stale (see below). No-op verified
on 7 of 7 testable bases, including odd `D` (`51_1`).

**The `Y2TWIST` prototype (`36ac71e`).** `find_y2_scales` cannot always pin the y2-scale from
sparse CM data, so `EquationsOfCovers` force-defers the cover (issue #36, `1768517`) and back-fill
usually cannot recover it. But the twist is decidable by machinery INDEPENDENT of the
Borcherds/Schofer path that produced the equation — the Eichler-Selberg point count that
`ModelVerification.m` runs as check [4]. Env-gated, off by default, and it accepts only when
exactly one squarefree twist survives at 3+ good primes, so it never trades a deferral for a guess.
Ground truth: at `22_5` it recovers `W={1,2,5,10}` at `d=1` and reproduces the committed
polynomial coefficient-for-coefficient. **So `models_22_5.m` is regenerable from current code with
the twist VERIFIED rather than trusted.** Scope, honestly: 1 cover of 4 withheld, 0 new equations,
`W={1}` still empty; the other three failed their SOLVES, so they are under-determined like
`14_3` — a different problem that this does not touch.

### THREE COMMITTED MODELS DO NOT REGENERATE — and now there is a test for it

Found while validating the vx fix, not looked for. `models_22_5.m`, `models_22_3.m` and
`models_15_2.m` do not reproduce from current code, all for one reason: the unpinned-y2-scale
guard (`1768517`, 2026-08-24 19:20) POSTDATES all three, so regeneration withholds covers they
contain. ⚠ **They are NOT wrong** — all three pass `ModelChecks` and their Guo-Yang comparisons.
They are *unreproducible*, which is a different failure and one nothing in the suite could see:
`ModelChecks` and `GuoYangEquations` read STORED models and never run the pipeline, and the
`X0_D_N.m` tests run it for only ~25 bases, each needing hand-written cover/AL data.

**`tests/_offline/ModelRegen.m` (`afa0412`) closes that gap** — auto-discovering, no per-base
authoring, works for all 81 models; regenerates and checks each committed cover is still produced
and still ISOMORPHIC. The three are listed in `MR_KNOWN_DRIFT`, reported rather than asserted away.
Two traps it cost: matching must be a **multiset** match (the first draft passed `22_3` clean while
it had LOST a cover, because two committed entries matched the same survivor), and selection must
be an **env var** (`MODELREGEN_BASES`) because `run_tests.m` `eval`s test files and a `name:=value`
argument is invisible there — it silently runs the default list instead.

**What `Y2TWIST=1` restores** (measured against the committed files): `22_5` FULLY (3/3,
coefficient-for-coefficient), `15_2` FULLY (12/12 keys, one cover differing in presentation but
`IsIsomorphic`), `22_3` 13/14. ⇒ The residual gap is **GENUS 0 by construction** —
`select_y2_twist` skips `X`g lt 1` because `HyperellipticCurve` needs degree >= 3, and both `22_3`
losses are conics. **Extending twist selection to conics is the next concrete step**, with `22_3`
as its regression target. This revises the scope note above: `Y2TWIST` is a REPRODUCIBILITY fix for
the model corpus, not the single-cover curiosity the commit message describes. Still 0 new
Guo-Yang equations.

### `26_3`: the Mobius anomaly is an `s` <-> `s~` SWAP

At discs `-267` and `-708` Guo-Yang's `s` sits in **our `s~` row**; the other 12 of 14 are correct,
and `s + s~ = 1` holds at every disc. The exact `z -> z/(z-1)` is just how an `s -> 1-s` swap looks
after the checker's cross-ratio normalisation — the involution was the shadow, not the cause.
⚠ NOT a CM-point selection ambiguity (the old framing): both values are the same point, and each
disc appears exactly once in our table. **Root cause: `s + s~ = 1`, the relation used to pin the
pair, is SYMMETRIC under exchanging them**, so it cannot resolve the ordering; the signs are forced,
only the labelling is free. Deliberately NOT fixed — what pins the ordering at the other 12 discs is
unidentified, and a tie-break without that invariant is a guess. Memory: `26-3-hauptmodul-swap`.

### Corrections made this session — do not re-derive these

* **`22_5` is not "unreproducible".** I claimed a fresh run drops `[1,2,5,10]` because the target
  cover sets differ and `08ce5fa` came from a lost path. **False** — `{1,2,5,10}` (label 7584,
  g=1) is in `GetHyperellipticCandidates()` and `Xstar`CoveredBy` today. The cover is withheld
  **on purpose** by the y2 guard, which postdates the model file by 12 hours.
* **And the committed entry is CORRECT**, measured: `VerifyModelSet` passes it 24/24 and
  discriminates the twist — six twists `d = -1,2,-2,5,-5,11` all fail with 3-5 failures. So that
  guard is CONSERVATIVE, which is what motivated `Y2TWIST`. "Do not overwrite `models_22_5.m`"
  still stands; the reason changed.
* **The vx crash is at `BorcherdsForms.m:771`, not `ShimuraQuotients.m:1116`.** My first candidate
  was the unnormalised `denom` in `IsHyperelliptic`; the traceback never reaches it. `:615` is
  exonerated by its own `assert minval eq -Minimum(...)`.
* **`111_1`/`119_1` were attempted 2026-09-03, AFTER the 66x speedup (`04f1d7b`, 08-29).** So
  "re-run them now that the basis step is faster" is NOT free progress.
* **`cmsupply`'s `CMVERD OK` does not apply to the full curve.** It iterates `Xstar`CoveredBy`
  (`ShimuraQuotients.m:1526`), the immediate covers, whose genera top out at 1 (`14_3`) and 2
  (`22_5`) — while GY's published curves there are genus 3 and 5. `OK margin 0` means "adequate for
  the easy targets, zero slack", nothing about `W={1}`.

### Runs in flight when this was written

**On lovelace, and note there are now TWO checkouts there.** The long runs use
`~/shimura/ShimuraCurveALQuotients` (at `05471c8`, behind `main`); the new work uses a SEPARATE
clone `~/shimura/vxfix` pinned to `36ac71e`, deliberately, because `AttachSpec` loads packages
lazily and pulling under a 17-hour run could swap code mid-flight. **Do not `git pull` the first
one while those jobs are alive.**

* `34_11` — ~17 h, inside `AllEquationsAboveCovers` (4 ambiguous-sign points, 16 combinations),
  RSS plateaued ~22 GB. Past `M0MultiplierExact` and `ValuesAtCMPoints` entirely. This is the run
  that would give a SECOND base ever to produce models.
* `10_61`, `14_43` — ~15 h, still in the absolute-values phase. No gate failures anywhere.
* `93_1`, `95_1`, `159_1` — the vx bases, in `~/shimura/vxfix`, output to `~/shimura/vxout`.
  ⚠ Ran locally first; that was a mistake — the Mac has 48 GB and this class peaked at 40.6 GB on
  `119_1`. Use lovelace for these.
* the 8-base regeneration sweep, `~/shimura/regen/`.

⚠ **Clearing the vx assert is necessary, not sufficient** — these bases may still die downstream,
and `95_1`/`159_1` sharing `93_1`'s cause is inherited from memory, not measured.

## Update — 2026-09-04: `tier1-models` is RETIRED; `main` is the only code branch

`main` was fast-forwarded to `tier1-models` (`475e72b`) — a clean FF, `main` was a strict
ancestor 40 commits behind with nothing of its own. **`main` is no longer "code only": it now
carries `paper/` too**, so every description of the split below is historical. The open SHIP
question "merge `paper/`, or write down that the split is deliberate" is answered by merging.

What that merge carried, beyond the paper: the per-coset `tau` fix in `M0MultiplierExact`
(`34_11` goes from failing to passing, and its multipliers match the Prop 9.15 closed form 9/9 —
see `PLAN.md` REPAIR), `tests/KudlaYangLocal.m`, and 34 Guo-Yang CM-value tables as offline tests
under `tests/_offline/` (verified after the merge to emit zero CI targets).

**Then `tier1-models` was retired**, closing that question: deleted local and remote (it was a
strict ancestor of `main`, so `git branch -d` accepted it — nothing lost), and
`worktrees/mainport` removed as redundant. The layout is now just:

    .                    main  (this checkout)
    worktrees/campaign   m0-theta-campaign

⚠ **Everything below describing a `main` / `tier1-models` split is HISTORICAL.** `main` is not
"code only" any more; commit to it directly and do not recreate `tier1-models`.

Also confirmed 2026-09-04 while checking: the branch housekeeping this file and `PLAN.md` list as
open is already DONE — only `main`, `m0-theta-campaign` and `whbasis-speedup` exist, and the nine
retired branches (`non_optimal`, `odd_DN`, `pointlessconics` and the six `SCRATCH:` ones) are all
preserved as `archive/<name>` tags **on `origin`**, so `fix-15-2-find-signs`'s 3 never-pushed
commits are safe.

## Update — 2026-09-03 (evening)

Not a full rewrite; the record below (08-30) still stands for the model-pipeline arc. This adds
what happened on the `A_m` theorem (`MAIN LINE`) and two pieces of housekeeping.

**Worktree layout changed.** `worktrees/campaign` and `worktrees/mainport`, nested inside this
checkout — not sibling directories (`-campaign`, `-mainport`) as everywhere below still calls
them. See `CLAUDE.md` for the convention on new worktrees. Committed `827266b`.

**New CI test, `tests/KudlaYangLocal.m` (`a2e888c`).** Extends the existing Prop 5.4 check
(`mu = 0`) to Prop 5.5 (nonzero isotropic, `N`-only-supported coset): the level-prime local
Whittaker *value* there is the constant polynomial `1` for every `m` tested, every base tested
(1440 checks, including `N = 2`). This is a load-bearing negative result, not a tidy-up: it rules
out "insert KY Prop 5.4/5.5 into a level-`N` analogue of Theorem 8.1" as the source of `A_m`.

**Three derivation routes for `A_m` closed, all with reasons now on record (see `PLAN.md`, MAIN
LINE, for the full writeup):**
1. KY Prop 5.4/5.5 insertion — refuted above.
2. Schwagenscheidt's oldform relation (`eq:oldform`, `sec:ident`) — **the paper itself already
   tried this and documents why it fails** at weight 3/2 (needs analytic continuation, hits a
   non-holomorphic term; numerically witnessed as `-5` where the identity forces `0`). Missed on
   first read this session; re-read `sec:ident` before re-attempting anything like it.
3. The `s`-law / genus-theta closed form (`sec:slaw`) — this is the derivation *behind*
   `prop:closedcoef` (`-a_E(m)`), i.e. the already-refuted scalar route under a different name.

**The `rem:gauge` ambiguity is real, and provably resolvable — but the resolution needs data that
doesn't exist yet.** `-a_E(m)` and "Table A" disagree at 6 indices on `X_0^{15}(2)`
(`m=1,2,3,10,15,30`) while both fit the 9-form panel (rank 4 of 6, 2-dim kernel; `gauge152.py`).
The 158-monomial data already computed in `cusp7_15_2.out` (campaign branch) gives that same
6-index matrix **full rank 6** — this specific ambiguity is not one of `sec:exact`'s "50
unremovable" directions. But a naive per-monomial solve is invalid: `c_{eta*}(0)` sums over every
cusp class, and individual eta-monomials — unlike genuine Borcherds principal parts, protected by
[GY, Lemma 24] / `prop:nohalf` — can carry a nonzero constant term at *intermediate* cusps that an
`(A_m, B_j)`-only model never sees. Restricting to combinations with zero constant term at the
intermediate classes `cusp7_15_2.out` actually recorded (`g=3,4,5,12,15,20`, from its `PP` dump)
still gives rank 6, but the numeric solve against real `c_{eta*}(0)` remains inconsistent — because
that dump is missing 4 more intermediate classes (`g=2,6,10,30`), never captured by `cusp7.m`'s
first-encountered-per-class logic. **Next step: re-run a `cusp7.m`-style pass that guarantees every
intermediate class gets dumped, then redo the reachable-subspace solve.**

**Worktree/branch housekeeping, surveyed but not executed** (see `PLAN.md`, HOUSEKEEPING, for the
full list): the three stale remote-only branches (`non_optimal`, `odd_DN`, `pointlessconics`) are
still there, plus six more local branches missed by the 09-02 cleanup because they're genuinely
unmerged (mostly `SCRATCH:`-prefixed, predate 09-02) — `fix-15-2-find-signs` notably has 3 commits
that never even reached its own remote. One untracked stray file in `worktrees/campaign`
(`polymake/nmzsolve.err`, harmless).

    main               1c53865   speedup + zero-skip + hoist + IntegralSolution + Targets
                                 + slash-constant tolerance + pointless-conics guard
                                 + models_58_5.m + 3 new tests + cache (394 files)
    tier1-models       (merged)  carries the paper work; `origin/main` merged IN on 2026-09-02,
                                 so it now has the full code side too. The two had diverged
                                 27/27 on a clean split: paper/ on this side, all code on main.
    m0-theta-campaign  9059bb0   research branch: triage results, probes, predictors
    odd-d-zeroskip     b7067c3   MERGED to main; branch kept as the CI-green record
    odd-d-invariant-hoist 4c29d1e MERGED to main (afb80b2); correctness/clarity, ~2%
    intsol-optin       969fa85   MERGED to main; CI green
    fix-pointless-conics-empty 6071772  MERGED to main (133de9c); CI green
    whbasis-speedup    624b68e   MERGED (cherry-picked as 04f1d7b); branch kept likewise

### Where the model pipeline actually stands (2026-09-02)

**Four verified models exist** — `data/models/models_58_5.m`, `ModelChecks` 48/0. That is the
*only* base that has produced models. Everything else attempted since failed:

    34_11   gates 1-3 OK, then class-constancy (GENUINE, 43% of scale -- see below)
    74_5    gates 1-3 OK, capped at 5 h in the CM-value stage, no verdict
    74_3    no longer crashes (guard merged), but yields 0 keys
    10_61   slash-constant check, at the 1e-15 calibration
    14_43   slash-constant check, at the 1e-15 calibration

**⚠ Gate 3 is NOT fully solved.** The absolute→relative fix unblocked `58_5` decisively, but
`10_61` and `14_43` still fail at `1e-15` — the siblings' own calibration. **Do not make a third
tolerance change**: either their constants genuinely disagree (real mathematics) or something
else is wrong. A `GATE3B=1` measurement on `10_61` was running on lovelace when this was written;
its output lands at `~/shimura/models/10_61.gate3b.log` there.

**Gate 4 is a genuine violation** and is characterised: `vvdata/weyl-campaign/gate4/` on the
campaign branch. Deviation tracks the CLASS, pointing at the cusp-class partition.

**Both cheap-predictor routes are closed by measurement** for the 81 unclassified bases: route A
still costs a full `WeaklyHolomorphicBasis` (20-24 min and rising), route B's `k = 3/2` phase is
known wrong. `deficit.m` has been REPAIRED and validated (`38_5` -> 1 at every pole order).

Worktrees: `-campaign`, `-mainport` (main), `-spanprobe` (**THROWAWAY**, carries the live
instrumentation — the template for the next profiling pass). `-whspeed`, `-oddd` and
`-diagnostic` were removed on 08-30; see "Housekeeping" below for what was rescued first.

---

## What landed on 2026-08-30

**1. A CI failure, fixed.** `95bd502` ("PROTOTYPE (do not merge): integrality as an acceptance
criterion") was an ancestor of the campaign tip. It replaces the divisor solve's solution and
rejects triples, so `fs[-1]` becomes a *different* form and the reference comparisons in
`SchoferIsometry.m` (Guo–Yang Table 45) and `VectorValuedForm.m` (15_2 multiplier) fail.
Reverted in `badfe5d`; preserved as `vvdata/weyl-campaign/intsol-acceptance-criterion.patch`.
The prototype's *finding* stands (33_2 does go integral under it) — what was reverted is
shipping it as unconditional pipeline behaviour.

**2. "Failed to find all Borcherds forms" is the genuine BORCHERDS OBSTRUCTION.** Not a bug, not
a too-small space. At `38_5` the rank deficit is exactly 1 and **invariant** under deepening the
pole order (bump 0→8: rows 164→172, cols 36→38, rank 35→37); the annihilator φ is stable under
enlargement (hence a fixed modular object) and φ(target) = −22 ≠ 0. Bases that **succeed** have
deficit 0 at every key and find all forms on the first triple (`34_3`, `38_7`) — so the models
in CI never meet this because they have no obstruction space at all. **Not** a level threshold:
`38_7` is larger than `38_5` in every dimension and still surjective.
⇒ **Neither more divisor triples nor deeper poles can ever help an obstructed base.**

Untested escape hatch: φ(target) is **even** and gcd(φ) = 1, and a double cover depends on its
branch divisor only mod 2, so an even correction kills the pairing without changing the cover.
Unresolved: whether that introduces an unramified quadratic twist, and the exact-divisor
`assert` would need relaxing.

**3. A 66× speedup, merged to `main`.** Two changes in `WeaklyHolomorphicBasis`: select a
spanning row subset mod p before the echelon, and skip zero terms in the basis reconstruction
(the latter was the real bulk). At `38_5`: echelon 755 s → 1.5 s, basis 810 s → 12.3 s,
**end-to-end `ppint.m` 856 s → 30 s with the verdict unchanged.**

**4. Wave-4 triage: 11 of 18 previously-unmeasured bases recovered.**

    9 form-failure  10_43 10_47 142_3 22_23 46_11 74_7 82_5 86_5 94_5
    1 NONINTEGRAL   14_37        1 assertion  115_2
    7 still TIMEOUT 65_2 6_73 77_2 85_2 91_2 119_2 146_3

Backlog tally (wave 3 → now): form-failure 19 → **28**, NONINTEGRAL 20 → 21,
TIMEOUT 18 → **7**, assertion 4 → 5, INTEGRAL 4, CM-starved 1.
Results in `vvdata/weyl-campaign/triage-wave4/` on the campaign branch.

**5. The solution cache is unified on `main`** (`f90c441`): 43 new Normaliz solves from wave 4
extend the frontier from **M ≤ 1212 to M = 2260**, plus the two M = 1236 files that had been
committed to `m0-theta-campaign` only. Cache on main: 333 → 376 files. It *was* split across
branches — with the Normaliz backend an uncached level silently costs a full solve rather than
erroring, so a split cache produces confusing "why is this slow" sessions.

---

## 6. The odd-D branch, PROFILED — and one line was 93% of it

`BorcherdsForms.m:817` on `main` was the **last** `T[i][j]*` recombination in the file lacking
the `| T[i][j] ne 0` guard (others: ~437, ~480, ~618). It sits in the **odd-D-only** 0-side
block — exactly why the 66× speedup helped even D and left odd D untouched.

    stage             65_2    77_2    85_2
    oo basis           9.3     8.4     4.0
    0-cusp basis       3.6     2.4     4.0    <- REFUTED as the cause, <=3%
    CM points          0.10    0.11    0.09   <- REFUTED, negligible
    everything after  1787    1789    1792

`T` there is a pure SELECTION matrix (measured `nnz(T) = Nrows(T)`, one 1 per row), so the
unguarded sum did 79–426× more `EtaQuot` arithmetic than needed. **Fixed in `b7067c3`, now
MERGED to main as `619051a`, measured ~140×**: `etarecomb` 1.515 s → 0.0106 s per call;
`zside` 1539 s / 3 passes → 45.0 s / 4 passes. Tests green (incl. `15_2`, odd D).

**It does not make `65_2` / `85_2` complete.** The dominant cost is now
`basis_of_weakly_holomorphic_forms(... : Zero)` — real work, steep in pole order: 2.25 s @130,
26.45 @325, 72.85 @455, **556.12 @845**; `65_2`'s last m implies pole order 8450. Constant
factor removed, ceiling unmoved.

**~~The larger win, NOT yet done~~ — DONE, and REFUTED as a lever** (`odd-d-invariant-hoist`,
`4c29d1e`). The invariance is real: the 0-side block was recomputed 336× per `m_idx` pass at
`65_2` (= 8·7·6 triples, one key each before the break) and 210× at `85_2`, and it now runs
once. **But what repeats is cheap.** Instrumenting the *pre-hoist* code directly at `65_2`:

    pass 1   336 executions    1.75 s   (mean 0.0052)
    pass 2   336 executions    5.97 s   (mean 0.0178)
    pass 3   336 executions   10.02 s   (mean 0.0298)

≈ 17.7 s over three passes against an 1800 s cap, versus ≈ 0.9 s hoisted — **about 2%.**

**Why the "336× lever" claim was wrong, and the lesson.** That judgment was formed when the
block cost 1539 s / 3 passes — but *that* cost was the unguarded recombination line, and the
zero-skip removed it (140×), collapsing the block to 1.75 s/pass. The redundancy framing
outlived the fix that made it irrelevant, because nobody re-measured the block after changing
it. **A multiplier (336×) is only a lever when multiplied by something expensive; re-measure
the multiplicand after any fix that touches it.**

**The `T`-shadowing worry was vacuous.** The ∞-side `T` was never *read* — on odd D the 0-side
kernel overwrote it immediately, on even D nothing below touches `T` at all. It was dead code,
now deleted; the 0-side matrix is renamed `T_ker0` so the question cannot recur. One real trap
found while moving it: the hoisted lines must stay together, since `ech_etas_0` is sliced out of
`ech_etas_all_0` and then *replaced in place* by its own recombination — hoisting the
recombination without the slice recombines an already-recombined list on the second key, a
silent wrong answer rather than a crash.

⇒ **The odd-D constant factors are now exhausted. Everything left is the
`basis_of_weakly_holomorphic_forms(... : Zero)` ceiling above.** Do not spend more time here.

**Reclassify**: `133_2` is **not** a TIMEOUT — it fails an assertion at 168 s.

---

## Housekeeping — 2026-08-30 (later)

`odd-d-zeroskip` went CI-green (1h37m, 0 failures) and is **merged to main as `619051a`**.
Verified before the merge: all four `T[i][j]` recombination sites (437, 480, 618, 826) now
carry the zero-skip guard, and the merge brought one commit touching one file — the
`tier1-models` merge trap did not apply, because this branch was cut from `main`.

Three worktrees were retired. Everything single-copy in them was rescued to
`vvdata/weyl-campaign/` on the campaign branch first (`4b752d8`, `9059bb0`):

    bfprof.m  dsmall.m  diag_15_2.m         drivers that existed nowhere else
    bfprof-instrumentation.patch            BFPROF/BFINV timers  -- see note-probes.md
    valuesatcmpoints-characterization.patch the non-rationality probe        "
    MISSING_TARGETS.txt                     351 bases -- see note-missing-targets.md
    note-probes.md  note-missing-targets.md the caveats, which matter more than the code

**Two caveats worth carrying forward** (both in `note-probes.md`): the bfprof patch does **not**
apply to current main — its hunks mix the timers with a superseded inline prototype of the
speedup, so lift the timers by hand; and the characterization probe **cannot have run as
written**, since it patches the two-argument `ValuesAtCMPoints` at `SchoferFormula.m:1498`,
which has no `Xstar` in scope while the added lines reference `` Xstar`N ``. Any conclusion
attributed to that probe is unevidenced.

**`-whspeed`'s 61 uncommitted files were discarded, having been shown worthless**: 43 were
polymake solutions byte-identical to main's, and the 17 `data/curves_after_*.dat` were an
*incomplete* pipeline re-run — same 18379 records, but **strictly fewer** `IsHyp`/`IsSubhyp`
determinations at every stage (−190 UpdateByGenus, −340 UpdateCurves1, −95 UpdateCurves6).
`git diff HEAD origin/main -- data/` was empty, so main was never affected.
*Method note:* `grep -v TestInWhichProved` does **not** strip the attribution — the string sits
on a continuation line, which makes that diff look like ~58k lines of content change when it is
almost entirely attribution. Parse by splitting on `*])` and keying on `CurveID`.

Also harvested: two M = 532 Normaliz solutions left behind in `-spanprobe` (`9cc771e`). The
other two at that level were already tracked, so the cache was partial exactly there — the
silent-full-resolve mode `f90c441` set out to close. Cache on main: **378 files**.

---

## NEXT — in this order

**~~1. The invariant hoist~~ — DONE and merged-pending on `odd-d-invariant-hoist` (`4c29d1e`).
Worth ~2%, not the lever. See section 6. Nothing further to do on odd-D constant factors.**

**~~1. Decide the even-correction escape hatch~~ — MEASURED, and BLOCKED on open theory.**
Full account and tooling: `vvdata/weyl-campaign/even-correction/` on the campaign branch
(`e8d68f5`). Three results:

* **The precondition holds 28/28.** `φ(target)` is EVEN with `gcd(φ) = 1` at *every* obstructed
  base. Nothing is out of reach on parity grounds. The probe aborts at the first failing key
  (the deficit is invariant across triples), turning each base from a 900–1700 s exhaustive
  failure into one key — `38_5` reproduced exactly in **29.8 s instead of 860 s**.
* **CORRECTION to "the deficit is exactly 1"**, which was measured at `38_5` alone: `166_3`,
  `22_19` and `74_7` have a **2-dimensional** obstruction space and need a simultaneous
  2-condition solve, not a single shift. Both values are even in all three.
* **The correction is constructible but unusable.** A positive control at `34_3` (whose baseline
  reproduces the committed model exactly) builds forms with `div_f` exactly `ram + <disc,2>` at
  every key, then dies in `ValuesAtCMPoints`. Diagnostic: baseline **0** non-rational cells,
  perturbed **17**.

**Why it is blocked.** The mechanism is the `KNOWN DEFECT` at `SchoferFormula.m:589` — `Kappa0`
returns a zero log-`N` coefficient at firing discriminants where it should return `A_m`. **The
preprint does NOT supply `A_m`**: `prop:closedcoef` gives the *scalar* `a_E(m)` for all `m`, but
it reproduces only **1 of 13** measured `A_m`, and structurally so — `a_E` carries the embedding
support rule and vanishes exactly where `A_m` is nonzero (`15_2` m=2; `21_2` m=2,6,18). Per
[[b-eisenstein-coefficients-solved]] the relation is `A_r = -b^{η*}_0(r)/4` with `b` the
**vector-valued** coefficient at a **nonzero isotropic coset** (support `N | r`) — a different
object, and one that no product of local densities reproduces under any convention.

⇒ **The next theorem is: general `m` at a nonzero isotropic coset.** The preprint has `m = 0`
there (`prop:kappa0`) and all-`m` for the scalar (`prop:closedcoef`); `A_m` needs the
intersection. Until that exists this hatch cannot be finished, so **do not re-attempt it as an
implementation task** — and note the same defect is what the `coprime_to_level` filter
(`ShimuraQuotients.m:1420`, self-described as "a blunt instrument") already works around.

**1. The NONINTEGRAL class is mapped, and blocked at a THIRD gate.** Full account and tooling:
`vvdata/weyl-campaign/intsol/` on the campaign branch (`fe39373`, `c27707c`, `144b20f`); code on
branch `intsol-optin` (`6a6267c` the opt-in parameter, `4cdf1fb` the `Targets` threading), CI
clean. What was established:

* `IntegralSolution := false` makes the reverted August `intsol` finding usable — the prototype
  only ever failed because it shipped *unconditionally* and changed `fs[-1]` on working bases.
  It rescues **7 of 18** measured bases, so **≈39% of the class was a choice artifact** (which
  point `Solution` returned from `sol + Kernel`), not a divisor defect.
* **But none of the seven yields a model.** Six die of CM starvation, and **`cmsupply.m`
  predicted every one** at `ppint` cost — a 7/7 validation. **Run `cmsupply` FIRST on this
  class.**
* Genus-capping via `Targets` **clears** the CM gate (`58_5`, `74_5` ran 18 and 37 min instead of
  dying at it), but then `34_11`, `58_5` and `74_5` all fail the **`M0MultiplierExact`
  slash-constant two-point check** — three bases by two independent routes.

    integrality  →  CM supply  →  M0MultiplierExact slash-constant check

* **GATE 3 IS FIXED, and it was a miscalibrated tolerance** (`79d4e89` on main). The
  slash-constant check compared two evaluations with an **absolute** `1e-30` while the other four
  guards in `M0MultiplierExact` are relative — the one site the merged
  `m0exact-relative-tolerance` work never reached. Since `absdiff = reldiff * |k|` and `|k|`
  spans ten orders, it failed on LARGE constants whose agreement was unchanged.

⇒ **FOUR VERIFIED MODELS EXIST**: `data/models/models_58_5.m` (`3a50fa6`), the first output of
this entire line of work. `ModelChecks`: **48 checks, 0 failures**. `X_0(58,5)*` needed all three
fixes together — `IntegralSolution`, the `g ≤ 2` genus cap, and the tolerance — and fails without
any one. It is a **partial set by construction** (4 of 7 covers; the header says so).

**⇒ THE NEXT GATE IS 4, NOT 3.** `34_11` clears gates 1–3 and then fails
**class-constancy** (`dev 0.0186, scale 0.0430` — a 43% deviation). That is **NOT** another
tolerance: the check's own comment records roundoff at `1e-22` for the deepest known base and
states that *a genuine violation is O(scale)*. Loosening it would manufacture a multiplier wrong
by 43%. Treat as a real defect in the m=0 assembly at that base.

    integrality → CM supply → slash-constant (FIXED) → class-constancy (open, real)

**THE LESSON, and I got this wrong twice.** *Achievable precision is base-dependent*: at the same
`Prec := 80` the two evaluation points agree to **33 digits at `58_5` but only 18 at `34_11`**
(longer eta products, more accumulated rounding). I first made the guard relative but kept
`1e-30`, calibrated on `58_5` alone — that passed `58_5` and still blocked `34_11` and `74_5`,
which were otherwise ready. The siblings' `1e-15` is the right calibration and still catches what
the guard is for (a wrong constant differs at O(1), not in the 19th digit).
**Do not re-tighten a tolerance on the evidence of one base.**

Also worth carrying: losing 60 of 80 digits at `34_11` is real precision attrition — the first
place to look if a model from these bases ever appears suspect.

Two caveats from the earlier work still stand: `IntegralSolution` is **not monotone** (`69_2`
gets four orders of magnitude worse) so it must stay per-base opt-in; and `74_3`'s failure is
**unlocalised** because the driver truncated the error and had no verbosity — a bug in the
`Targets` threading is not excluded there.

**~~2. Re-run wave 4b (the 122 never-started bases)~~ — DONE.** Predictor sweep on **lovelace**
(256 cores, idle; `galois`/`verne` irrelevant, `legendre` busy and has no Normaliz, `lava` needs a
jump host). Full account: `vvdata/weyl-campaign/sweep122/` on the campaign branch (`a5d018f`).

     81 CAPPED-1h   21 OBSTRUCTED   13 VX-ASSERT   3 ASSERT
      2 NONINTEGRAL  1 CM-STARVED    1 INTEGRAL     = 122

* **TWO runnable candidates of 122**: `10_61` (INTEGRAL, CM OK margin 0) and `14_43`
  (NONINTEGRAL — the *fixable* gate — CM OK margin 0). Both at margin 0, the position `34_11`
  was in when it cleared gates 1–3.
* **⚠ THE OBSTRUCTED CLASS IS 49, NOT 28** — 21 more bases fail "Failed to find all Borcherds
  forms", which neither more triples nor deeper poles can help. **This raises the priority of the
  theory item (5) substantially: it is now worth 49 bases.**
* **`ppint`/`cmsupply` are NOT cheap predictors on large bases.** 81 of 122 gave no verdict in a
  full hour — `ppint` must build Borcherds forms before it can speak. My earlier advice "run
  `cmsupply` first, it is `ppint`-cost" holds only for bases earlier waves already reached. A
  600 s cap was strictly worse than useless: **the cap was bounding the measurement itself.**
  The remaining 81 need a genuinely cheaper predictor, not more wall-clock.

**Remote-run notes** (lovelace): Magma 2.29-9, Normaliz 3.10.2 at `/usr/bin/normaliz` — verified
to produce lattice points identical to local 3.11.1. Clone at `~/shimura/ShimuraCurveALQuotients`.
`pkill -x magma` matches NOTHING there (the binary is `magma.exe` behind a wrapper) — and
verifying a kill with the same pattern used to kill reports false success, which caused a
double-launch here. **Verify with a different pattern.**
**Pre-solve the cache first**: of the 351 bases in `MISSING_TARGETS.txt`, 328 sit inside the
committed M ≤ 2260 frontier and 23 do not — but those 23 share only **11 distinct M** (the cache
key), so it is ~22–33 solves, a bounded batch to run *ahead* of the wave rather than a silent
per-base tax inside it. That cohort is also the high-genus tail (g up to 17, CM demand
`max(2g+5)` = 39), so run `cmsupply.m` over it first — see `note-missing-targets.md`.

**3. Re-run the 7 remaining TIMEOUT bases.** Low value, and expect it to *confirm* rather than
clear them: both odd-D constant-factor fixes are in and the ceiling is untouched.
`basis_of_weakly_holomorphic_forms(... : Zero)` is steep in pole order (556 s @ 845; `65_2`'s
last m implies 8450) and `77_2` is structurally out of reach.

**4. Route B's k = 3/2 phase.**

**5. The theory item, if the paper is the priority:** state and prove the general-`m` analogue of
`prop:kappa0` — the vector-valued weight-3/2 Eisenstein coefficient at a nonzero isotropic coset.
It is the one object standing between the obstructed class (**49 bases** — see item 2) and a model, and the exact
values are already known at `15_2`, `6_5`, `10_3`, `21_2` as a regression set. **No longer low
priority: at 49 bases this is the largest single blocker in the backlog.**

---

## TRAPS — recorded so they are not repeated

* **Wave 4b's numbers are confounded; discard them.** The cap was halved (2400 → 1200 s) *and*
  concurrency raised to 8 streams (load 30 on 14 cores), so its near-total timeout rate measures
  the scheduling as much as the bases.
* **REFUTED: "the odd-D eta-quotient explosion is the blocker."** Proposed from a parity pattern
  without checking the mechanism. Measured: `nsol` ≈ 12k on odd D exactly as on even (`133_2`
  11964, `65_2` 14346, `38_5` 12784) and `t_ip` is *instant*. No explosion.
* **`deficit.m` is EVEN-D ONLY** (guard now in the file). For odd D it omits the 0-side block
  joined into `coeffs_trunc`, so deficits are **overestimates**. Tell-tale: the value *drifts*
  with pole order instead of staying invariant (`65_2`: 5, 6, 6, 9) where a real one is constant
  (`38_5`: 1 everywhere). The three validated results are all even D and stand.
* **MERGE TRAP.** `whbasis-speedup` branched off `tier1-models`, so merging it into `main` would
  have brought **26 commits** including the entire unmerged paper rewrite
  (`level-prime-kappa.tex` +1244, the PDF, `gtsweep.m`). Cherry-pick instead, and always run
  `git log origin/main..origin/<branch>` before merging anything off `tier1-models`.
* `git -C <repo> worktree add <relative-path>` resolves the path against the **repo**, not your
  cwd — it will silently create a worktree *inside* the repo. Use absolute paths.
* `magma | head` / `| tail` can hang; redirect to a file instead.
* **The `PROBESPAN` printf in `-spanprobe` ran two `Rank()` calls per key.** Any timing taken in
  that worktree before 2026-08-30 is inflated. Now gated behind `PROBE_SPAN=1` (default off).
* **A multiplier is only a lever when the multiplicand is expensive.** The "336x redundancy"
  claim was formed when the 0-side block cost 1539 s / 3 passes, survived the zero-skip that
  collapsed it to 1.75 s/pass, and was still being quoted as "the biggest remaining lever" a
  session later. Re-measure the multiplicand after any fix that touches it.
* **Measure the thing you changed, not the whole pipeline.** The 817 fix was first tested
  end-to-end with a 2400 s cap: both bases timed out before *and* after, so the test could not
  have detected the 140× win it actually produced. For partial speedups use the `BFPROF`
  per-stage timers in `-spanprobe` against the recorded baselines.
* `ppint.m`'s first `printf` fires only *after* `BorcherdsForms` returns, so an empty log tells
  you nothing about where a run is — instrument if you need progress.

## Tools (campaign branch, `vvdata/weyl-campaign/`)

    spanprobe.m  deficit.m  matrank.m  dsize.m    route A (measured deficit)
    weildim.m  weildim2.m  dsmall.m              route B (Weil rep) + its #disc_grp table
    bfprof.m                                     per-stage odd-D profiler — USE THIS to
                                                 measure the invariant hoist
    diag_15_2.m                                  non-rationality characterization driver
    span-obstruction-probe.patch                 instrumentation — THROWAWAY WORKTREE ONLY
    bfprof-instrumentation.patch                 BFPROF/BFINV timers — DOES NOT APPLY to main
    valuesatcmpoints-characterization.patch      probe — CANNOT HAVE RUN as written
    note-probes.md                               the caveats on both patches. Read first.
    MISSING_TARGETS.txt  note-missing-targets.md 351-base target list + cache-frontier analysis
    retriage.sh  wave4_*.txt                     triage driver + stream lists
    ppint.m  cmsupply.m  genmodels.m  backlog.m  earlier triage tooling

`matrank.m` records a refuted shortcut: `coeffs_to_divisor_matrix` has **full column rank**, so
the deficiency is pure Borcherds duality, not a property of the divisor matrix.
`weildim2.m`'s O(d) trace formulas are cross-checked against explicit Weil matrices at `6_1`
(all ten traces agree), but its k = 3/2 phase is wrong by very nearly **−d/6** — **do not tune
that constant to fit**; get the half-integral convention right, then check it against the
measured deficits (`38_5` → 1, `38_7` → 0, `34_3` → 0).
