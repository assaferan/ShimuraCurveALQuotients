# `IntegralSolution`: the NONINTEGRAL class, measured

Result of making the 2026-08-29 `intsol` finding usable. Code landed on branch `intsol-optin`
(`6a6267c`, off `main`); the reverted prototype it descends from is
`../intsol-acceptance-criterion.patch`.

## What was actually wrong before

The finding was never in doubt — what blocked it was that it shipped **unconditionally**. It
replaces the solve's solution, so `fs[-1]` becomes a different form and the reference
comparisons in `tests/SchoferIsometry.m` (Guo–Yang Table 45) and `tests/VectorValuedForm.m`
(the 15_2 multiplier) fail, on bases that are already fine. The prototype's own comment says why
it was unconditional: **`assigned` cannot see a top-level variable from inside a package
intrinsic**, so a flag would never have fired. A named **parameter** threads correctly. That is
the entire fix.

Default `false` ⇒ inert, bit-for-bit unchanged. All three tests that call `BorcherdsForms` pass
(`SchoferIsometry` 25.5 s, `M0MultiplierExact` 18.6 s, `VectorValuedForm` 97.3 s). Threaded
through `EquationsOfCovers` and both `AllEquationsAboveCovers` forms so model generation can opt
in (`INTSOL=1` in `genmodels.m`).

## The measurement (`intsol.m`, both passes per base)

    RESCUED (integral, den 1)   33_2  34_11  46_3  57_2  58_5  74_3  74_5
    better, not cleared         14_37 (6->2)   39_2 (65536->32)   82_3 (3->2)
    unchanged                   22_17  46_5  51_2  62_3  62_7  86_3
    WORSE                       69_2 (64 -> 1048576)
    rejected every triple       26_7
    no verdict                  111_2 (capped at 2h32m)   87_2 (still running)

**≈39% of the class is a choice artifact** — the non-integrality was which point `Solution`
returned from `sol + Kernel`, not a property of the divisor.

## Two independent confirmations of earlier work

* **`74_3` is RESCUED.** [[nonintegrality-origin-echelonform]] used that very base to conclude
  the defect is in the SOLUTION the divisor solve picks, not the basis. Fixing the solution
  rescues it — exactly as predicted.
* **`39_2` improves but does NOT clear** (65536 → 32). [[nonintegral-forms-root-cause]] records
  `39_2` as a **basis-level** failure `intsol` cannot touch. Also exactly right.

Both were written from different evidence, months apart, and both hold.

## THE CAVEAT: it is not monotone — keep it opt-in

**`69_2` gets four orders of magnitude WORSE** (denominator 64 → 1048576). An integral *solution
vector* against a *non-integral basis* can produce worse principal parts than the arbitrary
rational solution did — the two objects are not the same, and `EchelonForm` is what makes the
pool non-integral ([[nonintegrality-origin-echelonform]], denominators ~3e28; `111_2`'s default
pass reports ~3.2e30, which is that pathology rather than ordinary arithmetic).

So **do not promote this to an automatic fallback.** Per-base opt-in is the right granularity.

`26_7` shows the other edge: rejecting non-integral triples can leave *nothing*, moving a base
from NONINTEGRAL to form-failure rather than to a model. Honest, but not obviously an
improvement.

## Necessary, not sufficient

Integral principal parts do not make a model. The seven still have to clear CM supply and the
`y2`-scale assembly — and `58_5` sits right beside the documented CM-starved `58_3`
([[cm-supply-second-triage-axis]]). Anything that does produce a model must go through
`ModelVerification` (genus, Weil-polynomial divisibility, trace-formula point counts), none of
which touches the Borcherds/Schofer path that produced it.

## Method note

The concurrency limit in the first sweep runner (`zsh`, `jobs -r | wc -l`) **silently did not
work** in a non-interactive script — all 16 jobs launched at once (load 33 on 14 cores). Memory
was fine and the verdicts are unaffected, because this sweep measures integrality and
denominators rather than durations — but **the timing column in these logs is confounded and
must not be quoted.** Use `xargs -P N` for a limit that actually holds.

## Files

    intsol.m           the two-pass driver (default vs IntegralSolution)
    sweep_<base>.log   all runs

---

# FOLLOW-UP: integrality was NOT the binding constraint — CM supply is

All seven rescued bases were run end to end with `IntegralSolution := true`
(`cmsupply/genmodels_<base>.log`). **None produced a model.**

    33_2 46_3 57_2 58_5 74_3 74_5   Could not find enough points, sorry!
    34_11                           M0MultiplierExact: slash constant failed its two-point check

`cmsupply.m` predicted this exactly, at `ppint` cost rather than full-pipeline cost:

    base    demand  margin  verdict          genera
    34_11     13       0    OK               0,1,1,3,4,4,4
    33_2       7      -1    SHORT            0,0,1,1,1,1,1
    57_2       9      -3    SHORT            0,1,1,1,2,2,2
    58_5      13      -3    SHORT            0,1,2,2,3,3,4
    46_3       9      -4    SHORT            0,1,1,1,1,1,2
    74_3      13      -4    SHORT            0,1,1,2,2,3,4
    74_5      15      -4    SHORT            1,1,2,3,3,4,5

Six SHORT, six died of CM starvation; the one called OK died of something else entirely.
**A clean validation of the predictor** — and it means this class can be triaged at `ppint` cost
from here on, instead of a full pipeline run per base.

**So `IntegralSolution` moves bases past a real gate but not to a model.** The finding stands and
the code is worth having, but the gate behind it is CM supply. Anyone continuing this class
should run `cmsupply.m` FIRST and only then think about integrality.

## The `Targets` lever, quantified

Demand is `max(2g+5)` over the RETAINED covers, so one high-genus sibling inflates it for
everyone. Capping the genus lowers demand while keeping most covers:

    74_3  g<=2  demand 13 ->  9  (supply  9)  keeps 5/7 covers
    74_5  g<=3  demand 15 -> 11  (supply 11)  keeps 5/7 covers
    58_5  g<=2  demand 13 ->  9  (supply 10)  keeps 4/7 covers
    33_2  g<=0  demand  7 ->  5  (supply  6)  keeps 2/7   -- not worth it
    46_3  g<=0  demand  9 ->  5  (supply  5)  keeps 1/7   -- not worth it
    57_2  g<=0  demand  9 ->  5  (supply  6)  keeps 1/7   -- not worth it

Only the first three retain enough covers to be worth running. Note two of them sit at **margin
0**, so they succeed only if `cmsupply`'s accounting is exact rather than approximately right —
a sharper test of the predictor than any of its SHORT calls.

Driver: `gencapped.m` (`GCAP` = genus cap; `Targets` restricted to covers with `g <= GCAP`;
`IntegralSolution := true`). Requires the `Targets` threading through
`AllEquationsAboveCovers` added in `4cdf1fb` — it must reach the `num_vals` computation, not
just the form search, or the cap changes nothing.

**Results are PARTIAL by construction**: 5/7, 5/7 and 4/7 covers. The high-genus covers are
deliberately abandoned and still need more CM points, not a different method.

## Method note (third occurrence tonight)

`genmodels.m` and the first `gencapped.m` did not call `SetVerbose`, so the `[n/6]` stage markers
in `EquationsCovers.m` (vprintf level 1) never printed and the logs stayed empty until the run
ended — making "is it progressing or stuck?" unanswerable. `gencapped.m` now sets verbosity 1 by
default (`VERB:=0` to mute). **Do not launch a long run without progress output.**

---

# CAPPED RUNS: the cap works, and a THIRD gate appears

`gencapped.m` on the three bases the `Targets` arithmetic said were worth running
(`capped/gencapped_<base>.log`). **No models.**

    58_5  g<=2  FAIL 1091 s  M0MultiplierExact: slash constant failed its two-point check
    74_5  g<=3  FAIL 2205 s  M0MultiplierExact: slash constant failed its two-point check
    74_3  g<=2  FAIL 4229 s  Maximum: Argument 1 is not non-empty   <- UNLOCALISED, see below

## The cap does what it was designed to do

`58_5` and `74_5` did **not** die of "Could not find enough points" — the failure that killed
them before. They cleared that gate and ran 18 and 37 minutes. The demand reduction worked
exactly as the arithmetic predicted, so genus-capping is validated as a technique even though it
did not produce models here.

## A third gate, reached by two different routes

`34_11`, `58_5` and `74_5` fail identically. `34_11` got there by having enough CM points
natively (`cmsupply` OK, margin 0); the other two by having their demand lowered. Three bases and
two routes make this a real gate rather than a quirk:

    integrality  ->  CM supply  ->  M0MultiplierExact slash-constant two-point check

That check is now the binding constraint for this class, and nothing in the integrality or
CM-supply direction will move it.

## `74_3` is UNLOCALISED, and that is an instrumentation failure of mine

The error `Maximum: Argument 1 is not non-empty` carries no location, because:

* `gencapped.m`'s `try`/`catch` truncated `e`Object` to 120 characters, discarding the traceback;
* verbosity was added to `gencapped.m` only AFTER these three had launched, so there are no
  `[n/6]` stage markers either.

So **I cannot rule out that this is a bug in the `Targets` threading** (`4cdf1fb`) rather than a
property of the base. What is ruled out: the `Maximum` at `EquationsCovers.m:734`, which the
`require not IsEmpty(aeac_keys)` guard protects. The leading candidate is
**`EquationsCovers.m:524`**, `Maximum(all_keys)` inside `EquationsAbovePointlessConics` —
assembly code that may assume the full cover set and would see an empty list once `Targets`
restricts it. **That is a hypothesis, not a finding.** Re-running with verbosity and an untruncated
error costs ~70 min and would not change the class verdict, since the class is blocked at the
`M0MultiplierExact` gate regardless.

Fixed for next time: `gencapped.m` now sets verbosity 1 by default. **The catch should also print
the full `e`Object` rather than a truncation** — a 120-character error is worth almost nothing
when the location is what you need.

## Net result of this line of work

**No new models.** What was gained instead:

* `IntegralSolution` as a correct, opt-in, CI-clean capability (7 of 18 bases past gate 1);
* `cmsupply.m` validated as an exact predictor of gate 2, at `ppint` cost (7 of 7);
* `Targets`/genus-capping threaded end to end and demonstrated to clear gate 2;
* gate 3 identified and named.

The failure chain for this class is now mapped, with a cheap predictor at gate 2. Anyone
resuming should start at the `M0MultiplierExact` slash-constant check, not at integrality.
